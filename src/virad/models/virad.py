"""VIRAD: multi-modal Viroformer mapping fθ(S, M, C, H)."""

from __future__ import annotations

from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

from virad.config import VIRADConfig
from virad.models.embeddings import ConservationEmbedding, HostEmbedding, SequenceEmbedding, SequenceEmbeddingWithMask
from virad.models.fusion import MultiModalFusion
from virad.models.gnn import StructureCNN2d, StructureGNN
from virad.models.heads import MotifHeads
from virad.models.losses import UncertaintyWeightedLoss, bce_logits, denoise_mlp_loss, dice_bce_loss, info_nce, kl_invariance
from virad.models.structure_module import StructureDecoder, StructureGate
from virad.models.viroformer import ViroformerEncoder
from virad.models.coevolution import JointEmbedding


class StructuralDenoiseDecoder(nn.Module):
    def __init__(self, cfg: VIRADConfig) -> None:
        super().__init__()
        self.fc1 = nn.Linear(cfg.d_model * 2, 256)
        self.fc2 = nn.Linear(256, 1)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        bsz, length, dim = z.shape
        out = z.new_empty(bsz, length, length)
        chunk = 32 if length > 32 else length
        for start in range(0, length, chunk):
            end = min(length, start + chunk)
            zi = z[:, start:end, :].unsqueeze(2).expand(bsz, end - start, length, dim)
            zj = z.unsqueeze(1).expand(bsz, end - start, length, dim)
            cat = torch.cat([zi, zj], dim=-1)
            out[:, start:end, :] = self.fc2(F.relu(self.fc1(cat))).squeeze(-1)
        return out


class VIRAD(nn.Module):
    def __init__(self, cfg: Optional[VIRADConfig] = None, use_cnn_structure: bool = False) -> None:
        super().__init__()
        self.cfg = cfg or VIRADConfig()
        self.seq_emb = SequenceEmbedding(self.cfg)
        self.seq_emb_mask = SequenceEmbeddingWithMask(self.cfg)
        self.cons_emb = ConservationEmbedding(self.cfg)
        self.host_emb = HostEmbedding(self.cfg)
        self.struct_gnn = StructureGNN(self.cfg)
        self.struct_cnn = StructureCNN2d(self.cfg) if use_cnn_structure else None
        self.use_cnn_structure = use_cnn_structure
        self.fusion = MultiModalFusion(self.cfg)
        self.encoder = ViroformerEncoder(self.cfg)
        self.struct_decoder = StructureDecoder(self.cfg)
        self.struct_gate = StructureGate(self.cfg)
        self.heads = MotifHeads(self.cfg)
        self.joint = JointEmbedding(self.cfg)
        self.uncertainty = UncertaintyWeightedLoss(6)
        self.denoise_dec = StructuralDenoiseDecoder(self.cfg)
        self.apply(self._init_weights)

    @staticmethod
    def _init_weights(module: nn.Module) -> None:
        if isinstance(module, nn.Linear):
            nn.init.xavier_uniform_(module.weight)
            if module.bias is not None:
                nn.init.zeros_(module.bias)

    def embed_modalities(
        self,
        sequence: torch.Tensor,
        conservation: torch.Tensor,
        structure: torch.Tensor,
        host: torch.Tensor,
        motif_mask: Optional[torch.Tensor] = None,
        zero_struct: bool = False,
        zero_cons: bool = False,
        zero_host: bool = False,
    ) -> dict:
        if motif_mask is not None:
            e_seq = self.seq_emb_mask(sequence, motif_mask)
        else:
            e_seq = self.seq_emb(sequence)
        e_cons = self.cons_emb(conservation)
        if zero_cons:
            e_cons = torch.zeros_like(e_cons)
        if self.use_cnn_structure and self.struct_cnn is not None:
            e_struct = self.struct_cnn(structure)
        else:
            e_struct = self.struct_gnn(e_seq, structure)
        if zero_struct:
            e_struct = torch.zeros_like(e_struct)
        e_host = self.host_emb(host, sequence.size(-1))
        if zero_host:
            e_host = torch.zeros_like(e_host)
        e_fused = self.fusion(
            e_seq,
            e_cons,
            e_struct,
            e_host,
            disable_struct_cross=zero_struct,
        )
        return {"e_seq": e_seq, "e_cons": e_cons, "e_struct": e_struct, "e_host": e_host, "e_fused": e_fused}

    def refine_structure(
        self,
        e_fused: torch.Tensor,
        conservation: torch.Tensor,
        pair_probs: torch.Tensor,
        e_host: torch.Tensor,
        attn_mask: Optional[torch.Tensor],
        disable_host_cross: bool = False,
        n_refine: Optional[int] = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        p = pair_probs
        z = e_fused
        n_refine = self.cfg.n_refine if n_refine is None else n_refine
        for _ in range(n_refine):
            z = self.encoder(
                z,
                conservation,
                p,
                e_host,
                attn_mask=attn_mask,
                disable_host_cross=disable_host_cross,
            )
            p_next = self.struct_decoder(z)
            delta = torch.norm(p_next - p, p="fro", dim=(-2, -1)).mean()
            p = p_next
            if float(delta.detach()) < self.cfg.refine_frobenius_tol:
                break
        g = self.struct_gate(e_fused, p)
        z_final = g * e_fused + (1.0 - g) * z
        return z_final, p

    def forward(
        self,
        sequence: torch.Tensor,
        conservation: torch.Tensor,
        structure: torch.Tensor,
        host: torch.Tensor,
        attn_mask: Optional[torch.Tensor] = None,
        motif_mask: Optional[torch.Tensor] = None,
        zero_struct: bool = False,
        zero_cons: bool = False,
        zero_host: bool = False,
        n_refine: Optional[int] = None,
        return_attn: bool = False,
    ) -> dict:
        emb = self.embed_modalities(
            sequence,
            conservation,
            structure,
            host,
            motif_mask=motif_mask,
            zero_struct=zero_struct,
            zero_cons=zero_cons,
            zero_host=zero_host,
        )
        cons_in = torch.zeros_like(conservation) if zero_cons else conservation
        struct_in = torch.zeros_like(structure) if zero_struct else structure
        z_final, p_final = self.refine_structure(
            emb["e_fused"],
            cons_in,
            struct_in,
            emb["e_host"],
            attn_mask,
            disable_host_cross=zero_host,
            n_refine=n_refine,
        )
        outs = self.heads(z_final)
        outs.update(emb)
        outs["z"] = z_final
        outs["p_refined"] = p_final
        if return_attn:
            _, attns = self.encoder(
                emb["e_fused"],
                conservation,
                p_final,
                emb["e_host"],
                attn_mask=attn_mask,
                disable_host_cross=zero_host,
                return_attn=True,
            )
            outs["attn"] = attns
        return outs

    def pooled_motif(self, z_motif: torch.Tensor, attn_mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        if attn_mask is None:
            return z_motif.mean(dim=1)
        w = attn_mask.unsqueeze(-1).clamp_min(0)
        return (z_motif * w).sum(dim=1) / w.sum(dim=1).clamp_min(1.0)

    def multitask_losses(
        self,
        outs: dict,
        batch: dict,
        positive_pairs: Optional[torch.Tensor] = None,
        z_mut: Optional[torch.Tensor] = None,
    ) -> dict:
        mask = batch.get("attention_mask")
        cons_target = (batch["conservation"] > self.cfg.cons_label_threshold).float()
        l_cons = bce_logits(outs["cons_logit"], cons_target, mask)
        l_motif = dice_bce_loss(outs["bound_logit"], batch["motif_mask"], mask)
        acc_mask = batch.get("shape_mask")
        acc_present = acc_mask is not None and bool(acc_mask.detach().any().item())
        if acc_present:
            l_acc = ((outs["accessibility"] - batch["accessibility"]) ** 2 * acc_mask).sum() / acc_mask.sum()
        else:
            l_acc = outs["accessibility"].new_zeros(())
        host_mask = batch.get("host_mask")
        host_present = host_mask is not None and bool(host_mask.detach().any().item())
        if host_present:
            l_host = bce_logits(outs["host_logit"], batch["host_label"], host_mask)
        else:
            l_host = outs["host_logit"].new_zeros(())
        z_pool = self.pooled_motif(outs["z_motif"], mask)
        contrast_present = positive_pairs is not None and int(positive_pairs.numel()) > 0
        if contrast_present:
            l_contrast = info_nce(z_pool, positive_pairs, self.cfg.contrastive_tau)
        else:
            l_contrast = z_pool.new_zeros(())
        mlm_mask = batch.get("mlm_mask")
        mlm_present = mlm_mask is not None and bool(mlm_mask.detach().any().item())
        if mlm_present:
            logits = outs["mlm_logit"]
            target = batch["mlm_target"]
            l_mlm = F.cross_entropy(logits[mlm_mask.bool()], target[mlm_mask.bool()])
        else:
            l_mlm = z_pool.new_zeros(())
        denoise_present = "p_clean" in batch
        if denoise_present:
            pred_p = self.denoise_dec(outs["z"])
            l_denoise = denoise_mlp_loss(pred_p, batch["p_clean"])
        else:
            l_denoise = z_pool.new_zeros(())
        l_aux = l_mlm + l_denoise
        losses = [l_cons, l_motif, l_acc, l_host, l_contrast, l_aux]
        present = [True, True, acc_present, host_present, contrast_present, mlm_present or denoise_present]
        total = self.uncertainty(losses, present)
        if z_mut is not None:
            total = total + self.cfg.lambda_kl * kl_invariance(outs["z"], z_mut)
        type_loss = F.cross_entropy(
            outs["type_logit"].reshape(-1, outs["type_logit"].size(-1)),
            batch["motif_type"].reshape(-1),
            ignore_index=-1,
            label_smoothing=self.cfg.label_smoothing,
        )
        total = total + type_loss
        presence_loss = bce_logits(outs["presence_logit"], batch["motif_mask"], mask)
        total = total + presence_loss
        return {
            "total": total,
            "cons": l_cons,
            "motif": l_motif,
            "acc": l_acc,
            "host": l_host,
            "contrast": l_contrast,
            "aux": l_aux,
            "type": type_loss,
            "presence": presence_loss,
        }
