"""Viroformer encoder blocks (Framework §VI)."""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from virad.config import VIRADConfig
from virad.models.attention import BiasedSelfAttention


class ViroformerBlock(nn.Module):
    def __init__(self, cfg: VIRADConfig) -> None:
        super().__init__()
        self.self_attn = BiasedSelfAttention(cfg)
        self.norm_attn = nn.LayerNorm(cfg.d_model)
        self.cross = nn.MultiheadAttention(
            embed_dim=cfg.d_model,
            num_heads=cfg.host_cross_heads,
            dropout=cfg.dropout,
            batch_first=True,
        )
        self.norm_cross = nn.LayerNorm(cfg.d_model)
        self.ff1 = nn.Linear(cfg.d_model, cfg.d_ffn)
        self.ff2 = nn.Linear(cfg.d_ffn, cfg.d_model)
        self.drop1 = nn.Dropout(cfg.dropout)
        self.drop2 = nn.Dropout(cfg.dropout)
        self.norm_ffn = nn.LayerNorm(cfg.d_model)
        self.stochastic_depth = cfg.stochastic_depth

    def forward(
        self,
        x: torch.Tensor,
        conservation: torch.Tensor,
        pair_probs: torch.Tensor,
        e_host: torch.Tensor,
        attn_mask: torch.Tensor | None = None,
        disable_host_cross: bool = False,
        return_attn: bool = False,
    ):
        if self.training and torch.rand(1, device=x.device).item() < self.stochastic_depth:
            if return_attn:
                return x, None
            return x
        attn_out = self.self_attn(
            self.norm_attn(x),
            conservation,
            pair_probs,
            attn_mask=attn_mask,
            return_weights=return_attn,
        )
        weights = None
        if return_attn:
            delta, weights = attn_out
        else:
            delta = attn_out
        z_attn = x + delta
        key_padding = None
        if attn_mask is not None:
            key_padding = attn_mask <= 0
        if disable_host_cross:
            z_cross = torch.zeros_like(z_attn)
        else:
            # PyTorch MultiheadAttention was the remaining BF16-sensitive path
            # after custom SDPA had been ruled out during production recovery.
            # Execute host cross-attention in FP32 and cast its output back so
            # the encoder trunk and FFN remain under BF16 autocast.
            if z_attn.device.type == "cuda":
                with torch.autocast(device_type="cuda", enabled=False):
                    z_cross, _ = self.cross(
                        z_attn.float(),
                        e_host.float(),
                        e_host.float(),
                        key_padding_mask=key_padding,
                        need_weights=False,
                    )
                z_cross = z_cross.to(dtype=z_attn.dtype)
            else:
                z_cross, _ = self.cross(
                    z_attn,
                    e_host,
                    e_host,
                    key_padding_mask=key_padding,
                    need_weights=False,
                )
        z_post = self.norm_cross(z_attn + z_cross)
        hidden = self.drop1(self.ff1(z_post))
        ffn = self.drop2(self.ff2(F.gelu(hidden)))
        out = self.norm_ffn(z_post + ffn)
        if return_attn:
            return out, weights
        return out


class ViroformerEncoder(nn.Module):
    def __init__(self, cfg: VIRADConfig) -> None:
        super().__init__()
        self.blocks = nn.ModuleList(ViroformerBlock(cfg) for _ in range(cfg.n_block))

    def forward(
        self,
        e_fused: torch.Tensor,
        conservation: torch.Tensor,
        pair_probs: torch.Tensor,
        e_host: torch.Tensor,
        attn_mask: torch.Tensor | None = None,
        disable_host_cross: bool = False,
        return_attn: bool = False,
    ):
        x = e_fused
        attns = []
        for block in self.blocks:
            if return_attn:
                x, w = block(
                    x,
                    conservation,
                    pair_probs,
                    e_host,
                    attn_mask=attn_mask,
                    disable_host_cross=disable_host_cross,
                    return_attn=True,
                )
                attns.append(w)
            else:
                x = block(
                    x,
                    conservation,
                    pair_probs,
                    e_host,
                    attn_mask=attn_mask,
                    disable_host_cross=disable_host_cross,
                )
        if return_attn:
            return x, attns
        return x
