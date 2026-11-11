"""Run inference windows of length 1024 / stride 256 (Framework §XIII)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from virad.config import VIRADConfig
from virad.data.annotations import parse_fasta
from virad.data.tracks import load_conservation
from virad.inference.pipeline import (
    InferenceEngine,
    attach_casper_jasper,
    load_structure_npz,
    load_structure_square,
    stitch_window_structures,
)
from virad.models.virad import VIRAD
from virad.seed import set_seed_from_config


def load_query_sequences(args) -> list[tuple[str, str, str]]:
    if args.fasta:
        recs = parse_fasta(Path(args.fasta).read_text())
        if not recs:
            raise SystemExit("empty FASTA")
        return recs
    return [(args.seq_id or "query", args.sequence, f">{args.seq_id or 'query'}")]


def _record_inputs(args, cfg, seq_id: str, n: int):
    cons_src = args.phylop or args.conservation
    cons = (
        torch.tensor(load_conservation(cons_src, seq_id, n), dtype=torch.float32)
        if cons_src
        else torch.zeros(n)
    )
    structures = {}
    host = torch.zeros(cfg.host_esm_dim)
    if args.processed:
        root = Path(args.processed)
        host_path = root / seq_id / "host.json"
        if host_path.is_file():
            host = torch.as_tensor(json.loads(host_path.read_text()), dtype=torch.float32)
        engine = "ufold" if str(cfg.structure_predictor).lower() == "ufold" else "spotrna"
        for path in sorted((root / seq_id / engine).glob("*/results.npz")):
            try:
                structures[int(path.parent.name)] = load_structure_square(path)
            except (ValueError, OSError, KeyError):
                continue
    if args.structure:
        structure = load_structure_npz(args.structure, n)
    elif structures:
        structure = structures
    else:
        structure = torch.zeros(n, n)
    return cons, structure, host


def run_infer(args) -> dict:
    cfg_path = Path(args.config)
    cfg = VIRADConfig.from_yaml(cfg_path) if cfg_path.is_file() else VIRADConfig()
    set_seed_from_config(cfg)
    model = VIRAD(cfg)
    if args.checkpoint:
        blob = torch.load(args.checkpoint, map_location="cpu")
        model.load_state_dict(blob["student"])
    engine = InferenceEngine(model, cfg, device=torch.device(args.device))
    out_dir = Path(args.out_dir) if args.out_dir else None
    jasper = Path(args.jasper) if args.jasper else (out_dir / "jasper33.json" if out_dir else Path("jasper33.json"))
    cluster_rows = {}
    if args.processed:
        manifest = Path(args.processed) / "window_manifest.json"
        if manifest.is_file():
            for row in json.loads(manifest.read_text()).get("records", []):
                cluster_rows.setdefault(row["seq_id"], row)
    all_motifs, all_guides, record_payloads = [], [], []
    records = load_query_sequences(args)
    for seq_id, sequence, header in records:
        n = len(sequence)
        cons, struct, host = _record_inputs(args, cfg, seq_id, n)
        if cfg.scale == "paper" and cfg.require_complete_structure and seq_id in cluster_rows:
            if isinstance(struct, dict):
                if not struct:
                    raise SystemExit(f"RNA structure matrices missing for {seq_id}.")
                for matrix in struct.values():
                    tensor = torch.as_tensor(matrix)
                    if tensor.ndim != 2 or tensor.shape[0] != tensor.shape[1] or not tensor.isfinite().all():
                        raise SystemExit(f"RNA structure matrix for {seq_id} has invalid shape/content.")
            else:
                tensor = torch.as_tensor(struct)
                if tensor.ndim != 2 or tensor.shape != (n, n) or not tensor.isfinite().all():
                    raise SystemExit(f"RNA structure matrix for {seq_id} has invalid shape/content.")
            # UFold binary contact maps may legitimately be all zero.
            if not torch.as_tensor(host).count_nonzero().item():
                raise SystemExit(
                    f"ESM-2 host vector for {seq_id} is missing or all zeros. "
                    "Paper inference does not substitute a zero host vector."
                )
        result = engine.predict_genome(sequence, cons, struct, host)
        attach_structure = (
            struct if torch.is_tensor(struct) else stitch_window_structures(struct, n)
        )
        family = args.family or next(
            (part.split("=", 1)[1] for part in header.split() if part.startswith("family=")),
            "unknown",
        )
        rec_payload = attach_casper_jasper(
            result,
            sequence,
            attach_structure,
            cfg,
            virus=args.virus or seq_id,
            family=family,
            jasper_path=None,
        )
        cluster = cluster_rows.get(seq_id, {})
        for motif in rec_payload["motifs"]:
            motif["seq_id"] = seq_id
            motif["cluster_id"] = cluster.get("cluster_id", seq_id)
            motif["family"] = family
        for guide in rec_payload["guides"]:
            guide["seq_id"] = seq_id
        record_payloads.append(
            {
                "seq_id": seq_id,
                "cluster_id": cluster.get("cluster_id", seq_id),
                "family": family,
                "n_motifs": len(rec_payload["motifs"]),
            }
        )
        all_motifs.extend(rec_payload["motifs"])
        all_guides.extend(rec_payload["guides"])
    payload = {
        "n_records": len(records),
        "records": record_payloads,
        "n_motifs": len(all_motifs),
        "motifs": all_motifs,
        "guides": all_guides,
        "jasper": str(jasper),
        "n_refine": cfg.n_refine,
    }
    jasper.parent.mkdir(parents=True, exist_ok=True)
    jasper.write_text(json.dumps(payload, indent=2))
    if out_dir is not None:
        from virad.training.artifacts import write_prediction_exports

        payload["exports"] = write_prediction_exports(out_dir, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", default=None)
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--seq-id", default=None)
    parser.add_argument("--checkpoint", default=None)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--config", default="configs/paper.yaml")
    parser.add_argument("--conservation", default=None, help="Alias for --phylop.")
    parser.add_argument("--phylop", default=None, help="phyloP wig/bedGraph/npy/npz/dir.")
    parser.add_argument("--structure", default=None, help="RNA structure results.npz with ss_prob_map.")
    parser.add_argument("--processed", default=None, help="Preprocess root with per-record structure/host inputs.")
    parser.add_argument("--jasper", default=None)
    parser.add_argument("--out-dir", default=None, help="Write predictions/JASPER/CASPER/CSV/FASTA here.")
    parser.add_argument("--virus", default=None)
    parser.add_argument("--family", default=None)
    args = parser.parse_args()
    if not args.sequence and not args.fasta:
        raise SystemExit("provide --sequence or --fasta")
    payload = run_infer(args)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
