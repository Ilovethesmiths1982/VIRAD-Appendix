"""Build model-ready windows, labels, conservation, structure, and host features."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from virad.config import VIRADConfig
from virad.data.prepare import preprocess


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", default="data/processed")
    parser.add_argument("--config", default="configs/paper.yaml")
    parser.add_argument(
        "--spot-rna",
        default=None,
        dest="spot_rna_script",
        help="Compatibility path to SPOT-RNA.py when structure_predictor=spotrna.",
    )
    parser.add_argument(
        "--spot-python",
        default=None,
        help="Compatibility Python interpreter for SPOT-RNA when structure_predictor=spotrna.",
    )
    parser.add_argument("--gff", default=None, help="Alias for --annotations (GFF/GFF3).")
    parser.add_argument("--annotations", default=None, help="Motif GFF, GFF3, BED, CSV, or TSV.")
    parser.add_argument(
        "--phylop",
        default=None,
        help="phyloP file or directory (wig/wigFix/bedGraph/npy/npz/bigWig).",
    )
    parser.add_argument(
        "--skip-spot-rna",
        action="store_true",
        help="Compatibility alias for --skip-structure.",
    )
    parser.add_argument(
        "--skip-structure",
        action="store_true",
        help="Skip the configured RNA-structure predictor.",
    )
    parser.add_argument("--motif-centered-only", action="store_true")
    parser.add_argument("--max-windows-per-record", type=int, default=None)
    parser.add_argument("--max-total-windows", type=int, default=None)
    parser.add_argument("--ann-coords", default=None, choices=("pdf", "gff", "bed"))
    parser.add_argument(
        "--esm",
        action="store_true",
        help="Embed host proteins with ESM-2 650M. Paper runs require non-zero host vectors.",
    )
    parser.add_argument("--host-fasta", default=None)
    args = parser.parse_args()

    cfg_path = Path(args.config)
    cfg = VIRADConfig.from_yaml(cfg_path) if cfg_path.is_file() else VIRADConfig()
    if args.motif_centered_only:
        cfg.motif_centered_only = True
    if args.max_windows_per_record is not None:
        cfg.max_windows_per_record = args.max_windows_per_record
    if args.max_total_windows is not None:
        cfg.max_total_windows = args.max_total_windows
        cfg.max_train_windows = args.max_total_windows
    if args.ann_coords:
        cfg.ann_coords = args.ann_coords
    if args.spot_python:
        cfg.spot_rna_python = args.spot_python

    ann = args.annotations or args.gff
    preprocess(
        fasta=Path(args.fasta),
        out=Path(args.out),
        cfg=cfg,
        annotations=Path(ann) if ann else None,
        phylop=Path(args.phylop) if args.phylop else None,
        skip_structure=args.skip_spot_rna or args.skip_structure,
        spot_rna_script=args.spot_rna_script,
        esm=args.esm,
        host_fasta=Path(args.host_fasta) if args.host_fasta else None,
    )


if __name__ == "__main__":
    main()
