"""Post-train validation: GFF recovery, hide-one-family prep, Rfam split, scaling manifests."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from virad.config import VIRADConfig
from virad.evaluation.metrics import stratified_paired_bootstrap_delta
from virad.evaluation.validate import (
    average_seed_cluster_f1,
    default_report_path,
    load_gff_spans,
    load_pred_spans,
    load_processed_label_spans,
    grouped_recovery_report,
    prepare_holdout,
    recovery_report,
    rfam_novel_split,
    sequences_for_partition,
    write_scaling_subsets,
    _norm_id,
)


def _maybe_infer(args) -> Path | None:
    if args.jasper and Path(args.jasper).is_file() and not args.run_infer:
        return Path(args.jasper)
    if not (args.checkpoint and args.fasta):
        return Path(args.jasper) if args.jasper else None
    import importlib.util
    from argparse import Namespace

    spec = importlib.util.spec_from_file_location("virad_infer_cli", ROOT / "scripts" / "infer.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    dest = Path(args.jasper or "jasper33.json")
    ns = Namespace(
        sequence=None,
        fasta=args.fasta,
        seq_id=args.seq_id,
        checkpoint=args.checkpoint,
        device=args.device,
        config=args.config,
        conservation=None,
        phylop=args.phylop,
        structure=args.structure,
        processed=args.processed,
        jasper=str(dest),
        out_dir=None,
        virus=args.seq_id,
        family=args.holdout_family or args.family,
    )
    payload = mod.run_infer(ns)
    dest.write_text(json.dumps(payload, indent=2))
    return dest


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Simple VIRAD validations (does not launch GPU trains; not a LOFO trainer)."
    )
    parser.add_argument("--jasper", default=None, help="Infer/JASPER JSON to score.")
    parser.add_argument("--report", default=None, help="Output JSON (default: next to --jasper).")
    parser.add_argument("--gff", default=None, help="Training/label GFF3.")
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--processed", default=None, help="preprocess output dir (window labels).")
    parser.add_argument("--phylop", default=None)
    parser.add_argument("--seq-id", default=None, help="Default seq id when JASPER omits seq_id.")
    parser.add_argument("--holdout-family", default=None)
    parser.add_argument("--prepare-holdout", action="store_true")
    parser.add_argument("--holdout-dir", default="data/raw/holdout")
    parser.add_argument("--rfam-cm", default=None)
    parser.add_argument(
        "--rfam-tblout",
        default=None,
        help="Fetch-time cmscan tblout. When this file exists, validation parses it and does not scan again.",
    )
    parser.add_argument("--fractions", default=None, help="Comma-separated e.g. 0.1,0.25,0.5,1.0")
    parser.add_argument("--scaling-dir", default="data/raw/scaling")
    parser.add_argument("--config", default="configs/significant-130.yaml")
    parser.add_argument("--checkpoint", default=None)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--structure", default=None)
    parser.add_argument("--run-infer", action="store_true")
    parser.add_argument("--family", default=None)
    parser.add_argument("--baseline-jasper", default=None, help="Locked strongest-baseline predictions for paired cluster bootstrap.")
    parser.add_argument("--aggregate-dir", default=None, help="Directory of seed-*/validate.json reports to average.")
    args = parser.parse_args()

    if args.aggregate_dir:
        tables = []
        for path in sorted(Path(args.aggregate_dir).glob("seed-*/validate.json")):
            payload = json.loads(path.read_text())
            per_cluster = (payload.get("cluster_validation") or {}).get("per_cluster") or {}
            if per_cluster:
                tables.append(per_cluster)
        report = average_seed_cluster_f1(tables)
        dest = Path(args.report) if args.report else Path(args.aggregate_dir) / "seed-average.json"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"ok": True, "n_seeds": report["n_seeds"], "macro_f1": report["macro_f1"]}))
        return

    cfg_path = Path(args.config)
    cfg = VIRADConfig.from_yaml(cfg_path) if cfg_path.is_file() else VIRADConfig()
    report: dict = {
        "notes": [
            "novel_to_labels and rfam_unmatched are not biologically novel.",
            "Trainer split stays 70/15/15; hide-one-family is an extra retrain you run yourself.",
            "P is SPOT-RNA v1.0. No RNAfold on train. No LOFO training split.",
        ],
        "max_total_windows": cfg.max_total_windows,
    }

    jasper_path = _maybe_infer(args)
    pred = []
    if jasper_path and jasper_path.is_file():
        pred = load_pred_spans(jasper_path, default_seq=args.seq_id or "")
        report["jasper"] = str(jasper_path)

    truth = []
    if args.gff:
        truth = load_gff_spans(args.gff)
        report["labels_from"] = args.gff
    elif args.processed:
        truth = load_processed_label_spans(args.processed)
        report["labels_from"] = args.processed

    if pred or truth:
        metadata = {}
        manifest_rows = []
        if args.processed:
            manifest = Path(args.processed) / "window_manifest.json"
            if manifest.is_file():
                manifest_rows = json.loads(manifest.read_text()).get("records", [])
                seen_part: dict[str, str] = {}
                for row in manifest_rows:
                    seq_id = _norm_id(row.get("seq_id") or "")
                    part = str(row.get("partition") or "")
                    if seq_id in seen_part and seen_part[seq_id] != part:
                        raise SystemExit(
                            f"{seq_id} is labeled both {seen_part[seq_id]} and {part}"
                        )
                    if seq_id:
                        seen_part[seq_id] = part
                    metadata.setdefault(
                        seq_id,
                        {
                            "cluster_id": row.get("cluster_id"),
                            "family": row.get("family"),
                            "partition": row.get("partition"),
                        },
                    )
        if str(cfg.methodology).startswith("significant-") and manifest_rows:
            test_ids = sequences_for_partition(manifest_rows, "test")
            truth = [span for span in truth if span.seq_id in test_ids]
            pred = [span for span in pred if span.seq_id in test_ids]
            metadata = {seq_id: row for seq_id, row in metadata.items() if seq_id in test_ids}
            report["scored_partition"] = "test"
        report["recovery"] = recovery_report(
            truth, pred, min_len=cfg.motif_min_len, max_len=cfg.motif_max_len
        )
        report["cluster_validation"] = grouped_recovery_report(truth, pred, metadata)
        if args.baseline_jasper:
            baseline_pred = load_pred_spans(args.baseline_jasper, default_seq=args.seq_id or "")
            baseline_grouped = grouped_recovery_report(truth, baseline_pred, metadata)
            virad_clusters = report["cluster_validation"]["per_cluster"]
            baseline_clusters = baseline_grouped["per_cluster"]
            common = sorted(set(virad_clusters) & set(baseline_clusters))
            report["baseline_cluster_validation"] = baseline_grouped
            report["paired_cluster_bootstrap"] = stratified_paired_bootstrap_delta(
                [virad_clusters[c]["f1"] for c in common],
                [baseline_clusters[c]["f1"] for c in common],
                [virad_clusters[c]["family"] for c in common],
                seed=cfg.split_seed,
            )

    if args.holdout_family and (args.prepare_holdout or (args.fasta and args.gff)):
        if not args.fasta or not args.gff:
            raise SystemExit("--holdout-family needs --fasta and --gff")
        hold = prepare_holdout(
            args.fasta,
            args.gff,
            args.holdout_family,
            out_dir=args.holdout_dir,
            phylop=args.phylop,
            max_total_windows=cfg.max_total_windows,
        )
        hold_ids = set(hold["holdout_ids"])
        hold_truth = [s for s in load_gff_spans(hold["holdout_gff"])]
        hold_pred = [
            s
            for s in pred
            if s.seq_id in hold_ids or str(s.meta.get("family") or "").lower() == args.holdout_family.lower()
        ]
        hold["recovery"] = recovery_report(
            hold_truth, hold_pred, min_len=cfg.motif_min_len, max_len=cfg.motif_max_len
        )
        hold["n_pred_holdout"] = len(hold_pred)
        if not jasper_path or not Path(jasper_path).is_file():
            hold["scoring"] = "no jasper passed; CPU path is prepare-holdout + score later"
        report["holdout"] = hold

    if pred and (args.gff or args.processed):
        report["rfam_novel"] = rfam_novel_split(
            pred,
            truth,
            fasta=args.fasta,
            rfam_cm=args.rfam_cm,
            tblout=args.rfam_tblout,
            min_len=cfg.motif_min_len,
            max_len=cfg.motif_max_len,
        )

    if args.fractions:
        if not args.fasta or not args.gff:
            raise SystemExit("--fractions needs --fasta and --gff")
        fracs = [float(x.strip()) for x in args.fractions.split(",") if x.strip()]
        report["scaling"] = write_scaling_subsets(
            args.fasta,
            args.gff,
            fracs,
            out_dir=args.scaling_dir,
            seed=cfg.seed,
            phylop=args.phylop,
            max_total_windows=cfg.max_total_windows,
        )

    dest = Path(args.report) if args.report else default_report_path(jasper_path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    from virad.training.artifacts import atomic_write_text

    atomic_write_text(dest, json.dumps(report, indent=2) + "\n")
    print(json.dumps({"report": str(dest), **{k: report[k] for k in report if k != "notes"}}, indent=2)[:4000])
    print("wrote", dest)


if __name__ == "__main__":
    main()
