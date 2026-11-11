#!/usr/bin/env python3
"""Validation-only interval reconstruction grid for VIRAD.

This script does not train the model and must not use the held-out test split.
It runs the frozen checkpoint over the validation partition, caches per-base
presence probabilities, and evaluates one-to-one interval F1 after thresholding,
fragment merging, and minimum-length filtering.

The inference cache is resumable. If execution is interrupted, rerunning this
script continues from the number of validation windows already stored in the
partial cache.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset


DEFAULT_REPO = Path("/workspace/runpod_local/VIRAD-2026")
DEFAULT_CONFIG = Path(
    "/workspace/VIRAD-2026/final_freeze/significant-130/significant-130.yaml"
)
DEFAULT_CHECKPOINT = Path(
    "/workspace/VIRAD-2026/final_freeze/significant-130/best.pt"
)
DEFAULT_PROCESSED = Path("/workspace/runpod_local/virad_processed_fast")
DEFAULT_OUT = Path("/workspace/VIRAD-2026/validation_interval_grid")
EXPECTED_VALIDATION_WINDOWS = 12852

THRESHOLDS = [
    0.50,
    0.60,
    0.65,
    0.70,
    0.75,
    0.80,
    0.82,
    0.84,
    0.86,
    0.88,
    0.90,
    0.92,
    0.94,
    0.96,
    0.98,
]
MERGE_GAPS = [0, 1, 2, 4, 8, 12, 16]
MIN_LENGTHS = [1, 2, 4, 6, 8, 12, 16, 20]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validation-only VIRAD interval postprocessing grid."
    )
    parser.add_argument("--repo", type=Path, default=DEFAULT_REPO)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--processed", type=Path, default=DEFAULT_PROCESSED)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--batch-size",
        type=int,
        default=16,
        help="Inference batch size. Drop to 8 if CUDA runs out of memory.",
    )
    parser.add_argument(
        "--cache-every",
        type=int,
        default=100,
        help="Save resumable partial probability cache every N inference batches.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help="Validation DataLoader workers.",
    )
    return parser.parse_args()


def spans_from_mask(mask: np.ndarray) -> list[tuple[int, int]]:
    x = np.asarray(mask, dtype=bool).reshape(-1)
    if x.size == 0:
        return []

    padded = np.concatenate(
        [np.array([False]), x, np.array([False])]
    )
    diff = np.diff(padded.astype(np.int8))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    return list(zip(starts.tolist(), ends.tolist()))


def merge_spans(
    spans: list[tuple[int, int]], gap: int
) -> list[tuple[int, int]]:
    if not spans:
        return []

    merged = [list(spans[0])]
    for start, end in spans[1:]:
        previous = merged[-1]
        separation = start - previous[1]
        if separation <= gap:
            previous[1] = max(previous[1], end)
        else:
            merged.append([start, end])

    return [(int(start), int(end)) for start, end in merged]


def min_length_filter(
    spans: list[tuple[int, int]], minimum: int
) -> list[tuple[int, int]]:
    return [(a, b) for a, b in spans if (b - a) >= minimum]


def safe_torch_load(path: Path):
    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    repo = args.repo.resolve()
    sys.path.insert(0, str(repo))
    sys.path.insert(0, str(repo / "src"))

    from data.dataset import (  # noqa: E402
        ViralRNADataset,
        collate_windows,
        load_processed_windows,
    )
    from virad.config import VIRADConfig  # noqa: E402
    from virad.evaluation.metrics import f1_score, presence_confusion  # noqa: E402
    from virad.models.virad import VIRAD  # noqa: E402
    from virad.training.trainer import Trainer  # noqa: E402

    cache_path = args.out / "validation_probabilities.pt"
    partial_cache_path = args.out / "validation_probabilities_partial.pt"

    cfg = VIRADConfig.from_yaml(args.config)
    cfg.micro_batch_size = args.batch_size
    cfg.num_workers = args.num_workers
    cfg.amp = "bf16"

    records = load_processed_windows(args.processed, cfg)
    validation_records = [
        record
        for record in records
        if str(record.get("partition")) == "val"
    ]

    print("validation windows:", len(validation_records), flush=True)
    if len(validation_records) != EXPECTED_VALIDATION_WINDOWS:
        raise SystemExit(
            "STOP: expected "
            f"{EXPECTED_VALIDATION_WINDOWS} validation windows, "
            f"got {len(validation_records)}"
        )

    if cache_path.exists():
        print("Loading complete probability cache:", cache_path, flush=True)
        cache = safe_torch_load(cache_path)
        probs = cache["probs"]
        truths = cache["truths"]
    else:
        dataset = ViralRNADataset(
            validation_records,
            cfg=cfg,
            length=cfg.n_res_train,
        )

        existing_probs: list[torch.Tensor] = []
        existing_truths: list[torch.Tensor] = []
        completed_windows = 0

        if partial_cache_path.exists():
            print(
                "Loading resumable partial cache:", partial_cache_path, flush=True
            )
            partial = safe_torch_load(partial_cache_path)
            partial_probs = partial["probs"]
            partial_truths = partial["truths"]
            completed_windows = int(partial_probs.shape[0])

            if partial_truths.shape[0] != completed_windows:
                raise SystemExit("STOP: partial cache tensor lengths disagree")
            if completed_windows > len(dataset):
                raise SystemExit("STOP: partial cache exceeds validation dataset size")

            existing_probs.append(partial_probs)
            existing_truths.append(partial_truths)
            print(
                f"Resuming after {completed_windows}/{len(dataset)} windows",
                flush=True,
            )

        remaining_indices = range(completed_windows, len(dataset))
        remaining_dataset = Subset(dataset, remaining_indices)
        loader = DataLoader(
            remaining_dataset,
            batch_size=args.batch_size,
            shuffle=False,
            num_workers=args.num_workers,
            collate_fn=collate_windows,
        )

        print(
            f"Running frozen model on {len(remaining_dataset)} remaining windows "
            f"with batch_size={args.batch_size}...",
            flush=True,
        )

        checkpoint_blob = safe_torch_load(args.checkpoint)
        model = VIRAD(cfg)
        model.load_state_dict(checkpoint_blob["student"], strict=True)
        trainer = Trainer(model, cfg, device=torch.device("cuda"))
        trainer.student.load_state_dict(checkpoint_blob["student"], strict=True)
        trainer.student.eval()

        new_prob_chunks: list[torch.Tensor] = []
        new_truth_chunks: list[torch.Tensor] = []

        with torch.no_grad():
            for batch_index, batch in enumerate(loader, start=1):
                moved = trainer._move_batch(batch)
                with trainer._amp_ctx():
                    output = trainer.student(
                        moved["sequence"],
                        moved["conservation"],
                        moved["structure"],
                        moved["host"],
                        attn_mask=moved.get("attention_mask"),
                        n_refine=cfg.n_refine,
                    )

                probability = (
                    output["p_pres"].squeeze(-1).detach().float().cpu()
                )
                truth = batch["motif_mask"].detach().bool().cpu()
                new_prob_chunks.append(probability)
                new_truth_chunks.append(truth)

                if batch_index % 50 == 0 or batch_index == len(loader):
                    processed_new = sum(x.shape[0] for x in new_prob_chunks)
                    processed_total = completed_windows + processed_new
                    print(
                        f"inference {batch_index}/{len(loader)} batches | "
                        f"windows {processed_total}/{len(dataset)}",
                        flush=True,
                    )

                if (
                    batch_index % args.cache_every == 0
                    or batch_index == len(loader)
                ):
                    all_prob_parts = existing_probs + new_prob_chunks
                    all_truth_parts = existing_truths + new_truth_chunks
                    partial_probs = torch.cat(all_prob_parts, dim=0)
                    partial_truths = torch.cat(all_truth_parts, dim=0)
                    torch.save(
                        {
                            "probs": partial_probs,
                            "truths": partial_truths,
                            "completed_windows": int(partial_probs.shape[0]),
                            "batch_size": args.batch_size,
                            "checkpoint": str(args.checkpoint),
                            "config": str(args.config),
                            "partition": "val",
                        },
                        partial_cache_path,
                    )
                    print(
                        "partial cache saved: "
                        f"{partial_probs.shape[0]}/{len(dataset)} windows",
                        flush=True,
                    )

        all_prob_parts = existing_probs + new_prob_chunks
        all_truth_parts = existing_truths + new_truth_chunks
        probs = torch.cat(all_prob_parts, dim=0)
        truths = torch.cat(all_truth_parts, dim=0)

        if probs.shape[0] != len(dataset):
            raise SystemExit(
                f"STOP: expected {len(dataset)} cached probability rows, "
                f"got {probs.shape[0]}"
            )

        torch.save(
            {
                "probs": probs,
                "truths": truths,
                "checkpoint": str(args.checkpoint),
                "config": str(args.config),
                "partition": "val",
                "n_windows": len(validation_records),
                "batch_size": args.batch_size,
            },
            cache_path,
        )
        print("Complete probability cache saved:", cache_path, flush=True)

        if partial_cache_path.exists():
            partial_cache_path.unlink()

    if probs.shape[0] != EXPECTED_VALIDATION_WINDOWS:
        raise SystemExit(
            "STOP: probability cache has wrong number of validation windows: "
            f"{probs.shape[0]}"
        )
    if truths.shape[0] != probs.shape[0]:
        raise SystemExit("STOP: probability and truth cache lengths disagree")

    print("probability tensor:", tuple(probs.shape), flush=True)
    print("truth tensor:", tuple(truths.shape), flush=True)

    probability_array = probs.numpy()
    truth_array = truths.numpy()

    truth_intervals = []
    positive_windows = 0
    true_interval_count = 0
    for truth_mask in truth_array:
        intervals = spans_from_mask(truth_mask)
        truth_intervals.append(intervals)
        true_interval_count += len(intervals)
        if intervals:
            positive_windows += 1

    print("positive validation windows:", positive_windows, flush=True)
    print("true intervals:", true_interval_count, flush=True)

    def evaluate(
        threshold: float, merge_gap: int, min_length: int
    ) -> dict[str, float | int]:
        tp = fp = fn = pred_count = 0

        for index, probability_row in enumerate(probability_array):
            raw = spans_from_mask(probability_row > threshold)
            reconstructed = merge_spans(raw, merge_gap)
            reconstructed = min_length_filter(reconstructed, min_length)

            tpi, fpi, fni = presence_confusion(
                reconstructed, truth_intervals[index]
            )
            tp += tpi
            fp += fpi
            fn += fni
            pred_count += len(reconstructed)

        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = f1_score(tp, fp, fn)

        return {
            "threshold": threshold,
            "merge_gap": merge_gap,
            "min_length": min_length,
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "pred_intervals": pred_count,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }

    print("\n=== BASELINE SANITY CHECK ===", flush=True)
    baseline = evaluate(threshold=0.50, merge_gap=0, min_length=1)
    for key, value in baseline.items():
        print(f"{key}: {value}", flush=True)

    print("Expected baseline validation F1 is approximately 0.1398.", flush=True)
    if abs(float(baseline["f1"]) - 0.1398) > 0.01:
        print(
            "WARNING: baseline differs materially from the previous "
            "validation result.",
            flush=True,
        )

    total = len(THRESHOLDS) * len(MERGE_GAPS) * len(MIN_LENGTHS)
    print("\n=== STARTING VALIDATION GRID ===", flush=True)
    print("combinations:", total, flush=True)

    results: list[dict[str, float | int]] = []
    counter = 0
    for threshold in THRESHOLDS:
        for gap in MERGE_GAPS:
            for min_len in MIN_LENGTHS:
                counter += 1
                result = evaluate(threshold, gap, min_len)
                results.append(result)

                if counter % 25 == 0 or counter == total:
                    best_so_far = max(
                        results,
                        key=lambda row: (
                            float(row["f1"]),
                            float(row["precision"]),
                        ),
                    )
                    print(
                        f"{counter}/{total} | current t={threshold:.2f} "
                        f"gap={gap} min={min_len} | "
                        f"best F1={float(best_so_far['f1']):.4f} "
                        f"(t={float(best_so_far['threshold']):.2f}, "
                        f"gap={best_so_far['merge_gap']}, "
                        f"min={best_so_far['min_length']})",
                        flush=True,
                    )

    results.sort(
        key=lambda row: (
            float(row["f1"]),
            float(row["precision"]),
            float(row["recall"]),
            -int(row["pred_intervals"]),
        ),
        reverse=True,
    )

    csv_path = args.out / "interval_grid_results.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "threshold",
                "merge_gap",
                "min_length",
                "tp",
                "fp",
                "fn",
                "pred_intervals",
                "precision",
                "recall",
                "f1",
            ],
        )
        writer.writeheader()
        writer.writerows(results)

    best = results[0]
    summary = {
        "status": "validation_only_interval_postprocessing",
        "checkpoint": str(args.checkpoint),
        "partition": "val",
        "n_validation_windows": len(validation_records),
        "true_intervals": true_interval_count,
        "baseline": baseline,
        "best": best,
        "top_25": results[:25],
    }
    summary_path = args.out / "grid_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print("\n" + "=" * 70, flush=True)
    print("BEST VALIDATION POSTPROCESSING", flush=True)
    print("=" * 70, flush=True)
    for key, value in best.items():
        print(f"{key}: {value}", flush=True)

    print("\n" + "=" * 70, flush=True)
    print("TOP 20 COMBINATIONS", flush=True)
    print("=" * 70, flush=True)
    print(
        "rank  thr   gap  min     TP     FP     FN   PredInt   Prec     Rec      F1",
        flush=True,
    )
    for rank, row in enumerate(results[:20], start=1):
        print(
            f"{rank:>4}  "
            f"{float(row['threshold']):.2f}  "
            f"{int(row['merge_gap']):>3}  "
            f"{int(row['min_length']):>3}  "
            f"{int(row['tp']):>5}  "
            f"{int(row['fp']):>5}  "
            f"{int(row['fn']):>5}  "
            f"{int(row['pred_intervals']):>7}  "
            f"{float(row['precision']):.4f}  "
            f"{float(row['recall']):.4f}  "
            f"{float(row['f1']):.4f}",
            flush=True,
        )

    absolute = float(best["f1"]) - float(baseline["f1"])
    relative = (
        absolute / float(baseline["f1"]) * 100.0
        if float(baseline["f1"])
        else 0.0
    )
    print("\n" + "=" * 70, flush=True)
    print("IMPROVEMENT OVER ORIGINAL", flush=True)
    print("=" * 70, flush=True)
    print("original F1:", baseline["f1"], flush=True)
    print("best F1:", best["f1"], flush=True)
    print("absolute improvement:", absolute, flush=True)
    print("relative improvement (%):", relative, flush=True)

    print("\nOutputs:", flush=True)
    print(cache_path, flush=True)
    print(csv_path, flush=True)
    print(summary_path, flush=True)


if __name__ == "__main__":
    main()
