"""Train VIRAD with the four-stage curriculum (Framework §XII)."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import torch
from torch.utils.data import DataLoader

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from data.dataset import ViralRNADataset, collate_windows, load_processed_windows
from virad.config import VIRADConfig
from virad.data.splits import stratified_split
from virad.models.virad import VIRAD
from virad.seed import set_seed_from_config
from virad.training.artifacts import RunDirectory
from virad.training.budget import write_json
from virad.training.checkpoint import resolve_resume_path
from virad.training.loop import run_curriculum
from virad.training.speed import (
    apply_fast_preset,
    apply_speed_settings,
    effective_batch_size,
    make_dataloader,
    maybe_compile,
    maybe_distributed_sampler,
)
from virad.training.trainer import Trainer


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-dir", default=None)
    parser.add_argument("--run-dir", default=None)
    parser.add_argument("--resume", default="none", help="auto | PATH | none")
    parser.add_argument("--checkpoint-every", type=int, default=None)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--config", default="configs/paper.yaml")
    parser.add_argument("--seed", type=int, default=None, help="Override config seed for one training.")
    parser.add_argument("--data", default="data/processed")
    parser.add_argument("--fast", action="store_true", help="TF32/AMP/fused Adam/workers; keeps B_eff=32 and paper epochs.")
    parser.add_argument("--compile", action="store_true", help="Optional torch.compile (eager fallback on failure).")
    parser.add_argument("--amp", default=None, help="auto | bf16 | fp16 | off")
    parser.add_argument("--num-workers", type=int, default=None)
    parser.add_argument("--micro-batch", type=int, default=None, help="Micro-batch; accum is set so B_eff stays 32.")
    parser.add_argument("--profile-steps", type=int, default=0)
    parser.add_argument("--profile-only", action="store_true")
    parser.add_argument("--no-tf32", action="store_true")
    parser.add_argument("--no-sdpa", action="store_true")
    return parser.parse_args(argv)


def _apply_cli_speed(cfg: VIRADConfig, args) -> None:
    if args.fast:
        apply_fast_preset(cfg)
    if args.compile:
        cfg.torch_compile = True
    if args.amp:
        cfg.amp = args.amp
    if args.num_workers is not None:
        cfg.num_workers = args.num_workers
    if args.no_tf32:
        cfg.allow_tf32 = False
    if args.no_sdpa:
        cfg.use_sdpa = False
    if args.micro_batch:
        cfg.micro_batch_size = int(args.micro_batch)
        target = cfg.batch_size if cfg.scale == "paper" else effective_batch_size(cfg)
        if cfg.scale == "paper":
            target = 32
        if target % cfg.micro_batch_size != 0:
            raise SystemExit(f"--micro-batch {cfg.micro_batch_size} must divide effective batch {target}")
        cfg.grad_accum_steps = target // cfg.micro_batch_size
    if args.checkpoint_every is not None:
        cfg.checkpoint_every = args.checkpoint_every


def _profile(trainer: Trainer, loader: DataLoader, n: int) -> dict:
    opt, _ = trainer.make_optimizer(1)
    times = []
    it = iter(loader)
    for _ in range(max(n, 1)):
        try:
            batch = next(it)
        except StopIteration:
            it = iter(loader)
            batch = next(it)
        t0 = time.perf_counter()
        trainer.train_batch(batch, opt, stage=1)
        if trainer.device.type == "cuda":
            torch.cuda.synchronize()
        times.append(time.perf_counter() - t0)
    skip = 1 if len(times) > 2 else 0
    core = times[skip:]
    mean_t = sum(core) / max(len(core), 1)
    return {
        "n": n,
        "seconds": times,
        "s_per_step": mean_t,
        "steps_per_sec": (1.0 / mean_t) if mean_t else None,
        "note": "First steps include warmup; band is large. Do not treat as a 187 h replacement.",
    }


def main(argv=None) -> int:
    args = parse_args(argv)
    cfg_path = Path(args.config)
    cfg = VIRADConfig.from_yaml(cfg_path) if cfg_path.is_file() else VIRADConfig()
    if args.seed is not None:
        cfg.seed = int(args.seed)
    _apply_cli_speed(cfg, args)
    set_seed_from_config(cfg)
    records = load_processed_windows(Path(args.data), cfg)
    if not records:
        print("No training windows found. Populate data/processed and retry.")
        return 1
    if cfg.cache_windows:
        cache = Path(args.data) / "windows.pt"
        if not cache.is_file():
            try:
                torch.save(records, cache)
            except OSError as exc:
                print("windows.pt cache not written:", exc)
    for rec in records:
        if "has_motif" not in rec:
            mask = rec.get("motif_mask") or []
            rec["has_motif"] = any(float(v) > 0 for v in mask)
    declared = {str(r.get("partition") or "") for r in records}
    frozen = declared and declared <= {"train", "val", "test"} and "" not in declared
    strict_paper = cfg.scale == "paper"
    if strict_paper:
        if not frozen:
            raise SystemExit(
                "paper windows are missing a frozen cluster partition. "
                "Training will not invent a second split or treat held-out rows as train."
            )
        for rec in records:
            host = rec.get("host") or []
            if not any(abs(float(value)) > 0.0 for value in host):
                name = rec.get("seq_id") or rec.get("name") or "window"
                raise SystemExit(
                    f"ESM-2 host vector for {name} is missing or all zeros. "
                    "Paper training does not substitute a zero host vector."
                )
            if cfg.require_complete_structure:
                matrix = rec.get("structure_matrix")
                if matrix is None:
                    name = rec.get("seq_id") or rec.get("name") or "window"
                    raise SystemExit(f"RNA structure matrix for {name} is missing.")
                tensor = torch.as_tensor(matrix)
                if tensor.ndim != 2 or tensor.shape[0] != tensor.shape[1] or not tensor.isfinite().all():
                    name = rec.get("seq_id") or rec.get("name") or "window"
                    raise SystemExit(f"RNA structure matrix for {name} has invalid shape/content.")
                # UFold binary maps may legitimately contain no predicted pairs.
    if frozen:
        split = {name: [r for r in records if r.get("partition") == name] for name in ("train", "val", "test")}
    else:
        split = stratified_split(
            records,
            train_frac=cfg.train_fraction,
            val_frac=cfg.val_fraction,
            seed=cfg.split_seed,
        )
    cluster_sets = {
        name: {str(r.get("cluster_id") or r.get("seq_id") or r.get("name")) for r in rows}
        for name, rows in split.items()
    }
    if any(cluster_sets[a] & cluster_sets[b] for a, b in (("train", "val"), ("train", "test"), ("val", "test"))):
        raise RuntimeError("cluster leakage detected across train/val/test")
    if strict_paper and not split["train"]:
        raise SystemExit("paper training split is empty")
    train_recs = split["train"] if strict_paper else (split["train"] or records)
    val_recs = split.get("val") or []
    test_recs = split.get("test") or []
    device = torch.device(args.device)
    speed = apply_speed_settings(cfg, device)
    train_ds = ViralRNADataset(train_recs, cfg=cfg, length=cfg.n_res_train)
    sampler = maybe_distributed_sampler(train_ds, shuffle=True)
    loader = make_dataloader(train_ds, cfg, shuffle=True, device=device, collate_fn=collate_windows, sampler=sampler)
    val_loader = None
    if val_recs:
        val_ds = ViralRNADataset(val_recs, cfg=cfg, length=cfg.n_res_train)
        val_loader = make_dataloader(val_ds, cfg, shuffle=False, device=device, collate_fn=collate_windows)

    run = RunDirectory(args.run_dir) if args.run_dir else None
    ckpt_dir = Path(args.checkpoint_dir) if args.checkpoint_dir else (
        run.checkpoints if run else Path("checkpoints/paper")
    )
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    if run:
        run.snapshot_config(cfg, cfg_path if cfg_path.is_file() else None)
        run.write_meta(
            ROOT,
            extra={
                "split": {"train": len(train_recs), "val": len(val_recs), "test": len(test_recs)},
                "speed": speed,
            },
        )
        run.write_refs({"processed": str(Path(args.data)), "checkpoint_dir": str(ckpt_dir)})

    model = VIRAD(cfg)
    trainer = Trainer(model, cfg, device=device)
    trainer.student = maybe_compile(trainer.student, cfg)
    resume_path = resolve_resume_path(args.resume, ckpt_dir)
    every = args.checkpoint_every if args.checkpoint_every is not None else cfg.checkpoint_every

    if args.profile_steps or args.profile_only:
        n = args.profile_steps or 5
        prof = _profile(trainer, loader, n)
        print(json.dumps({"profile": prof, "speed": speed}))
        if run:
            write_json(run.root / "profile.json", prof)
        if args.profile_only:
            return 0

    def log_row(row: dict) -> None:
        print(json.dumps(row))
        if run:
            run.append_metric(row)

    metrics = run_curriculum(
        trainer,
        loader,
        val_loader,
        ckpt_dir,
        cfg,
        resume_path=resume_path,
        checkpoint_every=every,
        metric_logger=log_row,
    )
    if run:
        run.write_summary(metrics)
        run.write_artifacts_index()
        run.mark_phase("train", "done", checkpoint_dir=str(ckpt_dir))
    print(json.dumps({"ok": True, "checkpoint_dir": str(ckpt_dir), "speed": speed, **{k: metrics[k] for k in metrics if k in ("stage", "global_epoch", "best_val_loss")}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
