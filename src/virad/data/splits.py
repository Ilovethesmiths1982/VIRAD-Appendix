"""Family- and motif-stratified 70/15/15 split with CD-HIT leakage guard (§IV.c)."""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Sequence


def stratified_split(
    records: Sequence[dict],
    train_frac: float = 0.70,
    val_frac: float = 0.15,
    seed: int = 42,
) -> Dict[str, List[dict]]:
    # Split indivisible genome/CD-HIT groups, never overlapping windows.
    clusters: Dict[str, List[dict]] = defaultdict(list)
    for index, rec in enumerate(records):
        cluster_id = str(rec.get("cluster_id") or rec.get("seq_id") or rec.get("name") or f"record-{index}")
        clusters[cluster_id].append(rec)
    grouped: Dict[tuple, List[tuple[str, List[dict]]]] = defaultdict(list)
    for cluster_id, members in clusters.items():
        families = {m.get("family") for m in members if m.get("family")}
        family = sorted(families)[0] if families else None
        grouped[(family, any(bool(m.get("has_motif")) for m in members))].append(
            (cluster_id, members)
        )

    rng_order = sorted(grouped.keys())
    train, val, test = [], [], []
    # Deterministic shuffle per stratum using the fixed seed without numpy.
    import random

    rng = random.Random(seed)
    for key in rng_order:
        items = list(grouped[key])
        rng.shuffle(items)
        n = len(items)
        n_train = int(round(train_frac * n))
        n_val = int(round(val_frac * n))
        train.extend(rec for _, members in items[:n_train] for rec in members)
        val.extend(rec for _, members in items[n_train : n_train + n_val] for rec in members)
        test.extend(rec for _, members in items[n_train + n_val :] for rec in members)
    return {"train": train, "val": val, "test": test}


def _cluster_ids(rows: Sequence[dict]) -> set[str]:
    return {str(row.get("cluster_id") or row.get("seq_id") or row.get("name")) for row in rows}


def family_holdout_split(
    records: Sequence[dict],
    holdout_family: str,
    train_frac: float = 0.85,
    val_frac: float = 0.15,
    seed: int = 42,
) -> Dict[str, List[dict]]:
    """Remove one family from train and validation.

    The other families are split into an inner train/validation pair for early
    stopping. Any rounding remainder from that inner split stays in train.
    Evaluation rows are only the held-out family's clusters.
    """
    target = holdout_family.strip().lower()
    held, rest = [], []
    for record in records:
        family = str(record.get("family") or "")
        (held if family.lower() == target else rest).append(record)
    inner = stratified_split(rest, train_frac=train_frac, val_frac=val_frac, seed=seed)
    train = list(inner["train"]) + list(inner["test"])
    val = list(inner["val"])
    test = list(held)
    for bucket in (train, val):
        leaked = [row for row in bucket if str(row.get("family") or "").lower() == target]
        if leaked:
            raise RuntimeError(f"held-out family {holdout_family!r} leaked into train/val")
    train_ids, val_ids, test_ids = _cluster_ids(train), _cluster_ids(val), _cluster_ids(test)
    if train_ids & val_ids or train_ids & test_ids or val_ids & test_ids:
        raise RuntimeError("family holdout reused a cluster across partitions")
    return {"train": train, "val": val, "test": test}


def leave_one_family_out(
    records: Sequence[dict],
    families: Sequence[str],
    seed: int = 42,
) -> Dict[str, Dict[str, List[dict]]]:
    """One split per family. Structure caches are not part of the split."""
    return {family: family_holdout_split(records, family, seed=seed) for family in families}


def split_assignments(records: Sequence[dict], **kwargs) -> Dict[str, str]:
    """Return a frozen cluster -> partition map and reject duplicate assignment."""
    split = stratified_split(records, **kwargs)
    assignments: Dict[str, str] = {}
    for partition, rows in split.items():
        for row in rows:
            cluster_id = str(row.get("cluster_id") or row.get("seq_id") or row.get("name"))
            prior = assignments.setdefault(cluster_id, partition)
            if prior != partition:
                raise ValueError(f"cluster {cluster_id!r} assigned to both {prior} and {partition}")
    return assignments
