"""Sliding-window generation and ensemble averaging (Framework §III.a, §IV.a, §XIII.a)."""

from __future__ import annotations

import math
from typing import List, Optional, Sequence

import torch


def n_windows(genome_length: int, window: int, stride: int) -> int:
    if genome_length <= 0 or window <= 0 or stride <= 0:
        return 0
    if genome_length <= window:
        return 1
    return int(math.ceil((genome_length - window) / stride)) + 1


def window_starts(genome_length: int, window: int, stride: int) -> List[int]:
    count = n_windows(genome_length, window, stride)
    starts = []
    for k in range(count):
        start = k * stride
        if start + window > genome_length:
            start = max(0, genome_length - window)
        starts.append(start)
    return starts


def fragment_at(sequence: str, start: int, window: int) -> tuple[int, str]:
    n = len(sequence)
    start = max(0, min(start, max(0, n - window))) if n > window else 0
    frag = sequence[start : start + window]
    if len(frag) < window:
        frag = frag + "N" * (window - len(frag))
    return start, frag


def extract_windows(sequence: str, window: int, stride: int) -> List[tuple[int, str]]:
    n = len(sequence)
    out = []
    for start in window_starts(n, window, stride):
        out.append(fragment_at(sequence, start, window))
    return out


def motif_centered_windows(
    sequence: str,
    window: int,
    intervals: Sequence[tuple[int, int]],
) -> List[tuple[int, str]]:
    """One L-nt window centered on each motif midpoint (Framework $80 crop)."""
    n = len(sequence)
    seen = set()
    out: List[tuple[int, str]] = []
    for a, b in intervals:
        mid = (int(a) + int(b)) // 2
        start = max(0, mid - window // 2)
        start, frag = fragment_at(sequence, start, window)
        if start in seen:
            continue
        seen.add(start)
        out.append((start, frag))
        if n <= window:
            break
    if not out:
        return extract_windows(sequence, window, max(window // 2, 1))
    return out


def cap_windows(windows: Sequence[tuple[int, str]], max_windows: Optional[int]) -> List[tuple[int, str]]:
    items = list(windows)
    if max_windows is None or max_windows <= 0 or len(items) <= max_windows:
        return items
    if max_windows == 1:
        return [items[0]]
    span = len(items) - 1
    idxs = []
    for i in range(max_windows):
        idxs.append(int(round(i * span / (max_windows - 1))))
    unique = []
    seen = set()
    for i in idxs:
        if i not in seen:
            seen.add(i)
            unique.append(items[i])
    return unique


def overlaps_interval(start: int, end: int, intervals: Sequence[tuple[int, int]]) -> bool:
    for a, b in intervals:
        if start < b and end > a:
            return True
    return False


def min_distance_to_motifs(start: int, end: int, intervals: Sequence[tuple[int, int]]) -> int:
    if not intervals:
        return 10**9
    mid = (start + end) // 2
    dist = 10**9
    for a, b in intervals:
        if start < b and end > a:
            return 0
        dist = min(dist, abs(start - b), abs(end - a), abs(mid - a), abs(mid - b))
    return dist


def ensemble_average(window_probs: Sequence[tuple[int, torch.Tensor]], genome_length: int) -> torch.Tensor:
    """Eq. 1: mean of window probabilities covering each genome position."""
    acc = torch.zeros(genome_length, dtype=torch.float32)
    counts = torch.zeros(genome_length, dtype=torch.float32)
    for start, probs in window_probs:
        length = int(probs.numel())
        end = min(genome_length, start + length)
        sl = end - start
        acc[start:end] += probs[:sl]
        counts[start:end] += 1.0
    return acc / counts.clamp_min(1.0)
