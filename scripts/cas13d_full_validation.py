#!/usr/bin/env python3
"""Leakage-controlled computational validation for VIRAD Cas13d guide panels.

Workflow
--------
1. develop: TRAIN+VALIDATION ONLY, source-excluded, family-specific evidence.
2. freeze: freeze guide sequences, panel composition, selection rule, and hashes.
3. test: evaluate the frozen panel once on the untouched TEST partition.

The development stage rejects test-derived guide candidates and partition-filters both
whole-genome and homologous-locus analyses. Test sequences therefore cannot influence
ranking, panel selection, or threshold choice through either occurrence-based coverage
or family MSAs.

Optional modules add human-transcript similarity screening, RNA accessibility and
energetics, complete-crRNA folding, external activity scores, matched controls, and
ranking-sensitivity summaries. All outputs are computational evidence; no efficacy
probability is fabricated.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

DNA_COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
DEV_PARTITIONS = {"train", "val", "validation"}
TEST_PARTITIONS = {"test"}


def norm(seq: str) -> str:
    return re.sub(r"[^ACGT]", "N", (seq or "").upper().replace("U", "T"))


def revcomp(seq: str) -> str:
    return norm(seq).translate(DNA_COMP)[::-1]


def fnum(x, default=float("nan")):
    try:
        return float(x)
    except Exception:
        return default


def inum(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fasta(path: Path):
    name = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, norm("".join(seq))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, norm("".join(seq))


def aligned_fasta(path: Path):
    name = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper().replace("U", "T")
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper().replace("U", "T")


def read_csv(path: Path):
    return list(csv.DictReader(path.open()))


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                fields.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def load_manifest(path: Path):
    rows = read_csv(path)
    out = {}
    for r in rows:
        sid = r.get("seq_id") or r.get("accession") or r.get("id")
        if sid:
            out[sid] = {
                "family": r.get("family", ""),
                "partition": (r.get("partition") or "").lower(),
                "cluster_id": r.get("cluster_id", ""),
            }
    return out


def pct_rank(vals):
    clean = [(i, v) for i, v in enumerate(vals) if not math.isnan(v)]
    out = [float("nan")] * len(vals)
    if not clean:
        return out
    order = sorted(clean, key=lambda x: x[1])
    n, i = len(order), 0
    while i < n:
        j = i
        while j + 1 < n and order[j + 1][1] == order[i][1]:
            j += 1
        r = ((i + j) / 2) / max(n - 1, 1)
        for k in range(i, j + 1):
            out[order[k][0]] = r
        i = j + 1
    return out


def best_window(query: str, genome: str):
    q = norm(query)
    L = len(q)
    if L == 0 or len(genome) < L:
        return None
    best = (L + 1, -1, "+", "")
    for strand, qq in (("+", q), ("-", revcomp(q))):
        for i in range(len(genome) - L + 1):
            s = genome[i:i + L]
            mm = sum(a != b for a, b in zip(qq, s))
            if mm < best[0]:
                best = (mm, i, strand, s)
                if mm == 0:
                    return best
    return best


def coverage_for_guide(guide, genomes, source_id=None):
    records = []
    for sid, seq in genomes:
        if source_id and sid == source_id:
            continue
        hit = best_window(guide, seq)
        if hit is None:
            records.append((sid, None, None, None, ""))
        else:
            mm, pos, strand, target = hit
            records.append((sid, mm, pos, strand, target))
    n = len(records)
    out = {"viral_genomes_n": n}
    for k, label in [(0, "exact"), (1, "one_mm"), (2, "two_mm")]:
        c = sum(mm is not None and mm <= k for _, mm, *_ in records)
        out[f"{label}_n"] = c
        out[f"{label}_pct"] = c / n if n else float("nan")
    return out, records


def shannon_entropy(chars):
    vals = [x for x in chars if x in "ACGT"]
    if not vals:
        return float("nan")
    n = len(vals)
    c = Counter(vals)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def family_msa_path(msa_dir: Path, family: str):
    candidates = [
        msa_dir / f"{family}.fa",
        msa_dir / f"{family}.fasta",
        msa_dir / f"{family}.aln.fa",
        msa_dir / f"{family}.aln.fasta",
    ]
    return next((p for p in candidates if p.exists()), None)


def ungapped_to_alignment_cols(aligned_source: str, start0: int, length: int):
    cols = []
    ungapped = 0
    wanted = set(range(start0, start0 + length))
    for col, ch in enumerate(aligned_source):
        if ch != "-":
            if ungapped in wanted:
                cols.append(col)
            ungapped += 1
    return cols if len(cols) == length else None


def locus_profile(guide_row, msa_dir: Path | None, manifest, allowed_partitions, exclude_source=False):
    """Evaluate only homologous loci from explicitly allowed partitions."""
    if not msa_dir:
        return {}, []
    fam = guide_row.get("family", "")
    src = guide_row.get("seq_id", "")
    start0 = inum(guide_row.get("guide_start"), -1)
    guide = norm(guide_row.get("sequence", ""))
    p = family_msa_path(msa_dir, fam)
    if not p or start0 < 0 or not guide:
        return {}, []
    aln = dict(aligned_fasta(p))
    if src not in aln:
        return {}, []
    cols = ungapped_to_alignment_cols(aln[src], start0, len(guide))
    if not cols:
        return {}, []

    variants = []
    for sid, aseq in aln.items():
        meta = manifest.get(sid)
        if not meta:
            continue
        if meta.get("family") != fam:
            continue
        if meta.get("partition") not in allowed_partitions:
            continue
        if exclude_source and sid == src:
            continue
        frag = "".join(aseq[c] if c < len(aseq) else "-" for c in cols)
        if "-" in frag:
            variants.append((sid, None, frag, []))
            continue
        frag = norm(frag)
        direct_m = [i for i, (a, b) in enumerate(zip(guide, frag)) if a != b]
        rc = revcomp(guide)
        rc_m = [i for i, (a, b) in enumerate(zip(rc, frag)) if a != b]
        mismatch_positions = direct_m if len(direct_m) <= len(rc_m) else rc_m
        variants.append((sid, len(mismatch_positions), frag, mismatch_positions))

    usable = [(sid, mm, frag, pos) for sid, mm, frag, pos in variants if mm is not None]
    n = len(usable)
    out = {"locus_genomes_n": n}
    for k, label in [(0, "locus_exact"), (1, "locus_one_mm"), (2, "locus_two_mm")]:
        c = sum(mm <= k for _, mm, _, _ in usable)
        out[f"{label}_n"] = c
        out[f"{label}_pct"] = c / n if n else float("nan")
    if usable:
        ent = [shannon_entropy([frag[i] for _, _, frag, _ in usable]) for i in range(len(guide))]
        finite = [x for x in ent if not math.isnan(x)]
        out["locus_mean_entropy"] = sum(finite) / len(finite) if finite else float("nan")
        out["locus_max_entropy"] = max(finite) if finite else float("nan")
        haps = Counter(frag for _, _, frag, _ in usable)
        out["locus_haplotype_count"] = len(haps)
        out["locus_major_haplotype_pct"] = max(haps.values()) / n
        pos_counts = Counter(i for _, _, _, positions in usable for i in positions)
        out["mismatch_position_counts_json"] = json.dumps(dict(sorted(pos_counts.items())))
        out["mismatch_count_distribution_json"] = json.dumps(dict(sorted(Counter(mm for _, mm, _, _ in usable).items())))
    return out, variants


def rnaplfold_accessibility(context: str, offset: int, L: int):
    exe = shutil.which("RNAplfold")
    if not exe or not context:
        return float("nan")
    seq = norm(context).replace("T", "U")
    with tempfile.TemporaryDirectory(prefix="virad_plfold_") as td:
        td = Path(td)
        cp = subprocess.run(
            [exe, "-W", str(min(150, len(seq))), "-L", str(min(100, max(1, len(seq)-1))), "-u", "28"],
            cwd=td, input=">x\n" + seq + "\n", text=True, capture_output=True,
        )
        if cp.returncode:
            return float("nan")
        fs = list(td.glob("*_lunp"))
        end1 = offset + L
        if not fs:
            return float("nan")
        for line in fs[0].read_text().splitlines():
            p = line.split()
            if p and p[0].isdigit() and int(p[0]) == end1 and len(p) > L:
                return fnum(p[L])
    return float("nan")


def rnafold_mfe(seq: str):
    exe = shutil.which("RNAfold")
    if not exe or not seq:
        return float("nan")
    cp = subprocess.run([exe, "--noPS"], input=seq.replace("T", "U") + "\n", text=True, capture_output=True)
    m = re.search(r"\(([-+]?\d+(?:\.\d+)?)\)\s*$", cp.stdout.strip())
    return float(m.group(1)) if m else float("nan")


def duplex_energy(target: str):
    exe = shutil.which("RNAduplex")
    if not exe or not target:
        return float("nan")
    antisense = revcomp(target)
    cp = subprocess.run([exe], input=antisense.replace("T", "U") + "\n" + target.replace("T", "U") + "\n", text=True, capture_output=True)
    m = re.search(r"\(([-+]?\d+(?:\.\d+)?)\)", cp.stdout)
    return float(m.group(1)) if m else float("nan")


def full_crrna_mfe(target: str, direct_repeat: str | None):
    if not direct_repeat:
        return float("nan")
    return rnafold_mfe(norm(direct_repeat) + revcomp(target))


def bowtie_screen(rows, transcriptome: Path | None, outdir: Path):
    bowtie = shutil.which("bowtie2")
    build = shutil.which("bowtie2-build")
    if not bowtie or not build or not transcriptome:
        return {}
    idx = outdir / "human_tx"
    if not Path(str(idx) + ".1.bt2").exists() and not Path(str(idx) + ".1.bt2l").exists():
        subprocess.run([build, str(transcriptome), str(idx)], check=True)
    fa = outdir / "guides.fa"
    sam = outdir / "human_offtargets.sam"
    with fa.open("w") as f:
        for i, r in enumerate(rows):
            f.write(f">g{i}\n{norm(r.get('sequence',''))}\n")
    subprocess.run([bowtie, "-x", str(idx), "-U", str(fa), "-f", "--end-to-end", "-N", "1", "-L", "15", "-a", "-k", "100", "--quiet", "-S", str(sam)], check=True)
    hits = defaultdict(list)
    for line in sam.open():
        if line.startswith("@"):
            continue
        p = line.rstrip().split("\t")
        if len(p) < 12 or int(p[1]) & 4:
            continue
        nm = next((int(x[5:]) for x in p[11:] if x.startswith("NM:i:")), None)
        hits[p[0]].append((p[2], nm))
    return hits


def matched_control(row, genome_seq: str, radius=500, gc_tol=0.05, seed=42):
    random.seed(seed + inum(row.get("guide_start"), 0))
    g = norm(row.get("sequence", ""))
    L = len(g)
    if L == 0 or len(genome_seq) < L:
        return ""
    target_gc = (g.count("G") + g.count("C")) / L
    center = inum(row.get("guide_start"), 0)
    lo, hi = max(0, center - radius), min(len(genome_seq) - L, center + radius)
    cands = []
    for s in range(lo, hi + 1):
        if abs(s - center) < L:
            continue
        q = genome_seq[s:s+L]
        qgc = (q.count("G") + q.count("C")) / L
        if abs(qgc - target_gc) <= gc_tol and "N" not in q:
            cands.append(q)
    return random.choice(cands) if cands else ""


def build_population(fasta_rows, manifest, partitions, family):
    allowed = set(partitions)
    return [
        (sid, seq) for sid, seq in fasta_rows
        if manifest.get(sid, {}).get("partition") in allowed
        and (not family or manifest.get(sid, {}).get("family") == family)
    ]


def activity_map(path: Path | None):
    out = {}
    if path and path.exists():
        for r in read_csv(path):
            out[norm(r.get("sequence", ""))] = fnum(r.get("cas13d_activity_score"))
    return out


def consensus_rank(rows):
    metrics = [
        ("exact_pct", 1), ("one_mm_pct", 1), ("locus_exact_pct", 1),
        ("target_accessibility", 1), ("motif_conservation", 1),
        ("motif_confidence", 1), ("cas13d_activity_score", 1),
        ("human_offtarget_alignments", -1), ("locus_mean_entropy", -1),
    ]
    ranks = []
    for key, direction in metrics:
        vals = [fnum(r.get(key)) for r in rows]
        ranks.append(pct_rank([v * direction if not math.isnan(v) else v for v in vals]))
    for i, r in enumerate(rows):
        av = [rr[i] for rr in ranks if not math.isnan(rr[i])]
        r["consensus_score"] = sum(av) / len(av) if av else float("nan")
    rows.sort(key=lambda r: (-fnum(r.get("consensus_score"), -1), -fnum(r.get("exact_pct"), -1), -fnum(r.get("one_mm_pct"), -1), r.get("sequence", "")))
    for i, r in enumerate(rows, 1):
        r["development_rank"] = i


def sensitivity_summary(rows):
    schemes = {
        "all": ["exact_pct", "one_mm_pct", "locus_exact_pct", "target_accessibility", "motif_conservation", "motif_confidence", "cas13d_activity_score"],
        "sequence_only": ["exact_pct", "one_mm_pct", "locus_exact_pct"],
        "without_virad_inputs": ["exact_pct", "one_mm_pct", "locus_exact_pct", "cas13d_activity_score"],
        "accessibility_plus_sequence": ["exact_pct", "one_mm_pct", "locus_exact_pct", "target_accessibility"],
    }
    out = []
    for name, keys in schemes.items():
        usable = []
        for r in rows:
            vals = [fnum(r.get(k)) for k in keys]
            vals = [v for v in vals if not math.isnan(v)]
            score = sum(vals) / len(vals) if vals else float("nan")
            usable.append((score, r.get("sequence", ""), r.get("family", ""), r.get("seq_id", "")))
        usable.sort(key=lambda x: (-fnum(x[0], -1), x[1]))
        for rank, (score, seq, fam, sid) in enumerate(usable, 1):
            out.append({"scheme": name, "rank": rank, "score": score, "sequence": seq, "family": fam, "seq_id": sid})
    return out


def develop(args):
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_guides = read_csv(args.guides)
    manifest = load_manifest(args.manifest)
    viral = list(fasta(args.viral_fasta))
    viral_map = dict(viral)
    motifs = {r.get("motif_id", ""): r for r in read_csv(args.motif_catalog)} if args.motif_catalog and args.motif_catalog.exists() else {}
    act = activity_map(args.activity_csv)

    # Hard firewall: only train/validation-derived guide candidates may enter development.
    guides, rejected_test = [], []
    for r in all_guides:
        src = r.get("seq_id", "")
        part = manifest.get(src, {}).get("partition") or (r.get("partition") or "").lower()
        if part in DEV_PARTITIONS:
            guides.append(r)
        else:
            rejected_test.append({"seq_id": src, "partition": part, "sequence": r.get("sequence", ""), "reason": "not_train_or_validation"})
    write_csv(args.outdir / "CAS13D_REJECTED_NONDEVELOPMENT_GUIDES.csv", rejected_test)

    hits = bowtie_screen(guides, args.human_transcriptome, args.outdir) if args.human_transcriptome else {}
    out, controls = [], []
    for i, r in enumerate(guides):
        q = dict(r)
        fam = r.get("family", "") or manifest.get(r.get("seq_id", ""), {}).get("family", "")
        src = r.get("seq_id", "")
        g = norm(r.get("sequence", ""))
        pop = build_population(viral, manifest, DEV_PARTITIONS, fam)
        cov, _ = coverage_for_guide(g, pop, source_id=src)
        q.update(cov)
        q["development_population"] = "train+validation_source_excluded"
        q["source_excluded"] = True
        q["source_partition"] = manifest.get(src, {}).get("partition", r.get("partition", ""))
        lp, _ = locus_profile(r, args.family_msa_dir, manifest, DEV_PARTITIONS, exclude_source=True)
        q.update(lp)
        mr = motifs.get(r.get("motif_id", ""), {})
        context = norm(mr.get("sequence", ""))
        off = inum(r.get("offset"), 0)
        q["target_accessibility"] = rnaplfold_accessibility(context, off, len(g)) if context else float("nan")
        q["guide_mfe"] = rnafold_mfe(revcomp(g))
        q["duplex_energy"] = duplex_energy(g)
        q["complete_crrna_mfe"] = full_crrna_mfe(g, args.direct_repeat)
        hh = hits.get(f"g{i}", [])
        q["human_offtarget_alignments"] = len(hh)
        q["closest_human_transcript"] = min(hh, key=lambda x: 99 if x[1] is None else x[1])[0] if hh else ""
        q["closest_human_mismatches"] = min((x[1] for x in hh if x[1] is not None), default="")
        q["cas13d_activity_score"] = act.get(g, float("nan"))
        out.append(q)

        if src in viral_map:
            cseq = matched_control(r, viral_map[src], seed=args.seed)
            if cseq:
                ccov, _ = coverage_for_guide(cseq, pop, source_id=src)
                controls.append({"source_sequence": g, "control_sequence": cseq, "seq_id": src, "family": fam, **ccov})

    consensus_rank(out)
    write_csv(args.outdir / "CAS13D_DEVELOPMENT_VALIDATION.csv", out)
    write_csv(args.outdir / "CAS13D_MATCHED_CONTROLS.csv", controls)
    write_csv(args.outdir / "CAS13D_RANKING_SENSITIVITY.csv", sensitivity_summary(out))
    summary = {
        "mode": "develop",
        "input_guides": len(all_guides),
        "development_guides": len(out),
        "rejected_non_development_guides": len(rejected_test),
        "population": "train+validation only",
        "source_exclusion": True,
        "msa_partition_filter": sorted(DEV_PARTITIONS),
        "test_used": False,
        "family_msa_locus_analysis": bool(args.family_msa_dir),
        "human_transcriptome_screened": bool(args.human_transcriptome),
        "complete_crrna_folded": bool(args.direct_repeat),
        "external_activity_scores": bool(args.activity_csv),
        "output": str(args.outdir / "CAS13D_DEVELOPMENT_VALIDATION.csv"),
    }
    (args.outdir / "CAS13D_DEVELOPMENT_SUMMARY.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def freeze(args):
    rows = read_csv(args.development_csv)
    byfam = defaultdict(list)
    for r in rows:
        if (r.get("source_partition") or "").lower() not in DEV_PARTITIONS:
            raise SystemExit("Development CSV contains a non-development guide; freeze aborted.")
        byfam[r.get("family", "")].append(r)
    selected = []
    for fam, rs in sorted(byfam.items()):
        rs.sort(key=lambda r: (inum(r.get("development_rank"), 10**9), -fnum(r.get("consensus_score"), -1), r.get("sequence", "")))
        selected.extend(rs[:args.panel_size_per_family])
    args.outdir.mkdir(parents=True, exist_ok=True)
    panel = args.outdir / "CAS13D_FROZEN_PANEL.csv"
    write_csv(panel, selected)
    meta = {
        "mode": "freeze",
        "development_csv": str(args.development_csv),
        "development_sha256": sha256_file(args.development_csv),
        "panel_csv": str(panel),
        "panel_sha256": sha256_file(panel),
        "panel_size_per_family": args.panel_size_per_family,
        "selection_rule": "lowest development_rank within family; deterministic sequence tie-break",
        "test_used_for_selection": False,
        "n_selected": len(selected),
    }
    (args.outdir / "CAS13D_PANEL_FREEZE.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


def test(args):
    args.outdir.mkdir(parents=True, exist_ok=True)
    panel = read_csv(args.frozen_panel)
    freeze_meta = json.loads(args.freeze_json.read_text())
    observed_hash = sha256_file(args.frozen_panel)
    if observed_hash != freeze_meta.get("panel_sha256"):
        raise SystemExit("Frozen-panel hash mismatch: test evaluation aborted.")
    manifest = load_manifest(args.manifest)
    viral = list(fasta(args.viral_fasta))
    out, family_union = [], defaultdict(set)

    for r in panel:
        q = dict(r)
        fam = r.get("family", "")
        g = norm(r.get("sequence", ""))
        pop = build_population(viral, manifest, TEST_PARTITIONS, fam)
        cov, recs = coverage_for_guide(g, pop, source_id=None)
        q.update({f"test_{k}": v for k, v in cov.items()})
        lp, _ = locus_profile(r, args.family_msa_dir, manifest, TEST_PARTITIONS, exclude_source=False)
        q.update({f"test_{k}": v for k, v in lp.items()})
        for sid, mm, *_ in recs:
            if mm is not None and mm <= args.panel_mismatch_threshold:
                family_union[fam].add(sid)
        out.append(q)

    write_csv(args.outdir / "CAS13D_FROZEN_TEST_EVALUATION.csv", out)
    panel_rows = []
    for fam in sorted({r.get("family", "") for r in panel}):
        denom = len(build_population(viral, manifest, TEST_PARTITIONS, fam))
        covered = len(family_union[fam])
        panel_rows.append({
            "family": fam,
            "test_genomes_n": denom,
            "panel_covered_n": covered,
            "panel_covered_pct": covered / denom if denom else float("nan"),
            "mismatch_threshold": args.panel_mismatch_threshold,
        })
    write_csv(args.outdir / "CAS13D_FROZEN_TEST_PANEL_COVERAGE.csv", panel_rows)
    summary = {
        "mode": "test",
        "frozen_panel_sha256": observed_hash,
        "n_guides": len(panel),
        "test_partition_only": True,
        "msa_partition_filter": sorted(TEST_PARTITIONS),
        "panel_mismatch_threshold": args.panel_mismatch_threshold,
        "guide_output": str(args.outdir / "CAS13D_FROZEN_TEST_EVALUATION.csv"),
        "panel_output": str(args.outdir / "CAS13D_FROZEN_TEST_PANEL_COVERAGE.csv"),
    }
    (args.outdir / "CAS13D_FROZEN_TEST_SUMMARY.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def parser():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("develop")
    d.add_argument("--guides", type=Path, required=True)
    d.add_argument("--viral-fasta", type=Path, required=True)
    d.add_argument("--manifest", type=Path, required=True, help="CSV with seq_id,family,partition[,cluster_id]")
    d.add_argument("--motif-catalog", type=Path)
    d.add_argument("--family-msa-dir", type=Path)
    d.add_argument("--human-transcriptome", type=Path)
    d.add_argument("--activity-csv", type=Path)
    d.add_argument("--direct-repeat", type=str)
    d.add_argument("--seed", type=int, default=42)
    d.add_argument("--outdir", type=Path, required=True)
    d.set_defaults(func=develop)

    f = sub.add_parser("freeze")
    f.add_argument("--development-csv", type=Path, required=True)
    f.add_argument("--panel-size-per-family", type=int, default=3)
    f.add_argument("--outdir", type=Path, required=True)
    f.set_defaults(func=freeze)

    t = sub.add_parser("test")
    t.add_argument("--frozen-panel", type=Path, required=True)
    t.add_argument("--freeze-json", type=Path, required=True)
    t.add_argument("--viral-fasta", type=Path, required=True)
    t.add_argument("--manifest", type=Path, required=True)
    t.add_argument("--family-msa-dir", type=Path)
    t.add_argument("--panel-mismatch-threshold", type=int, choices=[0, 1, 2], default=2)
    t.add_argument("--outdir", type=Path, required=True)
    t.set_defaults(func=test)
    return ap


if __name__ == "__main__":
    a = parser().parse_args()
    a.func(a)
