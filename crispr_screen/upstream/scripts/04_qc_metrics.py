#!/usr/bin/env python3
"""
04_qc_metrics.py — Dynamic Quality Control Metrics
====================================================

Computes ALL QC metrics dynamically from count data. No hardcoded
values — a critical fix from previous pipeline versions (v0.x).

Metrics per sample:
  • Total mapped reads (M)
  • Gini inequality index (0 = uniform, 1 = max inequality)
  • Zero-count sgRNA percentage
  • Number of detected sgRNAs (non-zero counts)
  • Detection rate (% of library detected)

QC grades assigned per sample:
  A — All metrics within ideal range
  B — 1–2 metrics in warning range
  C — >2 metrics in warning range
  D — >2 metrics in critical range

Output:
    results/reports/qc_metrics.csv
    results/reports/qc_report.txt
    results/reports/metadata_04_qc_metrics.json

Project: PRJ25-131 | Cingöz Lab | 2026

Usage:
    python 04_qc_metrics.py
    python 04_qc_metrics.py --count-file count/count_table_fixed.txt
    python 04_qc_metrics.py --min-reads 5000000
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _thesis_utils import (
    ThesisConfig,
    classify_condition,
    compute_gini,
    export_metadata,
    load_count_table,
    setup_logging,
)

logger = setup_logging("04_qc_metrics")

# QC grade colour codes for terminal output
_GRADE_COLORS = {"A": "✓", "B": "~", "C": "⚠", "D": "✗"}


def compute_qc_metrics(
    count_matrix: pd.DataFrame,
    *,
    min_reads: int = 9_000_000,
    ideal_gini: float = 0.2,
    ideal_zero_pct: float = 20.0,
    ideal_detected_pct: float = 90.0,
) -> pd.DataFrame:
    """
    Compute per-sample QC metrics from a count matrix.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Count matrix (rows = sgRNAs, columns = samples).
    min_reads : int
        Minimum mapped reads per sample (adequate depth).
    ideal_gini : float
        Maximum acceptable Gini index.
    ideal_zero_pct : float
        Maximum acceptable zero-count percentage.
    ideal_detected_pct : float
        Minimum acceptable detected sgRNA percentage.

    Returns
    -------
    pd.DataFrame
        Per-sample QC metrics with flags and grades.
    """
    total_sgrnas = len(count_matrix)
    rows: list[dict] = []

    for col in count_matrix.columns:
        vals = count_matrix[col].values.astype(np.float64)
        total_reads = vals.sum()
        nonzero = int((vals > 0).sum())
        zero_pct = float((vals == 0).mean() * 100)
        gini = compute_gini(vals)
        condition = classify_condition(col)

        flag_low_depth = total_reads < min_reads
        flag_high_gini = gini > ideal_gini
        flag_high_zero = zero_pct > ideal_zero_pct
        flag_low_detection = (nonzero / total_sgrnas * 100) < ideal_detected_pct

        n_critical = sum([flag_low_depth, flag_high_gini, flag_high_zero, flag_low_detection])
        if n_critical == 0:
            grade = "A"
        elif n_critical <= 2:
            grade = "B"
        elif n_critical <= 3:
            grade = "C"
        else:
            grade = "D"

        rows.append({
            "sample": col,
            "condition": condition,
            "total_reads": int(total_reads),
            "total_reads_M": round(total_reads / 1e6, 2),
            "nonzero_sgrnas": nonzero,
            "nonzero_sgrna_pct": round(nonzero / total_sgrnas * 100, 2),
            "zero_pct": round(zero_pct, 2),
            "gini": round(gini, 4),
            "total_sgrnas": total_sgrnas,
            "qc_grade": grade,
            "flag_low_depth": flag_low_depth,
            "flag_high_gini": flag_high_gini,
            "flag_high_zero": flag_high_zero,
            "flag_low_detection": flag_low_detection,
        })

    return pd.DataFrame(rows)


def write_qc_report(qc_df: pd.DataFrame, report_path: Path, config: ThesisConfig) -> None:
    """Write a human-readable QC report for the thesis supplement."""
    lines: list[str] = []
    w = lines.append

    w("CRISPR Screen Quality Control Report")
    w("=" * 60)
    w("Project:    PRJ25-131 — GBM TRAIL Metabolic CRISPR/Cas9 Screen")
    w("Library:    Sabatini Human Metabolic KO Library (Addgene #110066)")
    w(f"Date:       {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}")
    w(f"Samples:    {len(qc_df)}  (4 conditions: R0, Rpost, S0, Spost)")
    w("")
    w("QC Thresholds (reference):")
    w(f"  Reads per sample:      >{config.min_reads_per_sample / 1e6:.0f}M (adequate)")
    w(f"  Gini index:            <{config.ideal_gini} (uniform distribution)")
    w(f"  Zero-count sgRNAs:     <{config.ideal_zero_pct}%")
    w(f"  Detected sgRNAs:       >{config.ideal_detected_pct}% of library")
    w("")

    # Per-sample table
    w("Per-Sample Metrics:")
    w("-" * 60)
    hdr = f"{'Sample':<10} {'Cond':<7} {'Reads(M)':>8} {'Nonzero':>8} {'Zero%':>7} {'Gini':>7} {'Det%':>6} {'Grade':>6}"
    w(hdr)
    w("-" * len(hdr))
    for _, row in qc_df.iterrows():
        w(f"{row['sample']:<10} {row['condition']:<7} {row['total_reads_M']:>8.2f} "
          f"{row['nonzero_sgrnas']:>8} {row['zero_pct']:>7.2f} {row['gini']:>7.4f} "
          f"{row['nonzero_sgrna_pct']:>6.1f} {row['qc_grade']:>6}")
    w("")

    # Condition summaries
    w("Condition Summaries:")
    w("-" * 40)
    for cond in ["R0", "Rpost", "S0", "Spost"]:
        sub = qc_df[qc_df["condition"] == cond]
        if sub.empty:
            continue
        w(f"\n{cond} (n={len(sub)}):")
        w(f"  Reads:  {sub['total_reads_M'].mean():.2f}M "
          f"(range {sub['total_reads_M'].min():.2f}–{sub['total_reads_M'].max():.2f})")
        w(f"  Gini:   {sub['gini'].mean():.4f} "
          f"(range {sub['gini'].min():.4f}–{sub['gini'].max():.4f})")
        w(f"  Zero%:  {sub['zero_pct'].mean():.1f}%")
        w(f"  Grade:  {sub['qc_grade'].value_counts().to_dict()}")

    # Overall assessment
    w("")
    w("=" * 60)
    grade_counts = qc_df["qc_grade"].value_counts().to_dict()
    w(f"Overall QC Grade Distribution: {grade_counts}")

    if qc_df["qc_grade"].eq("D").any():
        w("")
        w("⚠  CRITICAL QC FAILURES — The following samples are unsuitable for")
        w("    definitive biological conclusions:")
        for _, row in qc_df[qc_df["qc_grade"] == "D"].iterrows():
            reasons = []
            if row["flag_low_depth"]:
                reasons.append(f"low depth ({row['total_reads_M']:.1f}M)")
            if row["flag_high_gini"]:
                reasons.append(f"high Gini ({row['gini']:.3f})")
            if row["flag_high_zero"]:
                reasons.append(f"high zero% ({row['zero_pct']:.1f}%)")
            if row["flag_low_detection"]:
                reasons.append(f"low detection ({row['nonzero_sgrna_pct']:.1f}%)")
            w(f"    {row['sample']} ({row['condition']}): {', '.join(reasons)}")
        w("")
        w("    This screen suffers from a severe wet-lab library bottleneck.")
        w("    Results should be treated as EXPLORATORY and require orthogonal")
        w("    validation (qPCR, Western blot, or independent CRISPR).")

    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    logger.info("QC report written to %s", report_path)


def main() -> None:
    cfg = ThesisConfig()

    parser = argparse.ArgumentParser(
        description="Compute CRISPR screen QC metrics dynamically from count data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--count-file", type=Path,
        default=Path(__file__).resolve().parent / "count" / "count_table_fixed.txt",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(__file__).resolve().parent / "reports",
    )
    parser.add_argument("--min-reads", type=int, default=cfg.min_reads_per_sample)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("STEP 4 — QC Metrics (Dynamic Computation)")
    logger.info("=" * 60)
    logger.info("Count table: %s", args.count_file)
    logger.info("Output:      %s", args.outdir)
    logger.info("")

    # ── Load data ───────────────────────────────────────────────────────────
    try:
        counts_raw, count_matrix, sample_groups = load_count_table(args.count_file)
    except FileNotFoundError as e:
        logger.error(str(e))
        return

    logger.info("Loaded: %d sgRNAs × %d samples", len(count_matrix), len(count_matrix.columns))

    # ── Compute QC ──────────────────────────────────────────────────────────
    qc_df = compute_qc_metrics(
        count_matrix,
        min_reads=args.min_reads,
        ideal_gini=cfg.ideal_gini,
        ideal_zero_pct=cfg.ideal_zero_pct,
        ideal_detected_pct=cfg.ideal_detected_pct,
    )

    # ── Print summary ───────────────────────────────────────────────────────
    logger.info("")
    logger.info("Per-Sample QC:")
    for _, row in qc_df.iterrows():
        grade_icon = _GRADE_COLORS.get(row["qc_grade"], "?")
        logger.info(
            "  %s %-6s  %-6s  %6.1fM reads  %5d sgRNAs  %5.1f%% zero  "
            "Gini=%.3f  Det=%.1f%%",
            grade_icon, row["sample"], row["condition"],
            row["total_reads_M"], row["nonzero_sgrnas"],
            row["zero_pct"], row["gini"], row["nonzero_sgrna_pct"],
        )

    # ── Grade distribution ──────────────────────────────────────────────────
    logger.info("")
    for grade in ["A", "B", "C", "D"]:
        n = int((qc_df["qc_grade"] == grade).sum())
        if n > 0:
            logger.info("  Grade %s: %d sample(s)", grade, n)

    critical = qc_df[qc_df["qc_grade"] == "D"]
    if not critical.empty:
        logger.warning("")
        logger.warning("CRITICAL QC FAILURES (%d samples):", len(critical))
        for _, row in critical.iterrows():
            logger.warning("  %s (%s): Gini=%.3f, zero=%.1f%%, detected=%d",
                          row["sample"], row["condition"],
                          row["gini"], row["zero_pct"], row["nonzero_sgrnas"])

    # ── Save ────────────────────────────────────────────────────────────────
    csv_path = args.outdir / "qc_metrics.csv"
    qc_df.to_csv(csv_path, index=False)
    logger.info("")
    logger.info("QC CSV: %s", csv_path)

    report_path = args.outdir / "qc_report.txt"
    write_qc_report(qc_df, report_path, cfg)

    meta = {
        "n_samples": len(qc_df),
        "n_total_sgrnas": int(qc_df["total_sgrnas"].iloc[0]),
        "grade_distribution": qc_df["qc_grade"].value_counts().to_dict(),
        "mean_gini": round(float(qc_df["gini"].mean()), 4),
        "mean_zero_pct": round(float(qc_df["zero_pct"].mean()), 2),
        "mean_detected_pct": round(float(qc_df["nonzero_sgrna_pct"].mean()), 2),
        "qc_thresholds": {
            "min_reads": args.min_reads,
            "ideal_gini": cfg.ideal_gini,
            "ideal_zero_pct": cfg.ideal_zero_pct,
            "ideal_detected_pct": cfg.ideal_detected_pct,
        },
    }
    export_metadata(args.outdir, step="04_qc_metrics", extra=meta, config=cfg)
    logger.info("Done.")


if __name__ == "__main__":
    main()
