#!/usr/bin/env python3
"""
08_scaffold_comparison.py — Scaffold Sequence Length Comparison
================================================================

Quantifies the effect of scaffold sequence length on sgRNA detection
in lentiCRISPR v1 libraries.

Background:
  • 20 bp scaffold (our analysis):  GTTTTAGAGCTAGAAATAGC
  • 19 bp scaffold (company used):   TTTTAGAGCTAGAAATAGC

The leading G in the 20 bp scaffold (from lentiCRISPRv1 U6 promoter
processing) enriches for sgRNAs ending in G. The 19 bp variant
captures a broader set but may lose specificity. This script
quantifies detection rates for both.

Method:
  For each test sample, sample N reads, grep for scaffold sequence,
  and compute detection rate as (reads_with_scaffold / total_reads).

Output:
    results/scaffold_results/scaffold_comparison.csv
    results/figures/SCAF1_scaffold_comparison.{pdf,png}
    results/scaffold_results/metadata_08_scaffold_comparison.json

Project: PRJ25-131 | Cingöz Lab | 2026

Usage:
    python 08_scaffold_comparison.py
    python 08_scaffold_comparison.py --n-reads 1000000  # full-depth scan
"""

from __future__ import annotations

import argparse
import subprocess
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from _thesis_utils import (
    ThesisConfig, ThesisStyle,
    export_metadata, save_figure, setup_logging,
)

logger = setup_logging("08_scaffold_comparison")
warnings.filterwarnings("ignore")

# Scaffold variants
SCAFFOLD_20BP = "GTTTTAGAGCTAGAAATAGC"   # Current (20 bp, includes leading G)
SCAFFOLD_19BP = "TTTTAGAGCTAGAAATAGC"    # Company (19 bp, without leading G)

# Representative samples (one per condition)
TEST_SAMPLES: dict[str, str] = {
    "R0_1":    "PRJ25-131-R01_L01_59_1.fq.gz",   # R0 replicate 1
    "Rpost_1": "PRJ25-131-R1_L01_65_1.fq.gz",    # Rpost replicate 1
    "S0_1":    "PRJ25-131-S01_L01_56_1.fq.gz",   # S0 replicate 1
    "Spost_3": "PRJ25-131-S3_L01_64_1.fq.gz",    # Spost replicate 3 (best S)
}


def count_scaffold_occurrences(
    fastq_gz: Path,
    scaffold: str,
    n_reads: int = 200_000,
) -> dict[str, int | float]:
    """
    Count reads containing the scaffold sequence via grep.

    Uses zcat + grep for speed. Each read is 4 FASTQ lines.

    Parameters
    ----------
    fastq_gz : Path
        Path to gzipped FASTQ file.
    scaffold : str
        DNA sequence to search for (case-sensitive exact match).
    n_reads : int
        Number of reads to sample.

    Returns
    -------
    dict
        {"total": int, "with_scaffold": int, "rate_pct": float}
    """
    if not fastq_gz.exists():
        return {"total": 0, "with_scaffold": 0, "rate_pct": 0.0}

    # Count total lines (4 lines per FASTQ read)
    result_total = subprocess.run(
        f"zcat '{fastq_gz}' | head -{n_reads * 4} | wc -l",
        shell=True, text=True, capture_output=True,
    )
    total_lines = int(result_total.stdout.strip()) if result_total.returncode == 0 else 0
    total_reads = total_lines // 4

    # Count reads containing scaffold
    result_scaffold = subprocess.run(
        f"zcat '{fastq_gz}' | head -{n_reads * 4} | grep -c '{scaffold}'",
        shell=True, text=True, capture_output=True,
    )
    with_scaffold = (
        int(result_scaffold.stdout.strip())
        if result_scaffold.returncode in (0, 1)  # grep returns 1 for no matches
        else 0
    )

    rate = with_scaffold / total_reads * 100 if total_reads > 0 else 0.0
    return {"total": total_reads, "with_scaffold": with_scaffold, "rate_pct": rate}


def main() -> None:
    cfg = ThesisConfig()

    parser = argparse.ArgumentParser(
        description="Compare 20 bp vs 19 bp scaffold detection rates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--fastq-dir", type=Path, default=cfg.fastq_dir)
    parser.add_argument("--n-reads", type=int, default=200_000,
                        help="Reads to sample per file (default: 200,000 for fast mode).")
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).resolve().parent / "scaffold_results")
    parser.add_argument("--fig-dir", type=Path,
                        default=Path(__file__).resolve().parent / "figures")
    args = parser.parse_args()

    ThesisStyle.apply(font_family=cfg.font_family, dpi=cfg.figure_dpi)
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("STEP 8 — Scaffold Sequence Comparison: 20 bp vs 19 bp")
    logger.info("=" * 60)
    logger.info("FASTQ directory: %s", args.fastq_dir)
    logger.info("Reads sampled:   %,d per file", args.n_reads)
    logger.info("Test samples:    %d", len(TEST_SAMPLES))
    logger.info("")

    results: list[dict] = []
    n_missing = 0

    for sample_label, fastq_name in TEST_SAMPLES.items():
        fastq_path = args.fastq_dir / fastq_name
        if not fastq_path.exists():
            logger.warning("  [MISSING] %s", fastq_name)
            n_missing += 1
            continue

        logger.info("Sample: %s (%s)", sample_label, fastq_name)

        for scaffold, scaffold_name in [
            (SCAFFOLD_20BP, "20 bp (our analysis)"),
            (SCAFFOLD_19BP, "19 bp (company)"),
        ]:
            res = count_scaffold_occurrences(fastq_path, scaffold, n_reads=args.n_reads)
            results.append({
                "sample": sample_label,
                "scaffold": scaffold_name,
                "scaffold_seq": scaffold,
                "total_reads": res["total"],
                "with_scaffold": res["with_scaffold"],
                "rate_pct": round(float(res["rate_pct"]), 2),
            })
            logger.info("  %-20s: %,d / %,d = %.2f%%",
                       scaffold_name, res["with_scaffold"],
                       res["total"], res["rate_pct"])

    if not results:
        logger.error("No FASTQ files found. Cannot perform comparison.")
        return

    # ── Save CSV ────────────────────────────────────────────────────────────
    df_results = pd.DataFrame(results)
    df_results.to_csv(args.outdir / "scaffold_comparison.csv", index=False)
    logger.info("\n[OK] scaffold_comparison.csv saved")

    # ── Plot ────────────────────────────────────────────────────────────────
    samples = df_results["sample"].unique()
    x = np.arange(len(samples))
    w = 0.30

    fig, ax = plt.subplots(figsize=(8, 4.5))
    df_20 = df_results[df_results["scaffold"].str.startswith("20")].set_index("sample")
    df_19 = df_results[df_results["scaffold"].str.startswith("19")].set_index("sample")

    rates_20 = [df_20.loc[s, "rate_pct"] if s in df_20.index else 0 for s in samples]
    rates_19 = [df_19.loc[s, "rate_pct"] if s in df_19.index else 0 for s in samples]

    bars_20 = ax.bar(x - w/2, rates_20, width=w, color="#0072B2",
                     label="20 bp scaffold (our analysis)", edgecolor="none")
    bars_19 = ax.bar(x + w/2, rates_19, width=w, color="#D55E00",
                     label="19 bp scaffold (company)", edgecolor="none")

    for bar in bars_20:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.15, f"{h:.1f}%",
                    ha="center", va="bottom", fontsize=7.5, color="#0072B2",
                    fontweight="bold")
    for bar in bars_19:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.15, f"{h:.1f}%",
                    ha="center", va="bottom", fontsize=7.5, color="#D55E00",
                    fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(samples, fontsize=9)
    ax.set_ylabel("% Reads Containing Scaffold")
    ax.set_title(
        "Scaffold Detection Rate: 20 bp vs 19 bp\n"
        f"(sampled {args.n_reads:,} reads per file; lentiCRISPR v1 backbone)",
        fontsize=10, fontweight="bold",
    )
    ax.legend(fontsize=8)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, "SCAF1_scaffold_comparison", args.fig_dir)

    # ── Summary statistics ──────────────────────────────────────────────────
    valid_20 = [r for r in rates_20 if r > 0]
    valid_19 = [r for r in rates_19 if r > 0]

    if valid_20 and valid_19:
        avg_20 = np.mean(valid_20)
        avg_19 = np.mean(valid_19)
        improvement = (avg_19 - avg_20) / avg_20 * 100 if avg_20 > 0 else 0
        logger.info("")
        logger.info("Average scaffold detection rate:")
        logger.info("  20 bp (our analysis): %.2f%%", avg_20)
        logger.info("  19 bp (company):      %.2f%%", avg_19)
        logger.info("  Relative difference:  %+.1f%%", improvement)

        if improvement > 5:
            logger.info("")
            logger.info("RECOMMENDATION: The 19 bp scaffold detects ~%.0f%% more reads.", improvement)
            logger.info("Consider switching for de novo extraction if mapping rates remain low.")
        else:
            logger.info("")
            logger.info("Both scaffolds perform similarly. 20 bp recommended for consistency.")

    meta = {
        "scaffold_20bp": SCAFFOLD_20BP,
        "scaffold_19bp": SCAFFOLD_19BP,
        "n_reads_sampled": args.n_reads,
        "n_test_samples": len(TEST_SAMPLES),
        "avg_rate_20bp": round(float(np.mean(valid_20)), 2) if valid_20 else None,
        "avg_rate_19bp": round(float(np.mean(valid_19)), 2) if valid_19 else None,
    }
    export_metadata(args.outdir, step="08_scaffold_comparison", extra=meta, config=cfg)
    logger.info("\nDone.")


if __name__ == "__main__":
    main()
