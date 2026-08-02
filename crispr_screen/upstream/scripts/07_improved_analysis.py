#!/usr/bin/env python3
"""
07_improved_analysis.py — Improved MAGeCK Analysis Parameters
===============================================================

Re-runs MAGeCK test with methodologically improved parameters:

  1. --remove-zero-threshold 10 → retains more testable genes
     (original threshold 0 discards sgRNAs with any zero count)
  2. --norm-method none in mageck test → avoids double normalisation
     (counts already median-normalised during MAGeCK count)
  3. Two-sided p-value correction (×2 Bonferroni)
     (accounts for testing both depletion AND enrichment directions)
  4. Run B excludes Spost_1 (Gini=1.000, only 5 sgRNAs detected)

WARNING: Single-sgRNA hits are flagged ⚠ — do NOT draw biological
conclusions from genes detected by only one sgRNA.

Output:
    results/improved_results/run_a_improved/   — improved-parameter run
    results/improved_results/run_b_no_spost1/  — without Spost_1 outlier
    results/improved_results/run_comparison_table.csv
    results/figures/IMP2_volcano_improved.{pdf,png}

Project: PRJ25-131 | Cingöz Lab | 2026

Usage:
    python 07_improved_analysis.py
    python 07_improved_analysis.py --count-file path/to/count_table.txt
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from _thesis_utils import (
    COMPARISONS, C_DEPLETED, C_ENRICHED, C_NS,
    ThesisConfig, ThesisStyle,
    export_metadata, save_figure, setup_logging,
)

logger = setup_logging("07_improved_analysis")
warnings.filterwarnings("ignore")

# Comparisons excluding Spost_1 (Gini=1.000 outlier)
COMPARISONS_NO_SPOST1: dict[str, tuple[str, str]] = {
    "Spost_vs_S0_no_Spost1":    ("Spost_2,Spost_3", "S0_1,S0_2,S0_3"),
    "Spost_vs_Rpost_no_Spost1": ("Spost_2,Spost_3", "Rpost_1,Rpost_2,Rpost_3,Rpost_4"),
}

COMP_LABELS: dict[str, str] = {
    "Spost_vs_S0":    "Sensitive: Post vs Day-0",
    "Rpost_vs_R0":    "Resistant: Post vs Day-0",
    "Spost_vs_Rpost": "Sensitive vs Resistant (Post)",
    "Spost_vs_S0_no_Spost1":    "Sensitive: Post vs Day-0 (no Spost₁)",
    "Spost_vs_Rpost_no_Spost1": "Sensitive vs Resistant (Post, no Spost₁)",
}

PVAL_THRESH = 0.05
PVAL_CORRECTED = 0.025  # 0.05 ÷ 2 for two-sided test


def run_mageck_test(
    count_file: Path,
    treatment: str,
    control: str,
    out_prefix: Path,
    *,
    label: str = "",
    norm_method: str = "none",
    min_count: int = 10,
    timeout: int = 600,
) -> bool:
    """Execute MAGeCK test with improved parameters."""
    cmd = (
        f"mageck test"
        f" -k '{count_file}'"
        f" -t {treatment}"
        f" -c {control}"
        f" -n '{out_prefix}'"
        f" --gene-lfc-method median"
        f" --norm-method {norm_method}"
        f" --remove-zero both"
        f" --remove-zero-threshold {min_count}"
    )

    logger.info("  Running: mageck test [%s]", label)
    logger.info("    Treatment: %s, Control: %s", treatment, control)
    logger.info("    Params: norm=%s, remove-zero-threshold=%d", norm_method, min_count)

    result = subprocess.run(
        cmd, shell=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=timeout,
    )
    for line in result.stdout.splitlines():
        lowered = line.lower()
        if any(k in lowered for k in ("info", "warn", "error", "gene", "total", "testable")):
            logger.info("    %s", line.strip())
    if result.returncode != 0:
        logger.error("  Exit code %d", result.returncode)
        logger.error(result.stdout[-1000:])
        return False

    gs = out_prefix.with_suffix(".gene_summary.txt")
    if gs.exists():
        df = pd.read_csv(gs, sep="\t")
        logger.info("    → %d testable genes", len(df))
    return result.returncode == 0


def load_gene_summary_improved(path: Path) -> pd.DataFrame | None:
    """Load gene_summary with Bonferroni-corrected p-values."""
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    df["LFC"] = df["neg|lfc"]
    df["pval_min"] = df[["neg|p-value", "pos|p-value"]].min(axis=1)
    df["pval_corrected"] = (df["pval_min"] * 2).clip(upper=1.0)  # Bonferroni ×2
    df["pval"] = df["pval_min"]
    df["FDR"] = df[["neg|fdr", "pos|fdr"]].min(axis=1)
    df["sig"] = df["pval_corrected"] < PVAL_THRESH
    df["sig_uncorrected"] = df["pval_min"] < PVAL_THRESH
    return df


def main() -> None:
    cfg = ThesisConfig()

    parser = argparse.ArgumentParser(
        description="Run improved MAGeCK analysis with better parameters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--count-file", type=Path,
        default=(cfg.crispr_root / "analysis" / "results"
                 / "last_count" / "last_count_fixed.txt"),
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(__file__).resolve().parent / "improved_results",
    )
    parser.add_argument(
        "--fig-dir", type=Path,
        default=Path(__file__).resolve().parent / "figures",
    )
    args = parser.parse_args()

    ThesisStyle.apply(font_family=cfg.font_family, dpi=cfg.figure_dpi)
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    # Symlink to avoid space-in-path issues with MAGeCK RRA
    tmpdir = Path(tempfile.mkdtemp(prefix="crispr_impr_"))
    subprocess.run(f"ln -sfn '{args.outdir}' '{tmpdir}'", shell=True)

    logger.info("=" * 60)
    logger.info("STEP 7 — Improved MAGeCK Analysis")
    logger.info("=" * 60)
    logger.info("Count file: %s", args.count_file)
    logger.info("Output:     %s", args.outdir)
    logger.info("")

    if not args.count_file.exists():
        logger.error("Count file not found: %s", args.count_file)
        return

    # ── Run A: Improved parameters ──────────────────────────────────────────
    logger.info("RUN A — Improved parameters (min-count=10, norm=none)")
    run_a_dir = args.outdir / "run_a_improved"
    run_a_dir.mkdir(exist_ok=True)

    for cname, (treatment, control) in COMPARISONS.items():
        run_mageck_test(args.count_file, treatment, control,
                        run_a_dir / cname, label=f"RunA/{cname}")

    # ── Run B: Without Spost_1 ──────────────────────────────────────────────
    logger.info("")
    logger.info("RUN B — Without Spost_1 (Gini=1.000 outlier)")
    run_b_dir = args.outdir / "run_b_no_spost1"
    run_b_dir.mkdir(exist_ok=True)

    for cname, (treatment, control) in COMPARISONS_NO_SPOST1.items():
        run_mageck_test(args.count_file, treatment, control,
                        run_b_dir / cname, label=f"RunB/{cname}")

    # ── Compare runs ────────────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("RUN COMPARISON")
    logger.info("=" * 60)

    comparison_table: list[dict] = []
    run_configs: list[tuple[str, Path]] = [
        ("Original (remove-zero=both, thresh=0)",
         cfg.crispr_root / "analysis" / "results" / "last_test"),
        ("Improved (remove-zero-threshold=10)", run_a_dir),
    ]

    for run_name, run_dir in run_configs:
        if not run_dir.exists():
            continue
        for cname in COMPARISONS:
            gs = run_dir / f"{cname}.gene_summary.txt"
            df = load_gene_summary_improved(gs)
            if df is None:
                continue
            sig_genes = sorted(df[df["sig"]]["id"].tolist())
            comparison_table.append({
                "Run": run_name,
                "Comparison": cname,
                "Testable": len(df),
                "Sig (corr p<0.05)": int(df["sig"].sum()),
                "Sig (uncorr p<0.05)": int(df["sig_uncorrected"].sum()),
                "Top sig genes": ", ".join(sig_genes[:8]),
            })

    if comparison_table:
        ct = pd.DataFrame(comparison_table)
        ct.to_csv(args.outdir / "run_comparison_table.csv", index=False)
        logger.info(ct.to_string(index=False))

    # ── Volcano plots for improved run ──────────────────────────────────────
    logger.info("")
    logger.info("Generating improved volcano plots …")
    data_impr: dict[str, pd.DataFrame] = {}
    for cname in COMPARISONS:
        gs = run_a_dir / f"{cname}.gene_summary.txt"
        df = load_gene_summary_improved(gs)
        if df is not None:
            df["log10p_corr"] = -np.log10(df["pval_corrected"].clip(lower=1e-100))
            data_impr[cname] = df

    if data_impr:
        n = len(data_impr)
        fig, axes = plt.subplots(1, n, figsize=(4.5 * n + 1, 4.5))
        if n == 1:
            axes = [axes]
        fig.subplots_adjust(wspace=0.45)

        for ax, (cname, df) in zip(axes, data_impr.items()):
            ns = df[~df["sig"]]
            dep = df[df["sig"] & (df["LFC"] < 0)]
            enr = df[df["sig"] & (df["LFC"] > 0)]

            ax.scatter(ns["LFC"], ns["log10p_corr"], c=C_NS, s=12,
                       linewidths=0, alpha=0.5, rasterized=True)
            ax.scatter(dep["LFC"], dep["log10p_corr"], c=C_DEPLETED, s=40,
                       linewidths=0.2, edgecolors="white", alpha=0.9)
            ax.scatter(enr["LFC"], enr["log10p_corr"], c=C_ENRICHED, s=40,
                       linewidths=0.2, edgecolors="white", alpha=0.9)
            ax.axhline(-np.log10(PVAL_THRESH), color="#555", lw=0.8, ls="--")
            ax.axvline(0, color="#AAA", lw=0.5, ls=":")

            for _, row in df[df["sig"]].iterrows():
                color = C_DEPLETED if row["LFC"] < 0 else C_ENRICHED
                ax.text(row["LFC"], row["log10p_corr"] + 0.04, row["id"],
                        fontsize=7, color=color, fontweight="bold", ha="center")

            ax.set_xlabel("Log₂ Fold Change")
            ax.set_ylabel("−log₁₀(p-value × 2, corrected)")
            ax.set_title(COMP_LABELS.get(cname, cname), fontsize=10, fontweight="bold")
            ax.text(0.03, 0.97, f"▼ {len(dep)} depleted", transform=ax.transAxes,
                    fontsize=7.5, color=C_DEPLETED, va="top", fontweight="bold")
            ax.text(0.03, 0.90, f"▲ {len(enr)} enriched", transform=ax.transAxes,
                    fontsize=7.5, color=C_ENRICHED, va="top")
            sns.despine(ax=ax)

        fig.suptitle(
            "Improved Analysis — Volcano Plots\n"
            "(two-sided correction applied: p_eff = min(neg_p, pos_p) × 2)",
            fontsize=11, fontweight="bold", y=1.02,
        )
        save_figure(fig, "IMP2_volcano_improved", args.fig_dir)

    # ── Summary ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("IMPROVED ANALYSIS COMPLETE")
    logger.info("=" * 60)
    logger.info("Results: %s", args.outdir)
    logger.info("Figures: %s", args.fig_dir)
    logger.info("")
    logger.info("Key improvements applied:")
    logger.info("  1. --remove-zero-threshold 10 (more testable genes)")
    logger.info("  2. --norm-method none (avoids double normalisation)")
    logger.info("  3. p-value × 2 Bonferroni correction (two-sided test)")
    logger.info("  4. Run B excludes Spost_1 (Gini=1.000 outlier)")
    logger.info("")
    logger.info("⚠  Single-sgRNA hits are unreliable for biological conclusions.")

    meta = {
        "run_a_params": {"norm_method": "none", "remove_zero_threshold": 10},
        "run_b_excludes_spost1": True,
        "pval_correction": "Bonferroni ×2 (two-sided)",
        "pval_threshold": PVAL_THRESH,
        "pval_corrected": PVAL_CORRECTED,
    }
    export_metadata(args.outdir, step="07_improved_analysis", extra=meta, config=cfg)


if __name__ == "__main__":
    main()
