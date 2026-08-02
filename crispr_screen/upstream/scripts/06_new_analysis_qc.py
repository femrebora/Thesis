#!/usr/bin/env python3
"""
06_new_analysis_qc.py — Deep QC Inspection (New Analysis)
===========================================================

Compares nf-core (R1 only) vs R1+R2 combined MAGeCK runs with
detailed quality-control figures and a text comparison report.

Figures generated:
  QC3 — Testable gene count comparison (bar chart)
  QC4 — LFC & p-value scatter between runs (×3 comparisons)

Output:
    results/figures/QC{3,4}_*.{pdf,png}
    results/reports/qc_comparison_report.txt
    results/reports/metadata_06_new_analysis_qc.json

Project: PRJ25-131 | Cingöz Lab | 2026

Usage:
    python 06_new_analysis_qc.py
    python 06_new_analysis_qc.py --mt-dir /path/to/nfcore --lt-dir /path/to/combined
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

from _thesis_utils import (
    COMPARISONS, ThesisConfig, ThesisStyle,
    export_metadata, load_gene_summary,
    save_figure, setup_logging,
)

logger = setup_logging("06_new_analysis_qc")
warnings.filterwarnings("ignore")

COMP_LABELS = {
    "Spost_vs_S0":    "Sensitive: Post vs Day-0",
    "Rpost_vs_R0":    "Resistant: Post vs Day-0",
    "Spost_vs_Rpost": "Sensitive vs Resistant (Post)",
}


def qc3_testable_genes(
    results: dict,
    outdir: Path,
) -> None:
    """Bar chart comparing testable gene counts: nf-core vs R1+R2 combined."""
    logger.info("QC3 — Testable genes comparison")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = np.arange(len(COMPARISONS))
    w = 0.32

    nf_counts = [
        len(results[c]["nfcore"]) if results[c]["nfcore"] is not None else 0
        for c in COMPARISONS
    ]
    comb_counts = [
        len(results[c]["combined"]) if results[c]["combined"] is not None else 0
        for c in COMPARISONS
    ]

    ax.bar(x - w/2, nf_counts, width=w, color="#2980B9",
           label="nf-core (R1 only)", edgecolor="none")
    ax.bar(x + w/2, comb_counts, width=w, color="#E67E22",
           label="R1+R2 Combined", edgecolor="none")

    for xi, (nf, cb) in enumerate(zip(nf_counts, comb_counts)):
        ax.text(xi - w/2, nf + 1, str(nf), ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color="#2980B9")
        ax.text(xi + w/2, cb + 1, str(cb), ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color="#E67E22")

    ax.set_xticks(x)
    ax.set_xticklabels([COMP_LABELS[c] for c in COMPARISONS], fontsize=8.5)
    ax.set_ylabel("Testable Genes (MAGeCK, --remove-zero both)")
    ax.set_title("Testable Gene Count: nf-core vs R1+R2 Combined\n"
                 "R1+R2 combination reduced testable genes in Spost_vs_S0 "
                 "(due to Spost_1 outlier)",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=8)
    ax.axhline(50, color="#888", lw=0.7, ls=":", label="Min for reliable FDR")
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, "QC3_testable_genes_comparison", outdir)


def qc4_lfc_comparison(
    results: dict,
    outdir: Path,
) -> None:
    """LFC & p-value scatter plots comparing nf-core vs combined runs."""
    for cname in COMPARISONS:
        nf = results[cname].get("nfcore")
        comb = results[cname].get("combined")
        if nf is None or comb is None:
            continue

        merged = (
            nf.set_index("id")[["LFC", "pval"]]
            .rename(columns={"LFC": "LFC_nf", "pval": "pval_nf"})
            .join(
                comb.set_index("id")[["LFC", "pval"]]
                .rename(columns={"LFC": "LFC_cb", "pval": "pval_cb"}),
                how="inner",
            )
            .dropna()
        )

        if len(merged) < 3:
            logger.warning("  [SKIP] %s: too few common genes (%d)", cname, len(merged))
            continue

        logger.info("QC4 — LFC comparison: %s (%d shared genes)", cname, len(merged))

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        fig.suptitle(
            f"LFC & p-value: nf-core vs R1+R2 Combined\n{COMP_LABELS[cname]}",
            fontsize=11, fontweight="bold",
        )

        # LFC scatter
        ax = axes[0]
        sig_either = (merged["pval_nf"] < 0.05) | (merged["pval_cb"] < 0.05)
        ax.scatter(merged.loc[~sig_either, "LFC_nf"],
                   merged.loc[~sig_either, "LFC_cb"],
                   s=12, c="#AAA", alpha=0.5, rasterized=True)
        ax.scatter(merged.loc[sig_either, "LFC_nf"],
                   merged.loc[sig_either, "LFC_cb"],
                   s=45, c="#E74C3C", alpha=0.85, edgecolors="white", linewidths=0.3)
        lim = max(abs(merged[["LFC_nf", "LFC_cb"]].values).max(), 1) * 1.1
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.7, alpha=0.45)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("LFC — nf-core run")
        ax.set_ylabel("LFC — R1+R2 combined run")
        ax.set_title("Log₂ Fold Change")
        rho, _ = spearmanr(merged["LFC_nf"], merged["LFC_cb"])
        ax.text(0.04, 0.95, f"ρ = {rho:.3f}", transform=ax.transAxes,
                fontsize=10, fontweight="bold")
        for gene, row in merged[sig_either].iterrows():
            ax.annotate(gene, (row["LFC_nf"], row["LFC_cb"]),
                        fontsize=6.5, xytext=(3, 3), textcoords="offset points",
                        color="#C0392B")
        sns.despine(ax=ax)

        # p-value scatter
        ax = axes[1]
        ax.scatter(merged.loc[~sig_either, "pval_nf"],
                   merged.loc[~sig_either, "pval_cb"],
                   s=12, c="#AAA", alpha=0.5, rasterized=True)
        ax.scatter(merged.loc[sig_either, "pval_nf"],
                   merged.loc[sig_either, "pval_cb"],
                   s=45, c="#E74C3C", alpha=0.85, edgecolors="white", linewidths=0.3)
        ax.axhline(0.05, color="#E74C3C", lw=0.7, ls="--", alpha=0.5)
        ax.axvline(0.05, color="#E74C3C", lw=0.7, ls="--", alpha=0.5)
        ax.set_xlabel("p-value — nf-core run")
        ax.set_ylabel("p-value — R1+R2 combined run")
        ax.set_title("p-values (min of neg/pos direction)")
        for gene, row in merged[sig_either].iterrows():
            ax.annotate(gene, (row["pval_nf"], row["pval_cb"]),
                        fontsize=6.5, xytext=(3, 3), textcoords="offset points",
                        color="#C0392B")
        sns.despine(ax=ax)

        fig.tight_layout()
        save_figure(fig, f"QC4_lfc_comparison_{cname}", outdir)


def main() -> None:
    cfg = ThesisConfig()

    parser = argparse.ArgumentParser(
        description="Deep QC inspection comparing nf-core vs R1+R2 combined runs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--mt-dir", type=Path,
                        default=cfg.crispr_root / "analysis" / "mageck_test")
    parser.add_argument("--lt-dir", type=Path,
                        default=cfg.crispr_root / "analysis" / "results" / "last_test")
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).resolve().parent / "figures")
    parser.add_argument("--report-dir", type=Path,
                        default=Path(__file__).resolve().parent / "reports")
    args = parser.parse_args()

    ThesisStyle.apply(font_family=cfg.font_family, dpi=cfg.figure_dpi)
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("STEP 6 — Deep QC Inspection (New Analysis)")
    logger.info("=" * 60)
    logger.info("nf-core results:  %s", args.mt_dir)
    logger.info("Combined results: %s", args.lt_dir)
    logger.info("")

    # ── Load ────────────────────────────────────────────────────────────────
    results: dict[str, dict[str, pd.DataFrame | None]] = {}
    for cname in COMPARISONS:
        nf = load_gene_summary(args.mt_dir / f"{cname}.gene_summary.txt")
        comb = load_gene_summary(args.lt_dir / f"{cname}.gene_summary.txt")
        results[cname] = {"nfcore": nf, "combined": comb}
        n_nf = len(nf) if nf is not None else 0
        n_comb = len(comb) if comb is not None else 0
        sig_nf = (nf["sig"].sum()) if nf is not None else 0
        sig_comb = (comb["sig"].sum()) if comb is not None else 0
        logger.info("%-22s: nf-core=%3d genes (%d sig), combined=%3d genes (%d sig)",
                    cname, n_nf, sig_nf, n_comb, sig_comb)

    # ── QC3: Testable genes comparison ──────────────────────────────────────
    qc3_testable_genes(results, args.outdir)

    # ── QC4: LFC comparison ─────────────────────────────────────────────────
    qc4_lfc_comparison(results, args.outdir)

    # ── Cross-run gene comparison report ────────────────────────────────────
    report_path = args.report_dir / "qc_comparison_report.txt"
    with open(report_path, "w") as f:
        f.write("CRISPR Screen QC Comparison Report\n")
        f.write("=" * 60 + "\n")
        f.write("nf-core (R1 only) vs R1+R2 Combined Run\n\n")

        for cname in COMPARISONS:
            nf = results[cname]["nfcore"]
            comb = results[cname]["combined"]
            nf_sig = set(nf[nf["sig"]]["id"]) if nf is not None else set()
            comb_sig = set(comb[comb["sig"]]["id"]) if comb is not None else set()
            common = nf_sig & comb_sig
            nf_only = nf_sig - comb_sig
            comb_only = comb_sig - nf_sig

            f.write(f"\n{cname} ({COMP_LABELS[cname]}):\n")
            f.write(f"  Both runs ({len(common)}):     {sorted(common)}\n")
            f.write(f"  nf-core only ({len(nf_only)}):  {sorted(nf_only)}\n")
            f.write(f"  Combined only ({len(comb_only)}): {sorted(comb_only)}\n")

            logger.info("")
            logger.info("%s:", cname)
            logger.info("  Both runs (%d):     %s", len(common), sorted(common)[:8])
            logger.info("  nf-core only (%d):  %s", len(nf_only), sorted(nf_only)[:8])
            logger.info("  Combined only (%d): %s", len(comb_only), sorted(comb_only)[:8])

    logger.info("\n[OK] QC comparison report: %s", report_path)

    meta = {
        "nfcore_dir": str(args.mt_dir),
        "combined_dir": str(args.lt_dir),
        "comparison_summary": {
            cname: {
                "nfcore_genes": len(results[cname]["nfcore"]) if results[cname]["nfcore"] is not None else 0,
                "combined_genes": len(results[cname]["combined"]) if results[cname]["combined"] is not None else 0,
            }
            for cname in COMPARISONS
        },
    }
    export_metadata(args.report_dir, step="06_new_analysis_qc", extra=meta, config=cfg)
    logger.info("Done.")


if __name__ == "__main__":
    main()
