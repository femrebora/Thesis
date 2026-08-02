#!/usr/bin/env python3
"""
03_mageck_test.py — MAGeCK RRA Statistical Testing
====================================================

Runs MAGeCK test (Robust Rank Aggregation, RRA) for three comparisons:

  1. Spost_vs_S0    — Sensitive post-TRAIL vs Day-0 baseline
  2. Rpost_vs_R0    — Resistant post-TRAIL vs Day-0 baseline
  3. Spost_vs_Rpost — Sensitive vs Resistant (both post-TRAIL)

The RRA algorithm aggregates sgRNA-level fold changes to produce
gene-level statistics, accounting for the variable number of sgRNAs
per gene. Two-sided p-values are derived from the minimum of the
negative- and positive-selection p-values.

Output (per comparison):
    {comparison}.gene_summary.txt   — gene-level statistics
    {comparison}.sgrna_summary.txt  — sgRNA-level statistics
    {comparison}.log                — MAGeCK run log

Project: PRJ25-131 | Cingöz Lab | 2026

Usage:
    python 03_mageck_test.py
    python 03_mageck_test.py --count-file count/count_table_fixed.txt
    python 03_mageck_test.py --gene-lfc-method mean --remove-zero none
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


from _thesis_utils import (
    COMP_GROUPS,
    COMPARISONS,
    ThesisConfig,
    export_metadata,
    load_gene_summary,
    run_shell,
    setup_logging,
)

logger = setup_logging("03_mageck_test")


def main() -> None:
    cfg = ThesisConfig()

    parser = argparse.ArgumentParser(
        description="Run MAGeCK RRA test for 3 CRISPR screen comparisons.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--count-file", type=Path,
        default=Path(__file__).resolve().parent / "count" / "count_table_fixed.txt",
        help="Path to fixed count table (from 02_mageck_count.py).",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(__file__).resolve().parent / "test",
        help="Output directory for test results.",
    )
    parser.add_argument(
        "--gene-lfc-method", type=str, default=cfg.gene_lfc_method,
        choices=["median", "mean"],
    )
    parser.add_argument(
        "--remove-zero", type=str, default=cfg.remove_zero,
        choices=["both", "control", "treatment", "none"],
    )
    parser.add_argument(
        "--remove-zero-threshold", type=int, default=cfg.remove_zero_threshold,
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip comparisons where output already exists.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.count_file.exists():
        logger.error("Count table not found: %s", args.count_file)
        logger.error("Run 02_mageck_count.py first.")
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("STEP 3 — MAGeCK RRA Test")
    logger.info("=" * 60)
    logger.info("Count table:      %s", args.count_file)
    logger.info("Output:           %s", args.outdir)
    logger.info("Comparisons:      %d", len(COMPARISONS))
    logger.info("Gene LFC method:  %s", args.gene_lfc_method)
    logger.info("Remove-zero:      %s (threshold=%d)",
                args.remove_zero, args.remove_zero_threshold)
    logger.info("")

    results_summary: dict[str, dict] = {}

    for comp_name, (treatment, control) in COMP_GROUPS.items():
        description = COMPARISONS[comp_name]
        out_prefix = args.outdir / comp_name

        if args.skip_existing and out_prefix.with_suffix(".gene_summary.txt").exists():
            logger.info("  [skip] %s — already exists", comp_name)
            df = load_gene_summary(out_prefix.with_suffix(".gene_summary.txt"))
            if df is not None:
                results_summary[comp_name] = {
                    "n_genes": len(df),
                    "n_sig": int(df["sig"].sum()),
                    "min_pval": float(df["pval"].min()),
                }
            continue

        logger.info("─" * 50)
        logger.info("Comparison: %s", comp_name)
        logger.info("  %s", description)
        logger.info("  Treatment: %s", treatment)
        logger.info("  Control:   %s", control)

        cmd = (
            f"mageck test"
            f" -k '{args.count_file}'"
            f" -t {treatment}"
            f" -c {control}"
            f" -n '{out_prefix}'"
            f" --gene-lfc-method {args.gene_lfc_method}"
            f" --remove-zero {args.remove_zero}"
            f" --remove-zero-threshold {args.remove_zero_threshold}"
        )

        ret = run_shell(cmd, timeout=900)
        if ret != 0:
            logger.error("MAGeCK test failed for %s", comp_name)
            continue

        # Summarise results
        df = load_gene_summary(out_prefix.with_suffix(".gene_summary.txt"))
        if df is not None:
            n_sig = int(df["sig"].sum())
            min_p = float(df["pval"].min())
            max_abs_lfc = float(max(abs(df["LFC"].max()), abs(df["LFC"].min())))
            results_summary[comp_name] = {
                "n_genes": len(df),
                "n_sig": n_sig,
                "min_pval": min_p,
                "max_abs_lfc": max_abs_lfc,
            }
            logger.info("  → %d testable genes, %d significant (p<%.2f), "
                        "min p=%.4f, max |LFC|=%.2f",
                        len(df), n_sig, cfg.pval_threshold, min_p, max_abs_lfc)

    # ── Print summary table ─────────────────────────────────────────────────
    if results_summary:
        logger.info("")
        logger.info("=" * 60)
        logger.info("SUMMARY")
        logger.info("=" * 60)
        header = f"{'Comparison':<22} {'Genes':>6} {'p<0.05':>7} {'Min p':>9} {'Max|LFC|':>9}"
        logger.info(header)
        logger.info("-" * len(header))
        for comp_name, info in results_summary.items():
            logger.info(
                "%-22s %6d %7d %9.4f %9.2f",
                comp_name, info["n_genes"], info["n_sig"],
                info["min_pval"], info["max_abs_lfc"],
            )

    # ── Export metadata ─────────────────────────────────────────────────────
    meta = {
        "count_file": str(args.count_file),
        "gene_lfc_method": args.gene_lfc_method,
        "remove_zero": args.remove_zero,
        "remove_zero_threshold": args.remove_zero_threshold,
        "pval_threshold": cfg.pval_threshold,
        "comparisons": results_summary,
    }
    export_metadata(args.outdir, step="03_mageck_test", extra=meta, config=cfg)

    logger.info("")
    logger.info("Output files in %s:", args.outdir)
    for f in sorted(args.outdir.iterdir()):
        if f.is_file():
            logger.info("  %s", f.name)
    logger.info("Done.")


if __name__ == "__main__":
    main()
