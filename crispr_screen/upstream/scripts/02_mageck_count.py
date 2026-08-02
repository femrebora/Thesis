#!/usr/bin/env python3
"""
02_mageck_count.py — MAGeCK sgRNA Counting
============================================

Runs MAGeCK count on combined R1+R2 FASTQ files against the Sabatini
human metabolic gene knockout library (Addgene #110066).

Performs:
  1. sgRNA sequence alignment against the reference library
  2. Read counting per sgRNA per sample
  3. Median normalisation
  4. Automatic duplicate column-name fix (nf-core MAGeCK artefact)

Input:
    results/combined_fastq/{label}_combined.fq.gz
    nf-core_crisprseq/library_sabatini_metko.tsv

Output:
    results/count/count_table.count.txt              — raw counts
    results/count/count_table.count_normalized.txt   — normalised counts
    results/count/count_table.countsummary.txt       — per-sample statistics
    results/count/count_table_fixed.txt               — deduplicated column names
    results/count/metadata_02_mageck_count.json

Project: PRJ25-131 | Cingöz Lab | 2026

Usage:
    python 02_mageck_count.py
    python 02_mageck_count.py --scaffold TTTTAGAGCTAGAAATAGC --sgrna-len 20
    python 02_mageck_count.py --norm-method total --no-fix-columns
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from _thesis_utils import (
    SAMPLE_MAP,
    ThesisConfig,
    export_metadata,
    run_shell,
    setup_logging,
)

logger = setup_logging("02_mageck_count")


def fix_duplicate_columns(raw_count: Path, out_path: Path) -> Path:
    """
    Fix duplicate column names in MAGeCK count output.

    nf-core MAGeCK labels all sample columns identically (a known bug).
    This function appends ``_{n}`` suffix to duplicate column headers.

    Parameters
    ----------
    raw_count : Path
        Path to the raw MAGeCK count table.
    out_path : Path
        Output path for the fixed table.

    Returns
    -------
    Path
        Path to the fixed count table.

    Raises
    ------
    FileNotFoundError
        If raw_count does not exist.
    """
    if not raw_count.exists():
        raise FileNotFoundError(f"Count table not found: {raw_count}")

    df = pd.read_csv(raw_count, sep="\t")
    cols: list[str] = list(df.columns)
    seen: dict[str, int] = {}
    new_cols: list[str] = []

    for c in cols:
        if c in ("sgRNA", "Gene"):
            new_cols.append(c)
        else:
            n = seen.get(c, 0)
            if n > 0:
                new_cols.append(f"{c}_{n + 1}")
            else:
                new_cols.append(c)
            seen[c] = n + 1

    df.columns = new_cols
    df.to_csv(out_path, sep="\t", index=False)
    logger.info("  → Fixed count table: %s (%s sgRNAs)", out_path.name, f"{len(df):,}")
    return out_path


def main() -> None:
    cfg = ThesisConfig()

    parser = argparse.ArgumentParser(
        description="Run MAGeCK count on concatenated R1+R2 FASTQ files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--fastq-dir", type=Path,
                        default=Path(__file__).resolve().parent / "combined_fastq")
    parser.add_argument("--library", type=Path, default=cfg.library_file)
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).resolve().parent / "count")
    parser.add_argument("--scaffold", type=str, default=cfg.scaffold_seq)
    parser.add_argument("--sgrna-len", type=int, default=cfg.sgrna_len)
    parser.add_argument("--norm-method", type=str, default=cfg.norm_method,
                        choices=["median", "total", "none"])
    parser.add_argument("--no-fix-columns", action="store_true",
                        help="Skip duplicate column name fix.")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # ── Validate inputs ─────────────────────────────────────────────────────
    if not args.library.exists():
        logger.error("Library file not found: %s", args.library)
        sys.exit(1)

    # Build list of available FASTQs
    fastq_paths: list[Path] = []
    missing_labels: list[str] = []
    for record in SAMPLE_MAP:
        fq = args.fastq_dir / record.combined_name
        if fq.exists():
            fastq_paths.append(fq)
        else:
            missing_labels.append(record.label)

    if missing_labels:
        logger.warning("Missing FASTQ for %d sample(s): %s",
                       len(missing_labels), ", ".join(missing_labels))

    if not fastq_paths:
        logger.error("No combined FASTQ files found in %s", args.fastq_dir)
        logger.error("Run 01_prepare_fastq.py first.")
        sys.exit(1)

    # ── Print run info ──────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("STEP 2 — MAGeCK Count (R1+R2 Combined)")
    logger.info("=" * 60)
    logger.info("FASTQ files:     %d/%d", len(fastq_paths), len(SAMPLE_MAP))
    logger.info("Library:         %s", args.library)
    logger.info("Output:          %s", args.outdir)
    logger.info("Scaffold:        %s (%d bp)", args.scaffold, len(args.scaffold))
    logger.info("sgRNA length:    %d bp", args.sgrna_len)
    logger.info("Normalisation:   %s", args.norm_method)
    logger.info("")

    # ── Build MAGeCK count command ──────────────────────────────────────────
    count_out = args.outdir / "count_table"
    label_str = ",".join(
        fq.stem.replace("_combined", "") for fq in fastq_paths
    )
    fastq_str = " ".join(f"'{fq}'" for fq in fastq_paths)

    cmd = (
        f"mageck count"
        f" -l '{args.library}'"
        f" -n '{count_out}'"
        f" --sample-label '{label_str}'"
        f" --fastq {fastq_str}"
        f" --sgrna-len {args.sgrna_len}"
        f" --norm-method {args.norm_method}"
    )

    ret = run_shell(cmd, timeout=1200)  # 20 min timeout for large FASTQs
    if ret != 0:
        logger.error("MAGeCK count failed with exit code %d", ret)
        sys.exit(1)

    # ── Fix duplicate columns ───────────────────────────────────────────────
    if not args.no_fix_columns:
        raw_count = count_out.with_suffix(".count.txt")
        fixed_count = args.outdir / "count_table_fixed.txt"
        try:
            fix_duplicate_columns(raw_count, fixed_count)
        except FileNotFoundError:
            logger.warning("Raw count table not found; skipping column fix.")

    # ── Print count summary ─────────────────────────────────────────────────
    summary_file = count_out.with_suffix(".countsummary.txt")
    if summary_file.exists():
        summary = pd.read_csv(summary_file, sep="\t")
        logger.info("")
        logger.info("Count Summary:")
        logger.info("  Total reads:     %s", f"{summary['Reads'].sum():,}")
        logger.info("  Total mapped:    %s", f"{summary['Mapped'].sum():,}")
        logger.info("  Mean mapping %%:  %.1f%%", summary["Percentage"].mean() * 100)
        logger.info("  Mean Gini index: %.3f", summary["GiniIndex"].mean())
        logger.info("  Samples:         %d", len(summary))

    # ── Export metadata ─────────────────────────────────────────────────────
    meta = {
        "n_fastq_input": len(fastq_paths),
        "n_samples": len(SAMPLE_MAP),
        "library": str(args.library),
        "scaffold": args.scaffold,
        "sgrna_len": args.sgrna_len,
        "norm_method": args.norm_method,
    }
    if summary_file.exists():
        summary = pd.read_csv(summary_file, sep="\t")
        meta["count_summary"] = {
            "total_reads": int(summary["Reads"].sum()),
            "total_mapped": int(summary["Mapped"].sum()),
            "mean_mapping_pct": round(float(summary["Percentage"].mean()) * 100, 2),
            "mean_gini": round(float(summary["GiniIndex"].mean()), 4),
        }
    export_metadata(args.outdir, step="02_mageck_count", extra=meta, config=cfg)

    logger.info("")
    logger.info("Output files:")
    for f in sorted(args.outdir.iterdir()):
        if f.is_file():
            logger.info("  %s", f.name)
    logger.info("Done.")


if __name__ == "__main__":
    main()
