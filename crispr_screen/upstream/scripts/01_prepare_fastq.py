#!/usr/bin/env python3
"""
01_prepare_fastq.py — FASTQ Preparation & R1+R2 Concatenation
=============================================================

Concatenates paired-end R1 and R2 reads per sample to recover sgRNAs
sequenced in both orientations. R1 reads carry ~34.6% of sgRNAs in
reverse-complement orientation; R2 provides forward orientation.
Concatenating both maximises sgRNA detection.

Input:
    CRISPR-Seq_Analysis/Fastq/PRJ25-131-{sample}_L01_{id}_{read}.fq.gz

Output:
    results/combined_fastq/{label}_combined.fq.gz
    results/combined_fastq/metadata_01_prepare_fastq.json

Project: PRJ25-131 — GBM TRAIL Metabolic CRISPR/Cas9 Screen
Author:  Emre Bora, Cingöz Lab
Version: 1.0.0 — 2026-06-09

Usage:
    python 01_prepare_fastq.py
    python 01_prepare_fastq.py --fastq-dir /path/to/Fastq
    python 01_prepare_fastq.py --skip-existing --threads 4
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Optional

from _thesis_utils import (
    SAMPLE_MAP,
    ThesisConfig,
    export_metadata,
    run_shell,
    setup_logging,
)

logger = setup_logging("01_prepare_fastq")


def _check_binary(name: str) -> Optional[str]:
    """Check if a binary is available on PATH."""
    path = shutil.which(name)
    if path is None:
        logger.warning("%s not found on PATH — falling back to Python", name)
    return path


def concatenate_pair(
    r1: Path,
    r2: Path,
    output: Path,
    *,
    pigz_path: Optional[str] = None,
    threads: int = 1,
) -> int:
    """
    Concatenate paired-end FASTQ files (gzipped).

    Uses pigz for parallel decompression if available, otherwise
    falls back to zcat (or Python gzip).

    Parameters
    ----------
    r1, r2 : Path
        Input gzipped FASTQ files.
    output : Path
        Output gzipped FASTQ file.
    pigz_path : str | None
        Path to pigz binary.
    threads : int
        Number of decompression threads (with pigz).

    Returns
    -------
    int
        0 on success, non-zero on failure.
    """
    if not r1.exists():
        logger.error("Missing R1: %s", r1)
        return 1
    if not r2.exists():
        logger.error("Missing R2: %s", r2)
        return 1

    if pigz_path:
        cmd = (
            f"({pigz_path} -d -c -p {threads} '{r1}'; "
            f"{pigz_path} -d -c -p {threads} '{r2}') "
            f"| {pigz_path} -c -p {threads} > '{output}'"
        )
    else:
        cmd = f"cat '{r1}' '{r2}' > '{output}'"

    return run_shell(cmd, timeout=300)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Project: PRJ25-131 | Cingöz Lab | 2026",
    )
    parser.add_argument(
        "--fastq-dir", type=Path,
        default=ThesisConfig().fastq_dir,
        help="Directory containing raw FASTQ files.",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(__file__).resolve().parent / "combined_fastq",
        help="Output directory for combined FASTQ files.",
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip samples where combined FASTQ already exists.",
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of parallel decompression threads.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    pigz = _check_binary("pigz")

    logger.info("=" * 60)
    logger.info("STEP 1 — Concatenate R1 + R2 per Sample")
    logger.info("=" * 60)
    logger.info("FASTQ directory: %s", args.fastq_dir)
    logger.info("Output:          %s", args.outdir)
    logger.info("Samples:         %d", len(SAMPLE_MAP))
    logger.info("Decompressor:    %s", pigz or "zcat (system)")
    logger.info("")

    n_created, n_skipped, n_missing = 0, 0, 0
    sample_sizes: dict[str, float] = {}

    for record in SAMPLE_MAP:
        r1 = args.fastq_dir / record.r1_file
        r2 = args.fastq_dir / record.r2_file
        output = args.outdir / record.combined_name

        if args.skip_existing and output.exists():
            logger.info("  [skip] %-6s (%s) — already exists", record.label, record.condition)
            n_skipped += 1
            sample_sizes[record.label] = output.stat().st_size / 1e6
            continue

        if not r1.exists() or not r2.exists():
            missing = record.r1_file if not r1.exists() else record.r2_file
            logger.warning("  [WARN] Missing: %s", missing)
            n_missing += 1
            continue

        logger.info("  Combining %s (%s) …", record.label, record.condition)
        ret = concatenate_pair(r1, r2, output, pigz_path=pigz, threads=args.threads)
        if ret == 0:
            size_mb = output.stat().st_size / 1e6
            sample_sizes[record.label] = size_mb
            logger.info("    → %s  (%.0f MB)", output.name, size_mb)
            n_created += 1
        else:
            n_missing += 1

    # Summary
    total_size = sum(sample_sizes.values())
    logger.info("")
    logger.info("Summary: %d created, %d skipped, %d missing", n_created, n_skipped, n_missing)
    logger.info("Total combined size: %.0f MB across %d files", total_size, len(sample_sizes))

    # Export metadata
    meta = {
        "n_created": n_created, "n_skipped": n_skipped, "n_missing": n_missing,
        "total_size_mb": round(total_size, 1),
        "sample_sizes_mb": {k: round(v, 1) for k, v in sample_sizes.items()},
        "n_samples_total": len(SAMPLE_MAP),
        "condition_replicates": {k: len(v) for k, v in {
            "R0": ["R01", "R02", "R03"],
            "Rpost": ["R1", "R2", "R3", "R4"],
            "S0": ["S01", "S02", "S03"],
            "Spost": ["S1", "S2", "S3"],
        }.items()},
    }
    metapath = export_metadata(args.outdir, step="01_prepare_fastq", extra=meta)
    logger.info("Metadata: %s", metapath)
    logger.info("Done.")


if __name__ == "__main__":
    main()
