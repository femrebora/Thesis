#!/usr/bin/env python3
"""
_thesis_utils.py — Shared Thesis Utilities for CRISPR-Seq Pipeline
===================================================================

Provides standardized configuration, styling, logging, and metadata
for all pipeline scripts. Imported by 01–08 scripts.

Project:    PRJ25-131 — GBM TRAIL Metabolic CRISPR/Cas9 Screen
Lab:        Cingöz Lab
Thesis:     "Metabolic Determinants of TRAIL Resistance in GBM"
Author:     Emre Bora
Version:    1.0.0 — 2026-06-09

Usage:
    from _thesis_utils import (
        ThesisConfig, ThesisStyle, setup_logging,
        save_figure, export_metadata, COMPARISONS, COND_COLORS,
        C_DEPLETED, C_ENRICHED, C_NS, SAMPLE_MAP, build_sample_groups,
    )
"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ═══════════════════════════════════════════════════════════════════════════════
# Project Metadata
# ═══════════════════════════════════════════════════════════════════════════════

PROJECT_ID    = "PRJ25-131"
PROJECT_TITLE = "GBM TRAIL Metabolic CRISPR/Cas9 Screen"
LAB           = "Cingöz Lab"
AUTHOR        = "Emre Bora"
THESIS_TITLE  = "Metabolic Determinants of TRAIL Resistance in Glioblastoma"
LIBRARY       = "Sabatini Human Metabolic Gene Knockout Library (Addgene #110066)"
LIBRARY_STATS = "29,790 sgRNAs targeting 2,981 metabolic genes"
SEQUENCER     = "Illumina NovaSeq 6000 (PE 2×150 bp)"
PIPELINE_VERSION = "1.0.0"

# ═══════════════════════════════════════════════════════════════════════════════
# Sample Mapping
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class SampleRecord:
    """Immutable sample metadata."""
    label: str
    condition: str
    r1_file: str
    r2_file: str

    @property
    def combined_name(self) -> str:
        return f"{self.label}_combined.fq.gz"


SAMPLE_MAP: list[SampleRecord] = [
    # Resistant Day-0 (3 biological replicates)
    SampleRecord("R01", "R0",    "PRJ25-131-R01_L01_59_1.fq.gz", "PRJ25-131-R01_L01_59_2.fq.gz"),
    SampleRecord("R02", "R0",    "PRJ25-131-R02_L01_60_1.fq.gz", "PRJ25-131-R02_L01_60_2.fq.gz"),
    SampleRecord("R03", "R0",    "PRJ25-131-R03_L01_61_1.fq.gz", "PRJ25-131-R03_L01_61_2.fq.gz"),
    # Resistant Post-TRAIL (4 biological replicates)
    SampleRecord("R1",  "Rpost", "PRJ25-131-R1_L01_65_1.fq.gz",  "PRJ25-131-R1_L01_65_2.fq.gz"),
    SampleRecord("R2",  "Rpost", "PRJ25-131-R2_L01_66_1.fq.gz",  "PRJ25-131-R2_L01_66_2.fq.gz"),
    SampleRecord("R3",  "Rpost", "PRJ25-131-R3_L01_67_1.fq.gz",  "PRJ25-131-R3_L01_67_2.fq.gz"),
    SampleRecord("R4",  "Rpost", "PRJ25-131-R4_L01_68_1.fq.gz",  "PRJ25-131-R4_L01_68_2.fq.gz"),
    # Sensitive Day-0 (3 biological replicates)
    SampleRecord("S01", "S0",    "PRJ25-131-S01_L01_56_1.fq.gz", "PRJ25-131-S01_L01_56_2.fq.gz"),
    SampleRecord("S02", "S0",    "PRJ25-131-S02_L01_57_1.fq.gz", "PRJ25-131-S02_L01_57_2.fq.gz"),
    SampleRecord("S03", "S0",    "PRJ25-131-S03_L01_58_1.fq.gz", "PRJ25-131-S03_L01_58_2.fq.gz"),
    # Sensitive Post-TRAIL (3 biological replicates)
    SampleRecord("S1",  "Spost", "PRJ25-131-S1_L01_62_1.fq.gz",  "PRJ25-131-S1_L01_62_2.fq.gz"),
    SampleRecord("S2",  "Spost", "PRJ25-131-S2_L01_63_1.fq.gz",  "PRJ25-131-S2_L01_63_2.fq.gz"),
    SampleRecord("S3",  "Spost", "PRJ25-131-S3_L01_64_1.fq.gz",  "PRJ25-131-S3_L01_64_2.fq.gz"),
]

# Condition → sample labels mapping
CONDITION_SAMPLES: dict[str, list[str]] = {
    "R0":    ["R01", "R02", "R03"],
    "Rpost": ["R1", "R2", "R3", "R4"],
    "S0":    ["S01", "S02", "S03"],
    "Spost": ["S1", "S2", "S3"],
}

# ═══════════════════════════════════════════════════════════════════════════════
# Comparisons
# ═══════════════════════════════════════════════════════════════════════════════

COMPARISONS: dict[str, str] = {
    "Spost_vs_S0":    "Sensitive: Post-TRAIL vs Day-0",
    "Rpost_vs_R0":    "Resistant: Post-TRAIL vs Day-0",
    "Spost_vs_Rpost": "Sensitive vs Resistant (Post-TRAIL)",
}

COMP_GROUPS: dict[str, tuple[str, str]] = {
    "Spost_vs_S0":    ("S1,S2,S3", "S01,S02,S03"),
    "Rpost_vs_R0":    ("R1,R2,R3,R4", "R01,R02,R03"),
    "Spost_vs_Rpost": ("S1,S2,S3", "R1,R2,R3,R4"),
}

# ═══════════════════════════════════════════════════════════════════════════════
# Colour Palette — Colourblind-Friendly (Wong 2011, Nature Methods)
# ═══════════════════════════════════════════════════════════════════════════════

C_DEPLETED  = "#0072B2"   # blue — depleted (essential for TRAIL sensitivity)
C_ENRICHED  = "#D55E00"   # vermillion — enriched (confers TRAIL resistance)
C_NS        = "#B0B0B0"   # grey — not significant
C_HIGHLIGHT = "#F0E442"   # yellow — highlight

COND_COLORS: dict[str, str] = {
    "R0":    "#0072B2",   # blue
    "Rpost": "#D55E00",   # vermillion
    "S0":    "#009E73",   # bluish green
    "Spost": "#CC79A7",   # reddish purple
}

# Functional annotation for key hit genes
FUNC_ANNOT: dict[str, str] = {
    "CYP2E1":   "CYP450 2E1 — xenobiotic metabolism, ROS generation",
    "NAT2":     "N-acetyltransferase 2 — drug/xenobiotic acetylation",
    "ASL":      "Argininosuccinate lyase — urea cycle, arginine biosynthesis",
    "GBA3":     "Cytosolic β-glucosidase — glycolipid catabolism",
    "TRPM7":    "TRPM7 channel — Mg²⁺/Ca²⁺ homeostasis, cell proliferation",
    "ATP6V0A4": "V-ATPase subunit a4 — lysosomal/endosomal acidification",
    "TPCN1":    "Two-pore channel 1 — lysosomal Ca²⁺ release, autophagy",
    "HSD11B1":  "11β-HSD1 — glucocorticoid activation (cortisol/corticosterone)",
    "CYB5A":    "Cytochrome b5A — electron transfer, lipid metabolism",
    "ALOXE3":   "Epidermal lipoxygenase-3 — ceramide/lipid signalling",
    "MTF1":     "Metal-regulatory transcription factor 1 — oxidative stress",
    "DHRSX":    "Dehydrogenase/reductase X — unknown substrate",
    "PLA2G4E":  "Group IVE cytosolic PLA₂ — arachidonate release, inflammation",
    "SLC1A1":   "EAAT3 — neuronal glutamate transporter, GSH synthesis",
    "ARSD":     "Arylsulfatase D — lysosomal sulfatase",
    "GRIN1":    "NMDA receptor subunit GluN1 — Ca²⁺ signalling",
    "RRM2":     "Ribonucleotide reductase M2 — dNTP synthesis, DNA repair",
    "PFKFB1":   "6-phosphofructo-2-kinase — glycolysis regulation (liver)",
    "SLC25A15": "Mitochondrial ornithine carrier — urea cycle",
    "SLC6A3":   "Dopamine transporter (DAT) — neurotransmitter reuptake",
    "SLC44A1":  "Choline transporter-like protein 1 — membrane synthesis",
    "PTGS1":    "Cyclooxygenase-1 (COX-1) — prostaglandin synthesis",
}

# Known depleted/enriched for annotation
KNOWN_DEPLETED = {"CYP2E1", "ASL", "CYB5A", "MTF1", "NAT2", "GBA3", "TRPM7",
                  "ATP6V0A4", "TPCN1", "HSD11B1", "PFKFB1"}
KNOWN_ENRICHED = {"PLA2G4E", "SLC1A1", "ARSD", "GRIN1", "RRM2", "SLC6A3",
                  "SLC44A1", "DHRSX", "ALOXE3", "PTGS1", "SLC25A15"}

# ═══════════════════════════════════════════════════════════════════════════════
# Thesis Configuration
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ThesisConfig:
    """Centralised, immutable pipeline configuration."""
    # Paths
    project_root: Path = field(default_factory=lambda: Path(__file__).resolve().parent)
    crispr_root: Path = field(default_factory=lambda:
        Path(__file__).resolve().parent.parent / "CRISPR-Seq_Analysis")
    fastq_dir: Path = field(default_factory=lambda:
        Path(__file__).resolve().parent.parent / "CRISPR-Seq_Analysis" / "Fastq")
    library_file: Path = field(default_factory=lambda:
        Path(__file__).resolve().parent.parent / "CRISPR-Seq_Analysis"
        / "nf-core_crisprseq" / "library_sabatini_metko.tsv")
    output_dir: Path = field(default_factory=lambda: Path(__file__).resolve().parent)

    # Analysis parameters
    pval_threshold: float = 0.05
    fdr_threshold: float = 0.25
    top_n_label: int = 10
    scaffold_seq: str = "GTTTTAGAGCTAGAAATAGC"   # 20 bp, lentiCRISPR v1
    sgrna_len: int = 20
    norm_method: str = "median"
    gene_lfc_method: str = "median"
    remove_zero: str = "both"
    remove_zero_threshold: int = 0
    random_seed: int = 42

    # Figure settings
    figure_dpi: int = 300
    figure_format: str = "pdf"   # primary format for thesis
    figure_secondary_format: str = "png"
    figure_width_single: float = 4.72    # inches (A4 single column, 120 mm)
    figure_width_double: float = 6.69    # inches (A4 double column, 170 mm)
    font_family: str = "DejaVu Serif"

    # QC thresholds
    min_reads_per_sample: int = 9_000_000
    ideal_mapping_pct: float = 60.0
    ideal_gini: float = 0.2
    ideal_zero_pct: float = 20.0
    ideal_detected_pct: float = 90.0

    def __post_init__(self):
        """Validate configuration after initialisation."""
        if not 0 < self.pval_threshold <= 1:
            raise ValueError(f"pval_threshold must be in (0, 1], got {self.pval_threshold}")
        if self.random_seed is not None:
            np.random.seed(self.random_seed)

    @property
    def project_id(self) -> str:
        return PROJECT_ID

    @property
    def pipeline_version(self) -> str:
        return PIPELINE_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Export config as JSON-serialisable dict (excludes Paths)."""
        d = asdict(self)
        for k, v in d.items():
            if isinstance(v, Path):
                d[k] = str(v)
        return d


# ═══════════════════════════════════════════════════════════════════════════════
# Thesis Figure Styling
# ═══════════════════════════════════════════════════════════════════════════════

class ThesisStyle:
    """
    Standardised thesis figure styling.

    Applies consistent font (DejaVu Serif), colour palette, DPI,
    and spine configuration to all figures. Use as context manager
    or call ThesisStyle.apply() once at script start.
    """

    # Wong 2011 colourblind-friendly palette (Nature Methods)
    PALETTE_8 = [
        "#000000", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
    ]

    @staticmethod
    def apply(*, font_family: str = "DejaVu Serif", dpi: int = 300,
              font_scale: float = 1.0) -> None:
        """Apply thesis-standard matplotlib rcParams."""
        sns.set_theme(style="ticks", font_scale=font_scale)
        plt.rcParams.update({
            # Resolution
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            # Font
            "font.family": font_family,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            # Layout
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            # Colour
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "grid.color": "#E0E0E0",
            "grid.linewidth": 0.5,
            # Legend
            "legend.frameon": False,
            "legend.handlelength": 1.5,
            "legend.handletextpad": 0.5,
            # SVG/PDF quality
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
        })

    @staticmethod
    def legend_patches() -> list[mpatches.Patch]:
        """Return standard legend patches for depleted/enriched/NS."""
        return [
            mpatches.Patch(color=C_DEPLETED, label="Depleted (essential for TRAIL sensitivity)"),
            mpatches.Patch(color=C_ENRICHED, label="Enriched (confers TRAIL resistance)"),
            mpatches.Patch(color=C_NS, label="Not significant"),
        ]

    @staticmethod
    def condition_legend() -> list[mpatches.Patch]:
        """Return legend patches for experimental conditions."""
        labels = {
            "R0": "Resistant Day-0",
            "Rpost": "Resistant Post-TRAIL",
            "S0": "Sensitive Day-0",
            "Spost": "Sensitive Post-TRAIL",
        }
        return [mpatches.Patch(color=c, label=labels.get(g, g))
                for g, c in COND_COLORS.items()]


# ═══════════════════════════════════════════════════════════════════════════════
# Logging
# ═══════════════════════════════════════════════════════════════════════════════

_LOGGER: logging.Logger | None = None


def setup_logging(
    name: str = "crispr_pipeline",
    level: int = logging.INFO,
    log_file: Path | None = None,
) -> logging.Logger:
    """
    Configure thesis-standard logging.

    Parameters
    ----------
    name : str
        Logger name.
    level : int
        Logging level (default: INFO).
    log_file : Path | None
        If provided, also write logs to this file.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    global _LOGGER

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # File handler (optional)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file, mode="a")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    _LOGGER = logger
    return logger


def get_logger() -> logging.Logger:
    """Return the module-level logger, creating a default one if needed."""
    global _LOGGER
    if _LOGGER is None:
        _LOGGER = setup_logging()
    return _LOGGER


# ═══════════════════════════════════════════════════════════════════════════════
# Figure I/O
# ═══════════════════════════════════════════════════════════════════════════════

def save_figure(
    fig: plt.Figure,
    name: str,
    outdir: Path,
    primary_fmt: str = "pdf",
    secondary_fmt: str = "png",
    dpi: int = 300,
    close: bool = True,
) -> list[Path]:
    """
    Save a figure in primary and secondary formats.

    Parameters
    ----------
    fig : plt.Figure
        Matplotlib figure to save.
    name : str
        Base filename (without extension).
    outdir : Path
        Output directory.
    primary_fmt : str
        Primary format (default: pdf for thesis).
    secondary_fmt : str
        Secondary format (default: png for quick preview).
    dpi : int
        Resolution.
    close : bool
        Close figure after saving.

    Returns
    -------
    list[Path]
        Paths to saved files.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []

    # Primary format
    p1 = outdir / f"{name}.{primary_fmt}"
    # PDF uses vector; DPI ignored
    pdf_dpi = dpi if primary_fmt != "pdf" else None
    fig.savefig(p1, bbox_inches="tight", dpi=pdf_dpi or dpi,
                pad_inches=0.05)
    paths.append(p1)

    # Secondary format
    if secondary_fmt and secondary_fmt != primary_fmt:
        p2 = outdir / f"{name}.{secondary_fmt}"
        fig.savefig(p2, bbox_inches="tight", dpi=dpi, pad_inches=0.05)
        paths.append(p2)

    if close:
        plt.close(fig)

    get_logger().info("✓  %s  (%s)", name, ", ".join(p.suffix for p in paths))
    return paths


# ═══════════════════════════════════════════════════════════════════════════════
# Metadata Export
# ═══════════════════════════════════════════════════════════════════════════════

def export_metadata(
    outdir: Path,
    *,
    step: str = "",
    extra: dict[str, Any] | None = None,
    config: ThesisConfig | None = None,
) -> Path:
    """
    Export pipeline run metadata as JSON for reproducibility tracking.

    Parameters
    ----------
    outdir : Path
        Output directory.
    step : str
        Pipeline step identifier (e.g., "02_mageck_count").
    extra : dict | None
        Additional key-value pairs to include.
    config : ThesisConfig | None
        If provided, include config snapshot.

    Returns
    -------
    Path
        Path to metadata JSON file.
    """
    meta: dict[str, Any] = {
        "project_id": PROJECT_ID,
        "project_title": PROJECT_TITLE,
        "author": AUTHOR,
        "lab": LAB,
        "library": LIBRARY,
        "pipeline_version": PIPELINE_VERSION,
        "step": step,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
    }
    if config is not None:
        meta["config"] = config.to_dict()
    if extra:
        meta.update(extra)

    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"metadata_{step}.json"
    with open(path, "w") as f:
        json.dump(meta, f, indent=2, default=str)
    return path


# ═══════════════════════════════════════════════════════════════════════════════
# Data Loading Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def load_count_table(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """
    Load and parse the MAGeCK count table.

    Returns
    -------
    tuple
        (counts_raw, count_matrix, sample_groups)
        - counts_raw: full count DataFrame with Gene column
        - count_matrix: numeric count DataFrame (sample columns only)
        - sample_groups: {condition: [sample_col_names]}
    """
    if not path.exists():
        raise FileNotFoundError(f"Count table not found: {path}")

    df = pd.read_csv(path, sep="\t", index_col=0)

    # Remove intergenic controls
    if "Gene" in df.columns:
        counts_raw = df[~df["Gene"].str.upper().str.contains("INTERGENIC", na=False)]
        count_matrix = counts_raw.drop(columns=["Gene"])
    else:
        counts_raw = df
        count_matrix = df

    # Build sample groups from column names
    sample_groups: dict[str, list[str]] = {}
    for col in count_matrix.columns:
        if col.startswith("R0"):
            sample_groups.setdefault("R0", []).append(col)
        elif col.startswith("S0"):
            sample_groups.setdefault("S0", []).append(col)
        elif col.startswith("R") and not col.startswith("R0"):
            sample_groups.setdefault("Rpost", []).append(col)
        elif col.startswith("S") and not col.startswith("S0"):
            sample_groups.setdefault("Spost", []).append(col)

    return counts_raw, count_matrix, sample_groups


def load_gene_summary(path: Path, pval_thresh: float = 0.05) -> pd.DataFrame | None:
    """
    Load MAGeCK gene_summary.txt with derived columns.

    Parameters
    ----------
    path : Path
        Path to gene_summary.txt.
    pval_thresh : float
        p-value threshold for significance flag.

    Returns
    -------
    pd.DataFrame | None
        DataFrame with columns: LFC, pval, FDR, log10p, sig, direction.
        Returns None if file not found.
    """
    if not path.exists():
        return None

    df = pd.read_csv(path, sep="\t")
    df["LFC"] = df["neg|lfc"]
    df["pval"] = df[["neg|p-value", "pos|p-value"]].min(axis=1)
    df["FDR"] = df[["neg|fdr", "pos|fdr"]].min(axis=1)
    df["log10p"] = -np.log10(df["pval"].clip(lower=1e-100))
    df["sig"] = df["pval"] < pval_thresh
    df["direction"] = np.where(df["LFC"] < 0, "Depleted", "Enriched")
    return df


def load_all_comparisons(test_dir: Path, pval_thresh: float = 0.05
                         ) -> dict[str, pd.DataFrame]:
    """
    Load all three MAGeCK comparisons from a test directory.

    Parameters
    ----------
    test_dir : Path
        Directory containing {comp}.gene_summary.txt files.
    pval_thresh : float
        p-value threshold.

    Returns
    -------
    dict[str, pd.DataFrame]
        {comparison_name: gene_summary_DataFrame}
    """
    data: dict[str, pd.DataFrame] = {}
    for cname in COMPARISONS:
        df = load_gene_summary(test_dir / f"{cname}.gene_summary.txt", pval_thresh)
        if df is not None:
            data[cname] = df
            logger = get_logger()
            logger.info("  %-22s: %4d genes, %2d significant",
                        cname, len(df), df["sig"].sum())
    return data


# ═══════════════════════════════════════════════════════════════════════════════
# Statistical Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def compute_gini(values: np.ndarray) -> float:
    """
    Compute the Gini inequality index.

    Gini = (2·Σ(i·xᵢ) / (n·Σxᵢ)) − (n+1)/n

    Parameters
    ----------
    values : np.ndarray
        Array of non-negative values.

    Returns
    -------
    float
        Gini index [0, 1]. 0 = perfect equality, 1 = maximum inequality.
    """
    v = values[values > 0]
    if len(v) == 0:
        return 0.0
    sv = np.sort(v)
    n = len(sv)
    total = np.sum(sv)
    if total == 0:
        return 0.0
    g = (2.0 * np.dot(np.arange(1, n + 1, dtype=np.float64), sv)
         / (n * total) - (n + 1) / n)
    return float(np.clip(g, 0.0, 1.0))


def classify_condition(col_name: str) -> str:
    """Classify a sample column name into experimental condition."""
    if col_name.startswith("R0"):
        return "R0"
    if col_name.startswith("S0"):
        return "S0"
    if col_name.startswith("R") and not col_name.startswith("R0"):
        return "Rpost"
    if col_name.startswith("S") and not col_name.startswith("S0"):
        return "Spost"
    return "Unknown"


def bonferroni_correct(pvals: np.ndarray, n_tests: int | None = None) -> np.ndarray:
    """
    Apply Bonferroni correction for multiple testing.

    Parameters
    ----------
    pvals : np.ndarray
        Array of p-values.
    n_tests : int | None
        Number of tests. If None, uses len(pvals).

    Returns
    -------
    np.ndarray
        Bonferroni-corrected p-values, clipped to [0, 1].
    """
    if n_tests is None:
        n_tests = len(pvals)
    return np.clip(pvals * n_tests, 0.0, 1.0)


# ═══════════════════════════════════════════════════════════════════════════════
# Shell Command Runner
# ═══════════════════════════════════════════════════════════════════════════════

def run_shell(cmd: str, *, timeout: int | None = 600) -> int:
    """
    Execute a shell command with logging.

    Parameters
    ----------
    cmd : str
        Shell command to execute.
    timeout : int | None
        Timeout in seconds (default: 600 for MAGeCK operations).

    Returns
    -------
    int
        Return code (0 = success).
    """
    logger = get_logger()
    display = cmd[:160] + "…" if len(cmd) > 160 else cmd
    logger.info("  $ %s", display)

    try:
        result = subprocess.run(
            cmd, shell=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        logger.error("  Command timed out after %ds: %s", timeout, display)
        return 1

    # Log key lines
    for line in result.stdout.splitlines():
        lowered = line.lower()
        if any(k in lowered for k in ("error", "warn", "mapped", "total",
                                        "testable", "gene", "sgrna", "read")):
            logger.info("    %s", line.strip())

    if result.returncode != 0:
        logger.error("  Exit code %d", result.returncode)
        logger.error(result.stdout[-2000:])

    return result.returncode


# ═══════════════════════════════════════════════════════════════════════════════
# AdjustText availability
# ═══════════════════════════════════════════════════════════════════════════════

HAS_ADJUSTTEXT: bool
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    warnings.warn("adjustText not installed. Gene labels may overlap. "
                  "Install: pip install adjustText")


def label_top_genes(
    ax: plt.Axes,
    subset: pd.DataFrame,
    x_col: str,
    y_col: str,
    n: int = 10,
) -> None:
    """
    Add gene name labels to the top N most significant points.

    Uses adjustText for automatic label placement when available.
    """
    if subset.empty:
        return
    top = subset.nsmallest(n, "pval") if "pval" in subset.columns else subset.head(n)
    texts = []
    for _, row in top.iterrows():
        color = C_DEPLETED if row.get("LFC", 0) < 0 else C_ENRICHED
        t = ax.text(
            row[x_col], row[y_col], row["id"],
            fontsize=6.5, color=color, fontweight="bold",
            ha="center", va="bottom", zorder=5,
        )
        texts.append(t)
    if HAS_ADJUSTTEXT and len(texts) > 1:
        adjust_text(
            texts, ax=ax,
            expand_text=(1.5, 1.5), expand_points=(1.8, 1.8),
            force_text=(0.4, 0.6), force_points=(0.2, 0.4),
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        )


# ── Module initialisation ────────────────────────────────────────────────────
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="seaborn")
