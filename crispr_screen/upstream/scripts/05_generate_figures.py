#!/usr/bin/env python3
"""
05_generate_figures.py — Publication-Quality CRISPR Screen Figures
==================================================================

Generates manuscript- and thesis-quality figures for the chapter
"CRISPR/Cas9 Metabolic Screen of TRAIL Resistance in GBM."

Design language (curated from established references):
  - SciencePlots / matplotlib_for_papers ........ clean rcParams, thin spines,
                                                   consistent typography, vector export
  - scientific-visualization-book ................ restraint: no chartjunk, no 3D,
                                                   no dark backgrounds, generous whitespace
  - ICWallis publication-ready-figures ........... tidy Matplotlib/Seaborn idioms
  - kevinblighe/EnhancedVolcano .................. volcano grammar: depleted / enriched /
                                                   non-significant groups + thresholds
  - jokergoo/ComplexHeatmap ...................... diverging palette centred at zero

Engineering principles (this refactor):
  - Analysis logic is UNCHANGED — every value comes from _thesis_utils loaders.
    This module only governs plotting, layout, typography, colour and export.
  - ONE reusable global style function: set_publication_style().
  - Colourblind-safe, print-friendly Wong (2011) palette, shared via _thesis_utils.
  - Consistent font sizes everywhere via the FONT registry.
  - Every figure exported as vector PDF (editable text) + raster PNG (300/600 dpi).
  - Sparse, prioritised gene labelling (adjustText when available) — never crowded.
  - Panel labels (A, B, C, …) on every multi-panel figure.

Figures (one file each — vector PDF + PNG):
  Fig1 — Volcano plots (3-panel)
  Fig2 — Gene rank plots (3-panel)
  Fig3 — LFC heatmap (multi-comparison hits)
  Fig4 — Top-gene bar charts
  Fig5 — QC metrics (4-panel)
  Fig6 — sgRNA count profiles
  Fig7 — Replicate correlations
  Fig8 — Cross-comparison dot plot
  Fig9 — Significant-gene summary table

Project: PRJ25-131 | Cingöz Lab | 2026
"""

from __future__ import annotations

import argparse
import sys
import textwrap
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _thesis_utils import (  # analysis + colour contract — DO NOT reimplement
    COMPARISONS,
    COND_COLORS,
    C_DEPLETED,
    C_ENRICHED,
    C_NS,
    FUNC_ANNOT,
    compute_gini,
    export_metadata,
    load_all_comparisons,
    load_count_table,
    setup_logging,
)

logger = setup_logging("05_figures")
warnings.filterwarnings("ignore")

# adjustText is an optional polish dependency (declared in environment.yml).
# If unavailable we fall back to deterministic offset labelling — never crash.
try:
    from adjustText import adjust_text

    _HAS_ADJUST = True
except ImportError:  # pragma: no cover - environment dependent
    _HAS_ADJUST = False
    logger.warning("adjustText not found — using offset label fallback.")


# ═══════════════════════════════════════════════════════════════════════════════
# Global design tokens — colours, typography, sizing (single source of truth)
# ═══════════════════════════════════════════════════════════════════════════════

# Colourblind-safe Wong (2011, Nature Methods) palette, imported from _thesis_utils
# so analysis scripts and figures never drift apart.
COL_DEPLETED = C_DEPLETED   # blue       — depleted / dropout
COL_ENRICHED = C_ENRICHED   # vermillion — enriched
COL_NS = C_NS               # neutral grey — not significant
COL_REF = "#555555"         # threshold / reference lines
COL_AXIS = "#FFFFFF"        # marker edge for contrast on print
COL_GRID = "#E9E9E9"        # faint grid / gridlines
COL_NA = "#EDEDED"          # "not testable" heatmap cells

# Diverging colormap for heatmaps — RdBu reversed = blue(low) → white(0) → red(high),
# perceptually balanced and colourblind-distinguishable, centred at zero (ComplexHeatmap idiom).
DIVERGING_CMAP = "RdBu_r"

# One consistent type scale across ALL figures (points).
FONT = {
    "base": 8,
    "tick": 7.5,
    "label": 8.5,
    "title": 9.5,
    "suptitle": 10.5,
    "panel": 11,      # bold A/B/C panel letters
    "annot": 6.5,     # in-plot annotations (gene labels, cell values)
    "legend": 8,
}

# Standard manuscript widths (inches): single ≈ 3.4, 1.5-col ≈ 5.0, double ≈ 7.0.
W_SINGLE, W_ONEHALF, W_DOUBLE = 3.4, 5.0, 7.0


def set_publication_style() -> None:
    """Apply the single global publication style (call once before plotting).

    Combines SciencePlots-style minimalism with matplotlib_for_papers typography
    rules and vector-safe export defaults. Purely cosmetic — affects no analysis.
    """
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            # --- resolution & export -------------------------------------------------
            "figure.dpi": 200,
            "savefig.dpi": 600,                # crisp rasterised scatter inside the PDF
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "savefig.transparent": False,
            "pdf.fonttype": 42,                # editable (TrueType) text in vector PDF
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            # --- typography (Arial-metric fallback chain; mathtext for sub/superscripts)
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial", "Helvetica", "Nimbus Sans",
                "Liberation Sans", "DejaVu Sans",
            ],
            "mathtext.default": "regular",     # math renders in the body sans font
            "font.size": FONT["base"],
            "axes.titlesize": FONT["title"],
            "axes.labelsize": FONT["label"],
            "xtick.labelsize": FONT["tick"],
            "ytick.labelsize": FONT["tick"],
            "legend.fontsize": FONT["legend"],
            "figure.titlesize": FONT["suptitle"],
            # --- axes & spines (clean, sparse) --------------------------------------
            "axes.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titlepad": 6,
            "axes.titleweight": "bold",
            "axes.edgecolor": "#333333",
            "axes.labelcolor": "#1A1A1A",
            "text.color": "#1A1A1A",
            "xtick.color": "#333333",
            "ytick.color": "#333333",
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            # --- backgrounds & grid (white, no chartjunk) ---------------------------
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.grid": False,
            "grid.color": COL_GRID,
            "grid.linewidth": 0.5,
            # --- legend --------------------------------------------------------------
            "legend.frameon": False,
            "legend.handletextpad": 0.5,
            "legend.columnspacing": 1.2,
            "legend.borderpad": 0.3,
        }
    )


def save_figure(fig: plt.Figure, name: str, outdir: Path, png_dpi: int = 600) -> None:
    """Export a figure as vector PDF (publication) and raster PNG (preview)."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{name}.pdf")                       # vector, editable text
    fig.savefig(outdir / f"{name}.png", dpi=png_dpi)          # 300/600 dpi preview
    plt.close(fig)
    logger.info("  ✓  %s  (PDF + PNG @ %d dpi)", name, png_dpi)


def add_panel_label(ax: plt.Axes, letter: str, *, x: float = -0.06, y: float = 1.06) -> None:
    """Place a bold panel label (A, B, C, …) in axes-fraction coordinates."""
    ax.text(
        x, y, letter, transform=ax.transAxes,
        fontsize=FONT["panel"], fontweight="bold",
        va="bottom", ha="right", clip_on=False,
    )


@dataclass
class FigContext:
    """Shared, read-only context handed to every figure function."""

    data: dict[str, pd.DataFrame]
    counts: pd.DataFrame | None
    count_matrix: pd.DataFrame | None
    sample_groups: dict[str, list[str]]
    sig_genes: list[str]
    outdir: Path
    pval_thresh: float = 0.05
    png_dpi: int = 600


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers — condition assignment & sparse gene labelling
# ═══════════════════════════════════════════════════════════════════════════════

def _condition_of(sample: str) -> str:
    """Map a sample column name to its experimental condition (display only)."""
    if sample.startswith("R0"):
        return "R0"
    if sample.startswith("S0"):
        return "S0"
    if sample.startswith("R"):
        return "Rpost"
    return "Spost"


def _label_top_genes(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    *,
    max_labels: int = 8,
) -> None:
    """Label only the most significant genes, avoiding overlap.

    Uses adjustText (with thin leader lines) when available; otherwise falls back
    to deterministic alternating offsets. Labels are coloured by direction so the
    annotation reinforces the depleted/enriched grouping.
    """
    sub = df.nsmallest(max_labels, "pval")
    if sub.empty:
        return

    texts = []
    for _, row in sub.iterrows():
        color = COL_DEPLETED if row["LFC"] < 0 else COL_ENRICHED
        texts.append(
            ax.text(
                row[x_col], row[y_col], row["id"],
                fontsize=FONT["annot"], fontweight="bold", color=color,
                ha="center", va="center",
            )
        )

    if _HAS_ADJUST:
        # Strong text–text repulsion with thin leader lines keeps tightly
        # clustered hits (e.g. co-enriched genes at one corner) from overprinting.
        adjust_text(
            texts, ax=ax,
            expand=(1.5, 2.0),
            force_text=(0.5, 1.0),
            force_static=(0.2, 0.5),
            force_pull=(0.005, 0.005),
            min_arrow_len=4,
            arrowprops=dict(arrowstyle="-", color="#999999", lw=0.4, shrinkA=3),
            only_move={"text": "xy", "static": "xy", "explode": "xy"},
        )
        return

    # Fallback: deterministic alternating offsets (no overlap for small N).
    offsets = [(8, 4), (-8, -4), (8, -4), (-8, 4), (12, 2), (-12, -2), (6, 6), (-6, -6)]
    for i, t in enumerate(texts):
        ox, oy = offsets[i % len(offsets)]
        t.set_position(
            ax.transData.inverted().transform(
                ax.transData.transform(t.get_position()) + np.array([ox, oy])
            )
        )
        t.set_bbox(dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))


def _volcano_legend(fig: plt.Figure) -> None:
    """Shared EnhancedVolcano-style legend placed below the panel row."""
    handles = [
        Line2D([], [], marker="o", ls="", ms=6, mec=COL_AXIS, mew=0.4,
               color=COL_DEPLETED, label="Depleted (sensitiser)"),
        Line2D([], [], marker="o", ls="", ms=6, mec=COL_AXIS, mew=0.4,
               color=COL_ENRICHED, label="Enriched (resistance)"),
        Line2D([], [], marker="o", ls="", ms=5, color=COL_NS, label="Not significant"),
    ]
    fig.legend(
        handles=handles, loc="lower center", ncol=3,
        bbox_to_anchor=(0.5, -0.04), frameon=False, fontsize=FONT["legend"],
    )


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Volcano plots  (EnhancedVolcano grammar)
# ═══════════════════════════════════════════════════════════════════════════════

def fig1_volcano(ctx: FigContext) -> None:
    """Three-panel volcano: depleted / enriched / non-significant + thresholds."""
    logger.info("Fig 1 — Volcano plots")

    fig, axes = plt.subplots(1, 3, figsize=(W_DOUBLE, 2.9))
    fig.subplots_adjust(wspace=0.42, bottom=0.22, top=0.86)

    for ax, letter, (cname, clabel) in zip(axes, "ABC", COMPARISONS.items()):
        df = ctx.data.get(cname)
        if df is None:
            ax.set_visible(False)
            continue

        thr = -np.log10(ctx.pval_thresh)

        ns = df[~df["sig"]]
        dep = df[df["sig"] & (df["LFC"] < 0)]
        enr = df[df["sig"] & (df["LFC"] > 0)]

        # Layer order: NS (rasterised, faint) → significant groups on top.
        ax.scatter(ns["LFC"], ns["log10p"], s=10, c=COL_NS, linewidths=0,
                   alpha=0.5, zorder=1, rasterized=True)
        ax.scatter(dep["LFC"], dep["log10p"], s=30, c=COL_DEPLETED,
                   edgecolors=COL_AXIS, linewidths=0.3, alpha=0.95, zorder=3)
        ax.scatter(enr["LFC"], enr["log10p"], s=30, c=COL_ENRICHED,
                   edgecolors=COL_AXIS, linewidths=0.3, alpha=0.95, zorder=3)

        # Threshold guides (subtle, dashed).
        ax.axhline(thr, color=COL_REF, lw=0.7, ls=(0, (4, 3)), zorder=2)
        ax.axvline(0, color="#CCCCCC", lw=0.5, zorder=1)

        # Sparse labelling only when the panel is not saturated.
        if len(df[df["sig"]]) <= 25:
            _label_top_genes(ax, df[df["sig"]], "LFC", "log10p", max_labels=8)

        # Compact count summary in the lower corners — the reliably empty region
        # of a volcano (significant-gene labels cluster at the top), so it never
        # collides with the adjustText labels above.
        ax.text(0.03, 0.04, f"↓ {len(dep)}", transform=ax.transAxes,
                fontsize=FONT["annot"], color=COL_DEPLETED, ha="left", va="bottom",
                fontweight="bold")
        ax.text(0.97, 0.04, f"↑ {len(enr)}", transform=ax.transAxes,
                fontsize=FONT["annot"], color=COL_ENRICHED, ha="right", va="bottom",
                fontweight="bold")

        ax.set_title(textwrap.fill(clabel, 26), fontsize=FONT["title"])
        ax.set_xlabel(r"Log$_2$ fold change")
        ax.set_ylabel(r"$-\log_{10}\,p$")
        add_panel_label(ax, letter)
        sns.despine(ax=ax)

    _volcano_legend(fig)
    save_figure(fig, "Fig1_volcano", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Gene rank plots
# ═══════════════════════════════════════════════════════════════════════════════

def fig2_rank(ctx: FigContext) -> None:
    """Ranked LFC ('hockey-stick') plots with significant genes highlighted."""
    logger.info("Fig 2 — Rank plots")

    fig, axes = plt.subplots(1, 3, figsize=(W_DOUBLE, 2.9))
    fig.subplots_adjust(wspace=0.42, bottom=0.22, top=0.86)

    for ax, letter, (cname, clabel) in zip(axes, "ABC", COMPARISONS.items()):
        df = ctx.data.get(cname)
        if df is None:
            ax.set_visible(False)
            continue

        ranked = df.sort_values("LFC").reset_index(drop=True)
        ranked["rank"] = np.arange(len(ranked))

        ns = ranked[~ranked["sig"]]
        dep = ranked[ranked["sig"] & (ranked["LFC"] < 0)]
        enr = ranked[ranked["sig"] & (ranked["LFC"] > 0)]

        ax.scatter(ns["rank"], ns["LFC"], s=8, c=COL_NS, linewidths=0,
                   alpha=0.45, zorder=1, rasterized=True)
        ax.scatter(dep["rank"], dep["LFC"], s=28, c=COL_DEPLETED,
                   edgecolors=COL_AXIS, linewidths=0.3, alpha=0.95, zorder=3)
        ax.scatter(enr["rank"], enr["LFC"], s=28, c=COL_ENRICHED,
                   edgecolors=COL_AXIS, linewidths=0.3, alpha=0.95, zorder=3)

        ax.axhline(0, color="#CCCCCC", lw=0.5, zorder=1)

        if len(ranked[ranked["sig"]]) <= 25:
            _label_top_genes(ax, ranked[ranked["sig"]], "rank", "LFC", max_labels=8)

        ax.set_title(textwrap.fill(clabel, 26), fontsize=FONT["title"])
        ax.set_xlabel("Gene rank")
        ax.set_ylabel(r"Log$_2$ fold change")
        add_panel_label(ax, letter)
        sns.despine(ax=ax)

    _volcano_legend(fig)
    save_figure(fig, "Fig2_rank_lfc", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — LFC heatmap  (ComplexHeatmap idiom: diverging, centred at 0)
# ═══════════════════════════════════════════════════════════════════════════════

def fig3_heatmap(ctx: FigContext) -> None:
    """Clean diverging LFC heatmap across comparisons; grey = not testable."""
    logger.info("Fig 3 — Heatmap")

    if not ctx.sig_genes:
        logger.warning("  [SKIP] No significant genes")
        return

    heat = pd.DataFrame(
        {
            clabel: df.set_index("id")["LFC"].reindex(ctx.sig_genes)
            for cname, clabel in COMPARISONS.items()
            if (df := ctx.data.get(cname)) is not None
        }
    )
    # Order genes by overall effect magnitude (most striking on top).
    heat = heat.loc[heat.abs().sum(axis=1).sort_values(ascending=False).index]

    vmax = float(min(np.nanmax(np.abs(heat.values)) + 0.5, 16.0))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    n_rows, n_cols = heat.shape
    row_h = max(0.26, 4.0 / max(n_rows, 1))
    fig_h = max(3.2, n_rows * row_h + 1.4)
    fig, ax = plt.subplots(figsize=(W_ONEHALF, fig_h))

    display = np.ma.masked_invalid(heat.values)  # NaN drawn via 'bad' colour
    cmap = plt.get_cmap(DIVERGING_CMAP).copy()
    cmap.set_bad(COL_NA)
    im = ax.imshow(display, aspect="auto", cmap=cmap, norm=norm)

    # In-cell LFC values, contrast-aware (skip NaN / zero).
    for r in range(n_rows):
        for c in range(n_cols):
            v = heat.values[r, c]
            if not np.isnan(v) and v != 0:
                tc = "white" if abs(v) > vmax * 0.55 else "#1A1A1A"
                ax.text(c, r, f"{v:.1f}", ha="center", va="center",
                        fontsize=FONT["annot"], color=tc, fontweight="bold")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([textwrap.fill(c, 16) for c in heat.columns],
                       fontsize=FONT["tick"], rotation=20, ha="right")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(heat.index, fontsize=FONT["tick"])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.03, shrink=0.7)
    cbar.set_label(r"Log$_2$ fold change", fontsize=FONT["label"])
    cbar.ax.tick_params(labelsize=FONT["tick"])
    cbar.outline.set_visible(False)

    ax.set_title(
        f"Significant genes — LFC across comparisons\n"
        f"(grey = not testable; $p$ < {ctx.pval_thresh})",
        fontsize=FONT["title"], pad=8,
    )
    fig.tight_layout()
    save_figure(fig, "Fig3_heatmap", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Top-gene bar charts
# ═══════════════════════════════════════════════════════════════════════════════

def fig4_barplot(ctx: FigContext) -> None:
    """Horizontal bars of the top genes per comparison, coloured by direction."""
    logger.info("Fig 4 — Barplots")

    comps = [(c, l) for c, l in COMPARISONS.items() if ctx.data.get(c) is not None]
    if not comps:
        logger.warning("  [SKIP] No comparisons")
        return

    n = len(comps)
    fig, axes = plt.subplots(n, 1, figsize=(W_ONEHALF + 1.0, 2.7 * n))
    if n == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.5)

    for ax, letter, (cname, clabel) in zip(axes, "ABCDEF", comps):
        df = ctx.data[cname]
        top = df.nsmallest(15, "pval").sort_values("log10p")
        colors = [COL_DEPLETED if v < 0 else COL_ENRICHED for v in top["LFC"]]

        ax.barh(range(len(top)), top["log10p"].values, color=colors,
                edgecolor="none", height=0.62)

        xmax = top["log10p"].max() * 1.42
        for i, (_, row) in enumerate(top.iterrows()):
            ax.text(row["log10p"] + xmax * 0.012, i, f"$p$={row['pval']:.3f}",
                    va="center", fontsize=FONT["annot"], color="#555555")

        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["id"].values, fontsize=FONT["tick"])
        ax.axvline(-np.log10(ctx.pval_thresh), color=COL_REF, lw=0.7, ls=(0, (3, 3)))
        ax.set_xlabel(r"$-\log_{10}\,p$")
        ax.set_title(clabel, fontsize=FONT["title"])
        ax.set_xlim(0, xmax)
        ax.grid(axis="x", lw=0.4, color=COL_GRID)
        ax.set_axisbelow(True)
        add_panel_label(ax, letter, x=-0.32)
        sns.despine(ax=ax, left=True)

    save_figure(fig, "Fig4_barplot", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — QC metrics  (thresholds clear, not noisy)
# ═══════════════════════════════════════════════════════════════════════════════

def fig5_qc(ctx: FigContext) -> None:
    """Four-panel library/sequencing QC with subtle reference thresholds."""
    logger.info("Fig 5 — QC metrics")

    if ctx.count_matrix is None:
        logger.warning("  [SKIP] No count data")
        return

    cm = ctx.count_matrix
    total_sgrnas = len(cm)

    rows = []
    for col in cm.columns:
        vals = cm[col].values.astype(float)
        rows.append(
            {
                "sample": col,
                "group": _condition_of(col),
                "total_M": vals.sum() / 1e6,
                "nonzero": int((vals > 0).sum()),
                "zero_pct": float((vals == 0).mean() * 100),
                "gini": compute_gini(vals),
            }
        )
    qc = pd.DataFrame(rows)
    colors = [COND_COLORS.get(g, "#888888") for g in qc["group"]]

    fig, axes = plt.subplots(2, 2, figsize=(W_DOUBLE, 5.4))
    axes = axes.flatten()
    fig.subplots_adjust(hspace=0.55, wspace=0.32, bottom=0.14, top=0.9)

    panels = [
        ("Total reads per sample", "Reads (millions)", qc["total_M"], None, None),
        ("Zero-count sgRNAs", "Zero-count sgRNAs (%)", qc["zero_pct"], 50, "50%"),
        ("Gini inequality index", "Gini index", qc["gini"], 0.2, "0.2"),
        ("Detected sgRNAs", "Non-zero sgRNAs", qc["nonzero"],
         0.90 * total_sgrnas, "90% of library"),
    ]

    for ax, letter, (title, ylabel, vals, ref, ref_label) in zip(axes, "ABCD", panels):
        ax.bar(range(len(vals)), vals, color=colors, edgecolor="none", width=0.72)
        if ref is not None:
            ax.axhline(ref, color=COL_REF, lw=0.8, ls=(0, (4, 2)))
            ax.text(0.99, ref, f" {ref_label}", transform=ax.get_yaxis_transform(),
                    ha="right", va="bottom", fontsize=FONT["annot"], color=COL_REF)
        ax.set_xticks(range(len(qc)))
        ax.set_xticklabels(qc["sample"], rotation=40, ha="right", fontsize=FONT["annot"])
        ax.set_ylabel(ylabel, fontsize=FONT["label"])
        ax.set_title(title, fontsize=FONT["title"])
        ax.grid(axis="y", lw=0.4, color=COL_GRID)
        ax.set_axisbelow(True)
        add_panel_label(ax, letter)
        sns.despine(ax=ax)

    handles = [mpatches.Patch(color=c, label=g) for g, c in COND_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.01), frameon=False, fontsize=FONT["legend"])

    save_figure(fig, "Fig5_qc_metrics", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — sgRNA count profiles
# ═══════════════════════════════════════════════════════════════════════════════

def fig6_sgrna(ctx: FigContext) -> None:
    """RPM trajectories for key hit genes; only sgRNAs with signal are drawn."""
    logger.info("Fig 6 — sgRNA profiles")

    if ctx.counts is None or ctx.count_matrix is None:
        logger.warning("  [SKIP] No count data")
        return

    key_genes = ["CYP2E1", "ASL", "PLA2G4E", "CYB5A", "MTF1"]
    available = [g for g in key_genes if g in ctx.counts["Gene"].values]
    if not available:
        logger.warning("  [SKIP] Key genes not found")
        return

    sample_totals = ctx.count_matrix.sum(axis=0)
    rpm = ctx.count_matrix.div(sample_totals, axis=1) * 1e6
    sample_order = list(ctx.count_matrix.columns)

    fig, axes = plt.subplots(len(available), 1,
                             figsize=(W_ONEHALF + 1.0, 1.9 * len(available)))
    if len(available) == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.6)

    sg_palette = sns.color_palette("colorblind", 8)  # colourblind-safe per-sgRNA

    for ax, letter, gene in zip(axes, "ABCDEFGH", available):
        rpm_gene = rpm.loc[ctx.counts[ctx.counts["Gene"] == gene].index]
        active = rpm_gene[rpm_gene.sum(axis=1) > 0]
        n_active, n_total = len(active), len(rpm_gene)

        for si, (_, sg_row) in enumerate(active.iterrows()):
            ax.plot(range(len(sample_order)), sg_row[sample_order].values.astype(float),
                    marker="o", ms=3, lw=1.1, color=sg_palette[si % 8],
                    label=f"sg{si + 1}", zorder=3)

        # Light condition shading bands (very faint — orientation, not decoration).
        prev, start = None, 0
        for xi in range(len(sample_order) + 1):
            curr = _condition_of(sample_order[xi]) if xi < len(sample_order) else None
            if curr != prev:
                if prev is not None:
                    ax.axvspan(start - 0.5, xi - 0.5, alpha=0.06,
                               color=COND_COLORS.get(prev, "#EEEEEE"), zorder=0)
                prev, start = curr, xi

        ax.set_xticks(range(len(sample_order)))
        ax.set_xticklabels(sample_order, rotation=40, ha="right", fontsize=FONT["annot"])
        ax.set_ylabel("RPM", fontsize=FONT["label"])
        flag = " ⚠" if n_active == 1 else ""
        ax.set_title(f"{gene}  ({n_active}/{n_total} sgRNAs with signal){flag}",
                     fontsize=FONT["title"])
        ax.set_xlim(-0.5, len(sample_order) - 0.5)
        if n_active:
            ax.legend(fontsize=FONT["annot"], loc="upper right", ncol=2,
                      frameon=True, framealpha=0.9, facecolor="white",
                      edgecolor="none")
        ax.grid(axis="y", lw=0.4, color=COL_GRID)
        ax.set_axisbelow(True)
        add_panel_label(ax, letter, x=-0.08)
        sns.despine(ax=ax)

    save_figure(fig, "Fig6_sgrna_profiles", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 7 — Replicate correlations
# ═══════════════════════════════════════════════════════════════════════════════

def fig7_replicates(ctx: FigContext) -> None:
    """Pairwise replicate scatter grid with Spearman ρ per panel."""
    logger.info("Fig 7 — Replicate correlations")

    if ctx.count_matrix is None:
        logger.warning("  [SKIP] No count data")
        return

    cm = ctx.count_matrix
    # log10(1 + RPM) — matches the axis/title label; +1 keeps zeros finite.
    log_rpm = np.log10(cm.div(cm.sum(axis=0), axis=1) * 1e6 + 1.0)

    pairs = [
        (grp, a, b)
        for grp, samples in ctx.sample_groups.items()
        for i, a in enumerate([s for s in samples if s in log_rpm.columns])
        for b in [s for s in samples if s in log_rpm.columns][i + 1:]
    ]
    if not pairs:
        logger.warning("  [SKIP] No replicate pairs")
        return

    ncols = min(3, len(pairs))
    nrows = int(np.ceil(len(pairs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.1 * ncols, 2.0 * nrows))
    axes_flat = np.atleast_1d(axes).flatten()

    for ax, (grp, s1, s2) in zip(axes_flat, pairs):
        x, y = log_rpm[s1].values, log_rpm[s2].values
        mask = (x > 0) | (y > 0)
        rho, _ = spearmanr(x[mask], y[mask])
        color = COND_COLORS.get(grp, "#888888")

        ax.scatter(x[mask], y[mask], s=4, color=color, alpha=0.3,
                   linewidths=0, rasterized=True)
        lim = max(x[mask].max(), y[mask].max()) * 1.05
        ax.plot([0, lim], [0, lim], color="#999999", ls="--", lw=0.6, alpha=0.6)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel(s1, fontsize=FONT["annot"])
        ax.set_ylabel(s2, fontsize=FONT["annot"])
        ax.set_title(f"{grp}   ρ = {rho:.3f}", fontsize=FONT["tick"], color=color)
        ax.set_aspect("equal", adjustable="box")
        sns.despine(ax=ax)

    for ax in axes_flat[len(pairs):]:
        ax.set_visible(False)

    fig.suptitle(r"Replicate correlation — $\log_{10}$(RPM), Spearman ρ",
                 fontsize=FONT["suptitle"], y=1.02)
    fig.tight_layout()
    save_figure(fig, "Fig7_replicate_corr", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 8 — Cross-comparison dot plot
# ═══════════════════════════════════════════════════════════════════════════════

def fig8_cross(ctx: FigContext) -> None:
    """Dot plot of gene LFC across comparisons (size ∝ significance)."""
    logger.info("Fig 8 — Cross-comparison")

    all_genes = sorted(set().union(*[set(df["id"]) for df in ctx.data.values()]))
    lfc_mat = pd.DataFrame(
        {c: ctx.data[c].set_index("id")["LFC"].reindex(all_genes) for c in ctx.data}
    )
    pval_mat = pd.DataFrame(
        {c: ctx.data[c].set_index("id")["pval"].reindex(all_genes) for c in ctx.data}
    )

    n_sig = (pval_mat < ctx.pval_thresh).sum(axis=1)
    hit_genes = n_sig[n_sig >= 1].sort_values(ascending=False).index[:30].tolist()
    if not hit_genes:
        logger.warning("  [SKIP] No genes")
        return

    lfc_sub, pval_sub = lfc_mat.loc[hit_genes], pval_mat.loc[hit_genes]
    comp_names = list(ctx.data.keys())
    comp_labels = [COMPARISONS[c] for c in comp_names]

    fig_h = max(4.5, len(hit_genes) * 0.26 + 1.5)
    fig, ax = plt.subplots(figsize=(W_DOUBLE, fig_h))

    for xi, cname in enumerate(comp_names):
        for yi, gene in enumerate(hit_genes):
            lfc, pv = lfc_sub.loc[gene, cname], pval_sub.loc[gene, cname]
            if np.isnan(lfc) or np.isnan(pv):
                continue
            sig = pv < ctx.pval_thresh
            ax.scatter(
                xi, yi,
                s=max(25, -np.log10(pv + 1e-100) * 50),
                c=COL_DEPLETED if lfc < 0 else COL_ENRICHED,
                alpha=0.9 if sig else 0.12,
                edgecolors=COL_AXIS, linewidths=0.2,
                zorder=4 if sig else 2,
            )
            if sig:
                tc = "white" if abs(lfc) > 8 else "#1A1A1A"
                ax.text(xi, yi, f"{lfc:+.1f}", ha="center", va="center",
                        fontsize=FONT["annot"] - 1, color=tc, fontweight="bold",
                        zorder=5)

    ax.set_xticks(range(len(comp_names)))
    ax.set_xticklabels([textwrap.fill(l, 22) for l in comp_labels],
                       rotation=15, ha="right", fontsize=FONT["tick"])
    ax.set_yticks(range(len(hit_genes)))
    ax.set_yticklabels(hit_genes, fontsize=FONT["tick"])
    ax.set_xlim(-0.6, len(comp_names) - 0.4)
    ax.set_ylim(-0.5, len(hit_genes) - 0.5)
    ax.grid(True, lw=0.3, color=COL_GRID, zorder=0)
    ax.set_axisbelow(True)

    # Size-legend keyed to significance — placed below the plot so it never
    # collides with the (centred, two-line) title at the top.
    size_handles = [
        Line2D([], [], marker="o", ls="", color="#888888", mec=COL_AXIS, mew=0.2,
               ms=np.sqrt(max(25, -np.log10(p) * 50)),
               label=f"$p$ = {p:g}")
        for p in (0.05, 0.01, 0.001)
    ]
    ax.legend(handles=size_handles, loc="upper center", bbox_to_anchor=(0.5, -0.06),
              ncol=3, fontsize=FONT["annot"], handletextpad=0.4, columnspacing=1.4,
              frameon=False, title="Bubble size  ($p$-value)",
              title_fontsize=FONT["annot"])

    ax.set_title(
        "Significant genes across comparisons\n"
        "(blue = depleted, vermillion = enriched; faded = n.s.)",
        fontsize=FONT["title"], pad=10,
    )
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, "Fig8_cross_comparison", ctx.outdir, ctx.png_dpi)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 9 — Significant-gene summary table
# ═══════════════════════════════════════════════════════════════════════════════

def fig9_table(ctx: FigContext) -> None:
    """Typeset summary table of significant genes; also written to CSV."""
    logger.info("Fig 9 — Summary table")

    sig_rows = [
        {
            "Comparison": clabel,
            "Gene": row["id"],
            "LFC": row["LFC"],
            "p-value": row["pval"],
            "Direction": "Depleted" if row["LFC"] < 0 else "Enriched",
            "Function": FUNC_ANNOT.get(row["id"], "—"),
        }
        for cname, clabel in COMPARISONS.items()
        if (df := ctx.data.get(cname)) is not None
        for _, row in df[df["sig"]].iterrows()
    ]
    if not sig_rows:
        logger.warning("  [SKIP] No significant genes")
        return

    tbl = pd.DataFrame(sig_rows).sort_values(["Comparison", "p-value"]).reset_index(drop=True)
    display = tbl.copy()
    display["LFC"] = display["LFC"].map(lambda v: f"{v:+.2f}")
    display["p-value"] = display["p-value"].map(lambda v: f"{v:.4f}")
    display["Direction"] = display["Direction"].map({"Depleted": "↓ Depleted",
                                                     "Enriched": "↑ Enriched"})
    # matplotlib table cells do not auto-wrap, so wrap the two wide text columns
    # manually to keep long comparison/function labels from being clipped.
    display["Comparison"] = display["Comparison"].map(lambda s: textwrap.fill(s, 22))
    display["Function"] = display["Function"].map(lambda s: textwrap.fill(s, 40))

    # Tallest wrapped cell in each row (in text lines) drives that row's height.
    row_lines = [
        max(str(c).count("\n") + 1 for c in display.iloc[i])
        for i in range(len(display))
    ]
    header_units = 1.6
    total_units = header_units + sum(n + 0.55 for n in row_lines)
    fig_h = 0.20 * total_units + 0.9     # title band + table body
    fig, ax = plt.subplots(figsize=(W_DOUBLE + 0.5, fig_h))
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.01)

    # bbox=[0,0,1,1] makes the table fill the (title-reserved) axes; per-row
    # heights below are RELATIVE — matplotlib rescales them to fit the bbox, so
    # multi-line wrapped cells are never clipped and rows never overflow.
    table = ax.table(
        cellText=display.values.tolist(),
        colLabels=list(display.columns),
        cellLoc="left",
        colWidths=[0.22, 0.07, 0.06, 0.09, 0.12, 0.44],
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(FONT["annot"] + 0.5)

    for j in range(len(display.columns)):
        table[0, j].set_height(header_units)
    for i, n_lines in enumerate(row_lines):
        for j in range(len(display.columns)):
            table[i + 1, j].set_height(n_lines + 0.55)

    header_fill = "#34495E"
    for j in range(len(display.columns)):
        cell = table[0, j]
        cell.set_facecolor(header_fill)
        cell.set_text_props(color="white", fontweight="bold")
        cell.set_edgecolor("white")

    for i in range(len(display)):
        bg = "#F5F7FA" if i % 2 else "white"
        depleted = tbl.loc[i, "Direction"] == "Depleted"
        for j, col in enumerate(display.columns):
            cell = table[i + 1, j]
            cell.set_facecolor(bg)
            cell.set_edgecolor("#E3E3E3")
            if col == "Direction":
                cell.set_text_props(color=COL_DEPLETED if depleted else COL_ENRICHED,
                                    fontweight="bold")
            elif col == "Gene":
                cell.set_text_props(fontweight="bold")

    ax.set_title(
        "Significant genes — GBM TRAIL metabolic CRISPR/Cas9 screen\n"
        f"($p$ < {ctx.pval_thresh}; FDR unreliable due to library bottleneck — "
        "results are EXPLORATORY)",
        fontsize=FONT["base"], fontweight="bold", pad=10,
    )
    save_figure(fig, "Fig9_summary_table", ctx.outdir, ctx.png_dpi)

    tbl.to_csv(ctx.outdir / "significant_genes.csv", index=False)


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

FIGURES = [
    (1, fig1_volcano),
    (2, fig2_rank),
    (3, fig3_heatmap),
    (4, fig4_barplot),
    (5, fig5_qc),
    (6, fig6_sgrna),
    (7, fig7_replicates),
    (8, fig8_cross),
    (9, fig9_table),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate publication-quality CRISPR figures.")
    parser.add_argument("--test-dir", type=Path,
                        default=Path(__file__).resolve().parent / "test")
    parser.add_argument("--count-file", type=Path,
                        default=Path(__file__).resolve().parent / "count" / "count_table_fixed.txt")
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).resolve().parent / "figures")
    parser.add_argument("--pval-threshold", type=float, default=0.05)
    parser.add_argument("--png-dpi", type=int, default=600, choices=[300, 600],
                        help="Raster PNG resolution for preview (PDF is always vector).")
    parser.add_argument("--figures-only", type=str, default=None,
                        help="Comma-separated figure numbers (e.g. '1,3,5').")
    args = parser.parse_args()

    set_publication_style()
    args.outdir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Publication-Quality CRISPR Screen Figures")
    logger.info("=" * 60)
    logger.info("Output : %s", args.outdir)
    logger.info("PNG dpi: %d   |   PDF: vector", args.png_dpi)
    logger.info("")

    # ---- DATA LOADING (analysis logic unchanged — delegated to _thesis_utils) ----
    data = load_all_comparisons(args.test_dir, args.pval_threshold)
    if not data:
        logger.error("No gene summaries found.")
        return

    counts = count_matrix = None
    sample_groups: dict[str, list[str]] = {}
    if args.count_file.exists():
        try:
            counts, count_matrix, sample_groups = load_count_table(args.count_file)
        except FileNotFoundError:
            logger.warning("Count table not found — some figures will be skipped.")

    all_sig = sorted(set().union(*[set(df[df["sig"]]["id"]) for df in data.values()]))
    logger.info("Significant genes: %d unique", len(all_sig))

    ctx = FigContext(
        data=data, counts=counts, count_matrix=count_matrix,
        sample_groups=sample_groups, sig_genes=all_sig,
        outdir=args.outdir, pval_thresh=args.pval_threshold, png_dpi=args.png_dpi,
    )

    only = None
    if args.figures_only:
        only = {int(x) for x in args.figures_only.split(",") if x.strip().isdigit()}

    for n, fn in FIGURES:
        if only is not None and n not in only:
            continue
        try:
            fn(ctx)
        except Exception as exc:  # never let one figure abort the run
            logger.error("Fig%d failed: %s", n, exc)
            import traceback
            traceback.print_exc()

    n_pdf = len(list(args.outdir.glob("*.pdf")))
    n_png = len(list(args.outdir.glob("*.png")))
    logger.info("")
    logger.info("Done — %d PDFs, %d PNGs in %s", n_pdf, n_png, args.outdir)

    export_metadata(args.outdir, step="05_generate_figures", extra={
        "n_figures": n_pdf,
        "n_significant_genes": len(all_sig),
        "pval_threshold": args.pval_threshold,
        "png_dpi": args.png_dpi,
    })


if __name__ == "__main__":
    main()
