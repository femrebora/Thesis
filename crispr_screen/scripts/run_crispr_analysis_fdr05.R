#!/usr/bin/env Rscript
# ============================================================================
# run_crispr_analysis_fdr05.R — rebuild the full CRISPR figure set under a
# Benjamini-Hochberg FDR < 0.05 acceptance criterion.
#
# The thesis baseline calls a gene significant on its NOMINAL MAGeCK RRA
# p-value (p < 0.05, uncorrected) and shows enrichment terms ranked by adjusted
# p with no significance cutoff. This driver re-runs the identical figure
# scripts with the acceptance criterion tightened to FDR < 0.05 — the
# conventional cutoff in the CRISPR-screen literature — and writes every output
# with a "_fdr05" suffix so no baseline figure is overwritten.
#
# Expect null results at the gene level: in this screen no gene reaches
# FDR < 0.05 in any of the three comparisons (smallest observed FDR = 0.425).
# The figure scripts render that outcome as an explicitly labelled null-result
# panel rather than a blank canvas.
#
# Usage (from anywhere inside the repo):
#   Rscript crispr_screen/scripts/run_crispr_analysis_fdr05.R
# ============================================================================

suppressPackageStartupMessages(library(here))

Sys.setenv(
  THESIS_SIG_METRIC = "fdr",
  THESIS_SIG_THRESH = "0.05",
  THESIS_SIG_SUFFIX = "_fdr05"
)

source(here::here("crispr_screen", "scripts", "run_crispr_analysis.R"))
