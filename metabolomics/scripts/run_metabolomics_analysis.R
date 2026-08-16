#!/usr/bin/env Rscript
# ============================================================================
# run_metabolomics_analysis.R — Metabolomics module entry point
# A172-S vs A172-R LC-MS metabolomics profiling
#
# PROVENANCE / UNPUBLISHED / NOT INCLUDED IN THE FINAL SUBMITTED THESIS RESULTS.
# This module is retained to document the broader research process. Do not treat
# its plots as final thesis figures.
#
# This pipeline analyses 104 metabolites profiled in A172 sensitive (S, n=4)
# and resistant (R, n=4) GBM cell lines via LC-MS.
#
# Three sub-scripts (run in order):
#   1. metabolomics_figures.R       — Statistical analysis: normalization,
#                                     fold-change, t-tests, volcano plot,
#                                     boxplots, heatmap, PCA/PLS-DA
#   2. metabolomics_enrichment_pathway.R — Pathway enrichment (SMPDB/KEGG ORA)
#   3. metabolomics_integration_figures.R — CRISPR + metabolomics integration
#
# NOTE: PCA and PLS-DA outputs are exploratory within this provenance module.
# See metabolomics/upstream/README.md for known runner/network limitations.
# See docs/figure_table_mapping.md for the submitted-thesis figure scope.
#
# Usage:
#   Rscript metabolomics/scripts/run_metabolomics_analysis.R
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(cli)
})

cli::cli_h1("Metabolomics Analysis Pipeline")
cli::cli_text("Module: {.file {here::here('metabolomics')}}")
cli::cli_text("Comparison: A172-S (Sensitive, n=4) vs A172-R (Resistant, n=4)")
cli::cli_text("Metabolites: 104 (LC-MS)")

scripts <- c(
  "metabolomics_figures.R",
  "metabolomics_enrichment_pathway.R",
  "metabolomics_integration_figures.R"
)

for (scr in scripts) {
  cli::cli_h2("Running: {.file {scr}}")
  script_path <- here::here("metabolomics", "scripts", scr)

  if (!file.exists(script_path)) {
    cli::cli_alert_warning("{.file {scr}} not found — skipping")
    next
  }

  tryCatch(
    source(script_path),
    error = function(e) {
      cli::cli_alert_danger("{.file {scr}} failed: {e$message}")
    }
  )
}

cli::cli_h1("Metabolomics pipeline complete")
cli::cli_text("Outputs:")
cli::cli_li("{.file metabolomics/figures/} — All figures (PDF, PNG, TIFF)")
cli::cli_li("{.file metabolomics/data/} — Processed data files")
cli::cli_text("Note: PCA/PLS-DA figures are generated but excluded from thesis.")

writeLines(capture.output(sessionInfo()),
           file.path(here::here("metabolomics", "figures"),
                     "run_metabolomics_session_info.txt"))
