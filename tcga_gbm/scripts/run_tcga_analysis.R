#!/usr/bin/env Rscript
# ============================================================================
# run_tcga_analysis.R — TCGA-GBM module entry point
# Thesis §3.10 / §4.x — 12-gene expression & survival analysis in TCGA-GBM
#
# This is a two-phase pipeline:
#   Phase 1 — 01_expression_survival.R
#     Downloads TCGA-GBM RNA-seq data via TCGAbiolinks, prepares expression
#     matrices, performs median-split Cox PH survival analysis for the 12-gene
#     panel, runs robustness checks (continuous Cox, multivariable Cox), and
#     exports summary tables to tcga_gbm/results/01_expression_survival/tables/.
#
#   Phase 2 — make_tcga_figures.R
#     Reads the precomputed survival tables and optional GDCdata cache, generates
#     publication-quality Şekil 4.13–4.15 (median-split Cox forest / HR panel and
#     continuous Cox sensitivity forest) with Turkish labels, plus an expression
#     summary plot that is generated but not assigned a Şekil number in the
#     canonical thesis map (docs/figure_table_mapping.md).
#
# Usage:
#   Rscript tcga_gbm/scripts/run_tcga_analysis.R          # both phases
#   Rscript tcga_gbm/scripts/run_tcga_analysis.R phase1   # data acquisition only
#   Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2   # figures only
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(cli)
})

cli::cli_h1("TCGA-GBM Analysis Pipeline")
cli::cli_text("Module: {.file {here::here('tcga_gbm')}}")

args <- commandArgs(trailingOnly = TRUE)
run_phase <- if (length(args) > 0) tolower(args[1]) else "all"

# ---- Phase 1: Data acquisition & survival analysis ----
if (run_phase %in% c("all", "phase1")) {
  cli::cli_h2("Phase 1: Expression & Survival Analysis")

  phase1_script <- here::here("tcga_gbm", "scripts", "01_expression_survival.R")
  if (file.exists(phase1_script)) {
    cli::cli_alert_info("Running {.file 01_expression_survival.R}...")
    source(phase1_script)
  } else {
    cli::cli_alert_warning(
      "{.file 01_expression_survival.R} not found — ",
      "this script is a 1479-line monolith pending decomposition. ",
      "The precomputed survival tables in ",
      "{.file tcga_gbm/results/01_expression_survival/tables/} ",
      "should already be present."
    )
    cli::cli_alert_info(
      "To regenerate from source, run the original script at ",
      "{.file Results/analysis/01_expression_survival.R} ",
      "in the source workspace."
    )
  }
}

# ---- Phase 2: Figure generation ----
if (run_phase %in% c("all", "phase2")) {
  cli::cli_h2("Phase 2: Figure Generation")

  fig_script <- here::here("tcga_gbm", "scripts", "make_tcga_figures.R")
  if (file.exists(fig_script)) {
    cli::cli_alert_info("Running {.file make_tcga_figures.R}...")
    source(fig_script)
  } else {
    cli::cli_abort("{.file make_tcga_figures.R} not found")
  }
}

cli::cli_h1("TCGA-GBM pipeline complete")
cli::cli_text("Outputs:")
cli::cli_li("{.file tcga_gbm/figures/tcga_*.pdf/png/tiff}")
cli::cli_li("{.file tcga_gbm/results/01_expression_survival/tables/}")
cli::cli_li("{.file tcga_gbm/figures/tcga_sagkalim_ozet_tablo.csv}")

writeLines(capture.output(sessionInfo()),
           file.path(here::here("tcga_gbm", "figures"),
                     "run_tcga_session_info.txt"))
