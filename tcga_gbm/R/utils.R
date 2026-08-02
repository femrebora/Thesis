# utils.R -- shared configuration, themes, palettes, and helpers for TCGA-GBM module
# Thesis §3.10 / §4.x — 12-gene TCGA expression & survival analysis

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(cli)
})

set.seed(42)

# -- Paths (here::here() anchored) -----------------------------------------------
.MODULE_DIR    <- here::here("tcga_gbm")
RESULTS_DIR    <- file.path(.MODULE_DIR, "results")
FIGURES_DIR    <- file.path(.MODULE_DIR, "figures")
CACHE_DIR      <- file.path(.MODULE_DIR, "cache", "GDCdata")
OUTPUT_DIR_01  <- file.path(RESULTS_DIR, "01_expression_survival")
FIGURES_DIR_01 <- file.path(OUTPUT_DIR_01, "figures")
TABLES_DIR_01  <- file.path(OUTPUT_DIR_01, "tables")

dir.create(CACHE_DIR,      showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR_01, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_DIR_01,  showWarnings = FALSE, recursive = TRUE)

# -- Gene lists (canonical — thesis Table 4.2) -----------------------------------
GENES_DEPLETED <- c("ASL", "GBA3", "ATP6V0A4", "TPCN1")
GENES_ENRICHED <- c("PLA2G4E", "SLC1A1", "RRM2", "ARSD", "GRIN1")
GENES_OTHER    <- c("NAT2", "SLC25A20", "ALOX15")
ALL_GENES      <- c(GENES_DEPLETED, GENES_ENRICHED, GENES_OTHER)

# -- Classifications (English, for data processing) -------------------------------
GENE_CLASS <- c(
  setNames(rep("Depleted", length(GENES_DEPLETED)), GENES_DEPLETED),
  setNames(rep("Enriched", length(GENES_ENRICHED)), GENES_ENRICHED),
  setNames(rep("Other",    length(GENES_OTHER)),    GENES_OTHER)
)

# -- Turkish labels (for figure output) -------------------------------------------
GENE_CLASS_TR <- c(
  setNames(rep("Azalmış", length(GENES_DEPLETED)), GENES_DEPLETED),
  setNames(rep("Artmış",  length(GENES_ENRICHED)), GENES_ENRICHED),
  setNames(rep("Diğer",         length(GENES_OTHER)),    GENES_OTHER)
)

# -- LFC values from thesis Table 4.2 (CRISPR Rpost vs R0) -----------------------
GENE_LFC_RESISTANT <- c(
  "PLA2G4E"   = 12.3,  "SLC1A1"  = 11.4,  "RRM2"     = 10.9,
  "ARSD"      = 10.4,  "GRIN1"   = 10.3,  "NAT2"     = 3.2,
  "SLC25A20"  = -0.0,  "TPCN1"   = -6.7,  "ATP6V0A4" = -7.3,
  "GBA3"      = -8.8,  "ASL"     = -11.5, "ALOX15"   = NA
)

# -- TCGA Project ----------------------------------------------------------------
TCGA_PROJECT <- "TCGA-GBM"
USE_CACHE    <- TRUE

# -- Analysis thresholds ---------------------------------------------------------
MIN_PATIENTS_KM      <- 20   # Minimum matched patients for survival analysis
MIN_NORMAL_SAMPLES   <- 5    # Minimum normal samples for tumor-normal comparison
MIN_GENES_HEATMAP    <- 8    # Minimum genes found for expression heatmap
MIN_EVENTS_PER_GROUP <- 3    # Minimum events per KM group

# -- Output formats --------------------------------------------------------------
FIGURE_FORMATS <- c("pdf", "png", "tiff")
FIGURE_DPI_PNG  <- 600
FIGURE_DPI_TIFF <- 600
TABLE_FORMATS   <- c("csv", "tsv", "xlsx")

# ============================================================================
# Themes
# ============================================================================

# Publication journal theme (colorblind-safe, clean)
theme_gbm_journal <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text             = element_text(family = "sans", color = "#333333"),
      plot.title       = element_text(face = "bold", size = base_size + 3,
                                      hjust = 0, margin = margin(b = 8)),
      plot.subtitle    = element_text(size = base_size, hjust = 0,
                                      color = "#555555", margin = margin(b = 10)),
      axis.title       = element_text(size = base_size + 1),
      axis.title.x     = element_text(margin = margin(t = 6)),
      axis.title.y     = element_text(margin = margin(r = 6)),
      axis.text        = element_text(size = base_size - 1, color = "#444444"),
      axis.ticks       = element_line(color = "#cccccc"),
      legend.position  = "bottom",
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.box.spacing = unit(0, "mm"),
      legend.margin    = margin(t = 0),
      panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "#cccccc", fill = NA),
      strip.background = element_rect(fill = "#f5f5f5", color = "#cccccc"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = margin(12, 12, 12, 12)
    )
}

# ============================================================================
# Color palettes (colorblind-safe, academic)
# ============================================================================

# Wong (2011) 7-color colorblind-safe palette
cb_palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#F0E442",
  "#CC79A7", "#56B4E9", "#E69F00"
)

# Gene class colors (TCGA convention)
# NOTE: TCGA palette convention is INVERTED relative to CRISPR module.
#   TCGA:     Enriched = blue (#0072B2), Depleted = vermillion (#D55E00)
#   CRISPR:   Enriched = orange (#D95F02), Depleted = blue (#2C7FB8)
# This reflects the different biological contexts: in CRISPR, "Enriched" means
# sgRNA representation increased (depletion → enriched in resistant tumours),
# while in TCGA "Enriched" means higher mRNA expression.
# Both conventions are documented in docs/methods_and_thresholds.md.
gene_class_palette <- c(
  "Depleted" = "#D55E00",  # vermillion
  "Enriched" = "#0072B2",  # blue
  "Other"    = "#999999"   # grey
)

# Turkish-label palette (matches gene_class_palette values).
# Keys pinned to the thesis legend wording — see make_tcga_figures.R.
class_colors_tr <- c(
  "Artmış"  = "#0072B2",  # blue
  "Azalmış" = "#D55E00",  # vermillion
  "Diğer"   = "#999999"   # grey
)

# KM group colors
km_group_palette <- c(
  "High expression" = "#D55E00",
  "Low expression"  = "#0072B2"
)

# ============================================================================
# Figure saving
# ============================================================================

# Publication save: PDF + PNG(600dpi) + TIFF(600dpi LZW)
save_figure_pub <- function(plot, filepath, width = 8, height = 6) {
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)

  ggsave(paste0(filepath, ".pdf"), plot,
         width = width, height = height, device = "pdf")
  ggsave(paste0(filepath, ".png"), plot,
         width = width, height = height, dpi = 600, device = "png")
  ggsave(paste0(filepath, ".tiff"), plot,
         width = width, height = height, dpi = 600,
         device = "tiff", compression = "lzw")

  cli::cli_alert_success("{.file {basename(filepath)}} → .pdf, .png, .tiff")
}

# Save KM plot (ggsurvplot object) to publication files
save_figure_km <- function(km_plot, filepath, width = 8, height = 7) {
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)

  pdf(paste0(filepath, ".pdf"), width = width, height = height)
  print(km_plot)
  dev.off()

  png(paste0(filepath, ".png"), width = width, height = height,
      units = "in", res = 600)
  print(km_plot)
  dev.off()

  tiff(paste0(filepath, ".tiff"), width = width, height = height,
       units = "in", res = 600, compression = "lzw")
  print(km_plot)
  dev.off()

  cli::cli_alert_success("{.file {basename(filepath)}} → .pdf, .png, .tiff")
}

# ============================================================================
# Table saving
# ============================================================================

save_table_pub <- function(df, filepath) {
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, paste0(filepath, ".csv"))
  readr::write_tsv(df, paste0(filepath, ".tsv"))
  openxlsx::write.xlsx(df, paste0(filepath, ".xlsx"))
  cli::cli_alert_success("{.file {basename(filepath)}} → .csv, .tsv, .xlsx")
}

# ============================================================================
# Formatting helpers
# ============================================================================

format_pval <- function(p) {
  ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p))
}

# Turkish decimal: period → comma
format_pval_tr <- function(p) {
  x <- ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p))
  gsub("\\.", ",", x)
}

format_hr_ci <- function(hr, lower, upper) {
  sprintf("HR = %.2f (95%% CI %.2f–%.2f)", hr, lower, upper)
}

# TCGA barcode: extract sample type code (positions 14-15)
# 01 = Primary Solid Tumor, 11 = Solid Tissue Normal, etc.
parse_tcga_sample_type <- function(barcodes) {
  sample_code <- substr(barcodes, 14, 15)
  sample_type <- dplyr::case_when(
    sample_code == "01" ~ "Primary Tumor",
    sample_code == "02" ~ "Recurrent Tumor",
    sample_code == "06" ~ "Metastatic",
    sample_code == "11" ~ "Solid Tissue Normal",
    TRUE                ~ paste0("Other (", sample_code, ")")
  )
  return(list(code = sample_code, type = sample_type))
}

# TCGA barcode: extract patient ID (first 12 chars)
parse_tcga_patient <- function(barcodes) {
  substr(barcodes, 1, 12)
}

message("[tcga_gbm/utils] Loaded ", length(ALL_GENES), " genes ",
        "(Depleted=", length(GENES_DEPLETED),
        " Enriched=", length(GENES_ENRICHED),
        " Other=", length(GENES_OTHER), ")")
message("[tcga_gbm/utils] Module dir: ", .MODULE_DIR)
