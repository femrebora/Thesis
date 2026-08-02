# analysis/utils.R
# Shared utilities for all modules
# v2.0 — extended with publication-quality functions

library(ggplot2)

# ============================================================================
# Themes
# ============================================================================

# --- Legacy theme (kept for backward compatibility with modules 02-05) ---
theme_gbm <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
}

# --- Publication journal theme (colorblind-safe, clean) ---
theme_gbm_journal <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "#333333"),
      plot.title = element_text(
        face = "bold", size = base_size + 3, hjust = 0,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        size = base_size, hjust = 0, color = "#555555",
        margin = margin(b = 10)
      ),
      axis.title = element_text(size = base_size + 1),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      axis.text = element_text(size = base_size - 1, color = "#444444"),
      axis.ticks = element_line(color = "#cccccc"),
      legend.position = "bottom",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.box.spacing = unit(0, "mm"),
      legend.margin = margin(t = 0),
      panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#cccccc", fill = NA),
      strip.background = element_rect(fill = "#f5f5f5", color = "#cccccc"),
      strip.text = element_text(face = "bold", size = base_size - 1),
      plot.margin = margin(12, 12, 12, 12)
    )
}

# ============================================================================
# Color palettes (colorblind-safe, academic)
# ============================================================================

# Wong (2011) colorblind-safe palette — 7 colors
cb_palette <- c(
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#009E73",  # green
  "#F0E442",  # yellow
  "#CC79A7",  # reddish purple
  "#56B4E9",  # sky blue
  "#E69F00"   # orange
)

# Gene class colors
gene_class_palette <- c(
  "Depleted" = "#D55E00",  # vermillion
  "Enriched" = "#0072B2",  # blue
  "Other"    = "#999999"   # grey
)

# KM group colors
km_group_palette <- c(
  "High expression" = "#D55E00",
  "Low expression"  = "#0072B2"
)

# ============================================================================
# Figure saving
# ============================================================================

# --- Legacy save (kept for backward compatibility) ---
save_figure <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height, device = "pdf")
  ggsave(paste0(filename, ".tiff"), plot, width = width, height = height, dpi = dpi, device = "tiff", compression = "lzw")
  message("[save] ", filename, ".{pdf,tiff}")
}

# --- Publication save: PDF + PNG(600dpi) + TIFF(600dpi LZW) ---
save_figure_pub <- function(plot, filepath, width = 8, height = 6) {
  # Ensure output directory exists
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)

  # PDF (vector, best for publications). Use the string "pdf" device so
  # ggsave omits the raster-only `dpi` argument (cairo_pdf has no `dpi`
  # parameter and errors when ggplot2 >= 3.5 forwards it).
  ggsave(paste0(filepath, ".pdf"), plot,
         width = width, height = height,
         device = "pdf")

  # PNG (600 dpi)
  ggsave(paste0(filepath, ".png"), plot,
         width = width, height = height,
         dpi = 600, device = "png")

  # TIFF (600 dpi, LZW compression)
  ggsave(paste0(filepath, ".tiff"), plot,
         width = width, height = height,
         dpi = 600, device = "tiff", compression = "lzw")

  message(cli::cli_alert_success("{.file {basename(filepath)}} → .pdf, .png, .tiff"))
}

# --- Save KM plot (ggsurvplot object) to publication files ---
save_figure_km <- function(km_plot, filepath, width = 8, height = 7) {
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)

  # PDF
  pdf(paste0(filepath, ".pdf"), width = width, height = height)
  print(km_plot)
  dev.off()

  # PNG
  png(paste0(filepath, ".png"), width = width, height = height, units = "in", res = 600)
  print(km_plot)
  dev.off()

  # TIFF
  tiff(paste0(filepath, ".tiff"), width = width, height = height,
       units = "in", res = 600, compression = "lzw")
  print(km_plot)
  dev.off()

  message(cli::cli_alert_success("{.file {basename(filepath)}} → .pdf, .png, .tiff"))
}

# ============================================================================
# Table saving
# ============================================================================

# --- Legacy save (kept for backward compatibility) ---
save_table <- function(df, filename) {
  write.csv(df, paste0(filename, ".csv"), row.names = FALSE)
  message("[save] ", filename, ".csv")
}

# --- Publication save: CSV + TSV + XLSX ---
save_table_pub <- function(df, filepath) {
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)

  # CSV
  readr::write_csv(df, paste0(filepath, ".csv"))

  # TSV
  readr::write_tsv(df, paste0(filepath, ".tsv"))

  # XLSX
  openxlsx::write.xlsx(df, paste0(filepath, ".xlsx"))

  message(cli::cli_alert_success("{.file {basename(filepath)}} → .csv, .tsv, .xlsx"))
}

# ============================================================================
# Formatting helpers
# ============================================================================

# --- Format p-value for display ---
format_pval <- function(p) {
  ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p))
}

# --- Format HR with CI for plot annotation ---
format_hr_ci <- function(hr, lower, upper) {
  sprintf("HR = %.2f (95%% CI %.2f–%.2f)", hr, lower, upper)
}

# --- TCGA barcode: extract sample type code (positions 14-15) ---
parse_tcga_sample_type <- function(barcodes) {
  # TCGA barcode format: TCGA-XX-XXXX-XXY-ZZ-YYYY
  # Position 14-15: sample type code
  # 01 = Primary Solid Tumor, 11 = Solid Tissue Normal, etc.
  sample_code <- substr(barcodes, 14, 15)
  sample_type <- dplyr::case_when(
    sample_code == "01" ~ "Primary Tumor",
    sample_code == "02" ~ "Recurrent Tumor",
    sample_code == "06" ~ "Metastatic",
    sample_code == "11" ~ "Solid Tissue Normal",
    TRUE ~ paste0("Other (", sample_code, ")")
  )
  return(list(code = sample_code, type = sample_type))
}

# --- TCGA barcode: extract patient ID (first 12 chars) ---
parse_tcga_patient <- function(barcodes) {
  substr(barcodes, 1, 12)
}

message("[utils] Loaded theme_gbm, theme_gbm_journal, save_figure, save_figure_pub, save_figure_km, save_table, save_table_pub, format_pval, format_hr_ci, parse_tcga_sample_type, parse_tcga_patient")
