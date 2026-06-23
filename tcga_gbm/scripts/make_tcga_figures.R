#!/usr/bin/env Rscript
# ============================================================================
# make_tcga_figures.R — TCGA-GBM thesis figures (Turkish labels)
# 12-gene expression & survival analysis for Şekil 4.15–4.17
#
# Input:  tcga_gbm/results/01_expression_survival/tables/*.csv (precomputed)
#         tcga_gbm/cache/GDCdata/ (STAR-counts TSVs, optional)
# Output: tcga_gbm/figures/tcga_*.{pdf,png,tiff} (600 dpi)
#
# Usage:  Rscript tcga_gbm/scripts/make_tcga_figures.R
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
})

# Source shared utilities (config, themes, palettes, helpers)
source(here::here("tcga_gbm", "R", "utils.R"))

# ---- Paths (all anchored via here::here()) ----
TABLES    <- TABLES_DIR_01
OUT       <- FIGURES_DIR
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DPI       <- 600

# ---- Turkish gene classification (for figure labels) ----
ALL_GENES_TR <- ALL_GENES
GENE_CLASS_TR <- c(
  setNames(rep("Tükenmiş",     length(GENES_DEPLETED)), GENES_DEPLETED),
  setNames(rep("Zenginleşmiş", length(GENES_ENRICHED)), GENES_ENRICHED),
  setNames(rep("Diğer",        length(GENES_OTHER)),    GENES_OTHER)
)

# ---- Theme ----
theme_tez <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.3, color = "grey90"),
    axis.title        = element_text(size = 13, face = "bold"),
    axis.text         = element_text(size = 11, color = "black"),
    legend.title      = element_text(face = "bold"),
    legend.position   = "bottom"
  )

# ---- Save helper ----
save_fig <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(OUT, paste0("tcga_", name, ".pdf")),  p,
         width = w, height = h, device = "pdf")
  ggsave(file.path(OUT, paste0("tcga_", name, ".png")),  p,
         width = w, height = h, dpi = DPI)
  ggsave(file.path(OUT, paste0("tcga_", name, ".tiff")), p,
         width = w, height = h, dpi = DPI, compression = "lzw")
  cat(sprintf("  ✔ tcga_%s.{pdf,png,tiff}\n", name))
}

# ============================================================================
# 1. COX FOREST PLOT — Median-split Cox PH (Şekil 4.16)
# ============================================================================
cat("\n--- 1. Cox forest plot ---\n")

surv <- read.csv(file.path(TABLES, "Table02_survival_summary.csv"),
                 check.names = FALSE)
surv <- surv[order(surv$cox_HR), ]
surv$gene  <- factor(surv$gene, levels = surv$gene)
surv$sinif <- GENE_CLASS_TR[as.character(surv$gene)]

p_forest <- ggplot(surv, aes(x = cox_HR, y = gene, color = sinif)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "#999999", linewidth = 0.6) +
  geom_point(aes(size = n_patients), alpha = 0.9) +
  geom_errorbarh(aes(xmin = cox_CI_lower, xmax = cox_CI_upper),
                 height = 0.2, linewidth = 0.8) +
  scale_x_log10(
    breaks = c(0.5, 0.7, 1.0, 1.3, 1.7),
    labels = c("0,5", "0,7", "1,0", "1,3", "1,7"),
    limits = c(0.55, 1.9)
  ) +
  scale_color_manual(values = class_colors_tr, name = "Gen sınıfı") +
  scale_size_continuous(range = c(3, 6), guide = "none") +
  labs(x = "HR (log ölçeği)", y = NULL) +
  theme_tez +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  ) +
  geom_text(
    aes(label = sprintf("p = %s", format_pval_tr(cox_p)),
        x     = cox_CI_upper),
    hjust = -0.3, size = 3.2, color = "#444444"
  )

save_fig(p_forest, "cox_orman", w = 10.5, h = 6.5)

# ============================================================================
# 2. CONTINUOUS COX FOREST PLOT — Robustness analysis (Şekil 4.17)
# ============================================================================
cat("\n--- 2. Continuous Cox forest plot ---\n")

rob <- read.csv(file.path(TABLES, "Table04_survival_robustness.csv"),
                check.names = FALSE)
rob <- rob[!is.na(rob$cont_HR_per_SD), ]
rob <- rob[order(rob$cont_HR_per_SD), ]
rob$gene  <- factor(rob$gene, levels = rob$gene)
rob$sinif <- GENE_CLASS_TR[as.character(rob$gene)]

# FDR significance star
rob$yildiz <- ifelse(rob$cont_adj_p < 0.05, "*", "")

p_cont <- ggplot(rob, aes(x = cont_HR_per_SD, y = gene, color = sinif)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "#999999", linewidth = 0.6) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_errorbarh(aes(xmin = cont_CI_SD_lower, xmax = cont_CI_SD_upper),
                 height = 0.2, linewidth = 0.8) +
  scale_x_log10() +
  scale_color_manual(values = class_colors_tr, name = "Gen sınıfı") +
  labs(x = "HR (1 SS ifade artışı başına, log ölçeği)", y = NULL) +
  theme_tez +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  ) +
  geom_text(
    aes(label = sprintf("p = %s%s", format_pval_tr(cont_p), yildiz),
        x     = cont_CI_SD_upper),
    hjust = -0.3, size = 3.2, color = "#444444",
    fontface = ifelse(rob$cont_adj_p < 0.05, 2, 1)
  )

save_fig(p_cont, "surekli_cox_orman", w = 10.5, h = 6.5)

# ============================================================================
# 3. EXPRESSION BOXPLOT — From GDCdata cache (Şekil 4.15)
# ============================================================================
cat("\n--- 3. Expression boxplot ---\n")

build_expr_from_cache <- function(cache_dir, genes) {
  dirs <- list.dirs(cache_dir, full.names = TRUE, recursive = FALSE)
  if (length(dirs) == 0) stop("GDCdata cache directory is empty")

  cat(sprintf("  Scanning %d sample directories...\n", length(dirs)))

  expr_list <- list()
  barcodes  <- character()

  for (i in seq_along(dirs)) {
    tsv_files <- list.files(dirs[i], pattern = "\\.tsv$", full.names = TRUE)
    if (length(tsv_files) == 0) next

    dat <- tryCatch(
      read.delim(tsv_files[1], header = TRUE, comment.char = "#",
                 stringsAsFactors = FALSE, quote = ""),
      error = function(e) NULL
    )
    if (is.null(dat)) next
    if (!all(c("gene_name", "tpm_unstranded") %in% colnames(dat))) next

    idx <- which(dat$gene_name %in% genes)
    if (length(idx) == 0) next

    vals <- setNames(dat$tpm_unstranded[idx], dat$gene_name[idx])
    expr_list[[length(expr_list) + 1]] <- vals
    barcodes <- c(barcodes, basename(dirs[i]))
  }

  if (length(expr_list) == 0) stop("No target genes found in any sample")

  all_genes   <- unique(unlist(lapply(expr_list, names)))
  found_genes <- intersect(genes, all_genes)

  mat <- matrix(NA_real_, nrow = length(found_genes), ncol = length(expr_list))
  rownames(mat) <- found_genes
  colnames(mat) <- barcodes

  for (j in seq_along(expr_list)) {
    for (g in found_genes) {
      if (g %in% names(expr_list[[j]])) {
        mat[g, j] <- expr_list[[j]][g]
      }
    }
  }

  # log2(TPM + 1) transformation
  mat_log2 <- log2(mat + 1)
  cat(sprintf("  Expression matrix: %d genes × %d samples\n",
              nrow(mat_log2), ncol(mat_log2)))
  cat(sprintf("  Found genes: %s\n",
              paste(found_genes, collapse = ", ")))

  missing <- setdiff(genes, found_genes)
  if (length(missing) > 0) {
    cat(sprintf("  Missing genes: %s\n",
                paste(missing, collapse = ", ")))
  }

  mat_log2
}

expr_try <- tryCatch(
  build_expr_from_cache(CACHE_DIR, ALL_GENES),
  error = function(e) { cat(sprintf("  WARNING: %s\n", e$message)); NULL }
)

if (!is.null(expr_try) && nrow(expr_try) >= 8) {
  genes      <- rownames(expr_try)
  n_patients <- ncol(expr_try)

  expr_long <- data.frame(
    gen   = factor(rep(genes, each = n_patients), levels = genes),
    hasta = rep(colnames(expr_try), times = length(genes)),
    ifade = as.vector(t(expr_try)),
    sinif = rep(GENE_CLASS_TR[genes], each = n_patients),
    stringsAsFactors = FALSE
  )

  p_box <- ggplot(expr_long, aes(x = gen, y = ifade, fill = sinif)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, width = 0.65) +
    scale_fill_manual(values = class_colors_tr, name = "Gen sınıfı") +
    labs(x = NULL, y = expression(log[2] * "(TPM + 1) Expression Level")) +
    theme_tez +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  save_fig(p_box, "ifade_kutu", w = 10, h = 6)

} else {
  cat("  Expression boxplot skipped (insufficient data). "
      "Falling back to summary bar plot from Table01...\n")

  expr_sum <- read.csv(file.path(TABLES, "Table01_expression_gene_summary.csv"),
                       check.names = FALSE)
  expr_sum <- expr_sum[order(-expr_sum$mean_expression), ]
  expr_sum$gene  <- factor(expr_sum$gene, levels = expr_sum$gene)
  expr_sum$sinif <- GENE_CLASS_TR[as.character(expr_sum$gene)]

  p_bar <- ggplot(expr_sum, aes(x = gene, y = mean_expression, fill = sinif)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = mean_expression - sd_expression,
                      ymax = mean_expression + sd_expression),
                  width = 0.25, linewidth = 0.5) +
    scale_fill_manual(values = class_colors_tr, name = "Gen sınıfı") +
    labs(x = NULL,
         y = expression(log[2] * "(TPM + 1) Mean Expression")) +
    theme_tez +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  save_fig(p_bar, "ifade_cubuk", w = 10, h = 6)
}

# ============================================================================
# 4. SURVIVAL SUMMARY TABLE — Turkish CSV
# ============================================================================
cat("\n--- 4. Survival summary table (Turkish) ---\n")

tablo <- surv[, c("gene", "gene_class", "n_patients", "n_events",
                   "median_OS_high_months", "median_OS_low_months",
                   "cox_HR", "cox_CI_lower", "cox_CI_upper",
                   "cox_p", "adjusted_p", "PH_assumption_p")]
names(tablo) <- c("Gen", "Sınıf", "Hasta_sayısı", "Olay_sayısı",
                  "Medyan_OS_yüksek_ifade_ay", "Medyan_OS_düşük_ifade_ay",
                  "HR", "HR_GA_alt", "HR_GA_üst",
                  "p_değeri", "Düzeltilmiş_p_FDR", "PH_varsayımı_p")

tablo$Sınıf <- GENE_CLASS_TR[as.character(tablo$Gen)]

# Interpretation column
tablo$Yorum <- ifelse(
  tablo$p_değeri < 0.05 & tablo$HR > 1,
  "Yüksek ifade → daha kötü sağkalım (nominal)",
  ifelse(tablo$p_değeri < 0.05 & tablo$HR < 1,
         "Yüksek ifade → daha iyi sağkalım (nominal)",
         "Anlamlı ilişki yok")
)
tablo$Yorum[tablo$PH_varsayımı_p < 0.05] <- paste(
  tablo$Yorum[tablo$PH_varsayımı_p < 0.05],
  "[PH varsayımı ihlal edildi]"
)

write.csv(tablo, file.path(OUT, "tcga_sagkalim_ozet_tablo.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("  ✔ tcga_sagkalim_ozet_tablo.csv\n")

# ============================================================================
# 5. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("TCGA figures complete.\n")
cat("Output directory:", OUT, "\n")
cat("Files created:\n")
cat("  tcga_cox_orman.{pdf,png,tiff}         — Cox forest plot (Şekil 4.16)\n")
cat("  tcga_surekli_cox_orman.{pdf,png,tiff} — Continuous Cox robustness (Şekil 4.17)\n")
if (!is.null(expr_try) && nrow(expr_try) >= 8) {
  cat("  tcga_ifade_kutu.{pdf,png,tiff}       — Expression boxplot (Şekil 4.15)\n")
} else {
  cat("  tcga_ifade_cubuk.{pdf,png,tiff}      — Expression bar plot (summary fallback)\n")
}
cat("  tcga_sagkalim_ozet_tablo.csv          — Survival summary table (Turkish)\n")
cat("========================================\n")

# Record package versions for reproducibility
cat("\n--- Session Info ---\n")
writeLines(capture.output(sessionInfo()),
           file.path(OUT, "tcga_figures_session_info.txt"))
