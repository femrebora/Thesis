#!/usr/bin/env Rscript
# Fig09 -- Top 12 most significant genes ranked by effect size (TÜSEB table style).
# Wide layout: gene x comparison Log2 FC matrix. Unmeasured cells are blank.
# (never left empty or border-only). FC values are coloured by sign:
# enriched (orange) / depleted (blue) / blank (grey).
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(gridExtra)
  library(grid)
})

pval_thresh <- 0.05
top_n       <- 12

# Her karşılaştırma için tam gen tablosu (LFC + anlamlılık).
gene_data <- lapply(names(COMPARISONS), function(cn) load_gene_summary(cn, pval_thresh))
names(gene_data) <- names(COMPARISONS)

# Per-gene minimum p-value across comparisons (ranking + display).
gene_minp <- bind_rows(lapply(gene_data, function(d) d %>% select(id, pval, sig))) %>%
  group_by(id) %>%
  summarise(min_pval = min(pval), any_sig = any(sig), .groups = "drop")

# Top 12 most significant genes by minimum p-value.
top_genes <- gene_minp %>% filter(any_sig) %>% arrange(min_pval) %>%
  slice_head(n = top_n) %>% pull(id)

# Wide Log2 FC matrix -- unmeasured cells left blank.
lfc_wide <- sapply(names(COMPARISONS), function(cn) {
  v <- setNames(gene_data[[cn]]$LFC, gene_data[[cn]]$id)[top_genes]
  v
})
lfc_wide <- matrix(lfc_wide, nrow = length(top_genes),
                   dimnames = list(top_genes, names(COMPARISONS)))

# Sort rows by descending LFC in the first comparison: enriched on top, depleted below
# (same ordering as the heatmap).
ord       <- order(lfc_wide[, 1], decreasing = TRUE)
top_genes <- top_genes[ord]
lfc_wide  <- lfc_wide[ord, , drop = FALSE]

min_pval_top <- gene_minp$min_pval[match(top_genes, gene_minp$id)]

# -- Display data frame ---------------------------------------------------------
col_headers <- vapply(names(COMPARISONS),
                      function(cn) wrap_labels(COMPARISONS[cn], 14), character(1))

disp <- data.frame(Gen = top_genes, check.names = FALSE, stringsAsFactors = FALSE)
for (cn in names(COMPARISONS)) {
  disp[[wrap_labels(COMPARISONS[cn], 14)]] <- sprintf("%.1f", lfc_wide[, cn])
  disp[[wrap_labels(COMPARISONS[cn], 14)]][is.na(lfc_wide[, cn])] <- ""
}
disp[["p-değeri"]] <- formatC(min_pval_top, format = "e", digits = 1)
disp[["İşlev"]]    <- wrap_labels(vapply(top_genes, get_annotation_tr, character(1)), 34)

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
write_csv(disp, file.path(TAB_DIR, "top12_significant_genes.csv"))

# Editable Excel export (v3)
tryCatch({
  suppressPackageStartupMessages(library(openxlsx))
  wb <- createWorkbook()
  addWorksheet(wb, "Oncelikli 12 Aday Gen")
  writeData(wb, sheet = 1, disp)
  setColWidths(wb, sheet = 1, cols = 1:ncol(disp), widths = "auto")
  headerStyle <- createStyle(textDecoration = "bold")
  addStyle(wb, sheet = 1, headerStyle, rows = 1, cols = 1:ncol(disp), gridExpand = TRUE)
  saveWorkbook(wb, file.path(FIG_DIR, "Fig09_summary_table_revised_v3.xlsx"),
               overwrite = TRUE)
  message("Fig09 Excel saved")
}, error = function(e) message("Excel export skipped: ", conditionMessage(e)))

# -- tableGrob -----------------------------------------------------------------
th <- ttheme_minimal(
  base_size = 8.5, base_family = "sans",
  core = list(
    fg_params = list(col = "#1A1A1A", fontfamily = "sans"),
    bg_params = list(fill = rep(c("white", "#F4F6F9"), length.out = nrow(disp)))
  ),
  colhead = list(
    fg_params = list(fontface = "bold", col = "#1A1A1A", fontfamily = "sans"),
    bg_params = list(fill = "#E3E8EF")
  )
)

tbl <- tableGrob(disp, rows = NULL, theme = th)

# Gen adlarını italik yap (1. sütun, gövde satırları).
gene_cells <- which(grepl("core-fg", tbl$layout$name) & tbl$layout$l == 1)
for (gi in gene_cells) {
  tbl$grobs[[gi]]$gp <- gpar(fontface = "italic", fontfamily = "sans",
                             fontsize = 8.5, col = "#1A1A1A")
}

# KD sütunlarını işarete göre renklendir (turuncu / mavi / gri) + kalın.
for (cj in seq_along(COMPARISONS)) {
  body_col <- cj + 1                     # 1. sütun = Gen
  for (ri in seq_len(nrow(disp))) {
    v  <- lfc_wide[ri, cj]
    if (is.na(v)) next
    cc <- if (v > 0) C_ENRICHED else if (v < 0) C_DEPLETED else "#999999"
    cell <- which(grepl("core-fg", tbl$layout$name) &
                  tbl$layout$t == (ri + 1) & tbl$layout$l == body_col)
    if (length(cell)) {
      tbl$grobs[[cell]]$gp <- gpar(col = cc, fontface = "bold",
                                   fontfamily = "sans", fontsize = 8.5)
    }
  }
}

# Başlık altına ince ayraç çizgisi.
tbl <- gtable::gtable_add_grob(
  tbl,
  grobs = grid::segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                             gp = grid::gpar(lwd = 1.2, col = "#888888")),
  t = 1, b = 1, l = 1, r = ncol(tbl))

title <- grid::textGrob(
  "Öncelikli 12 aday gen",
  gp = grid::gpar(fontsize = 12, fontface = "bold", fontfamily = "sans"),
  x = 0.01, hjust = 0)

caption <- grid::textGrob(
  "Değerler Log2 kat değişimini göstermektedir; boş hücreler ölçülmeyen karşılaştırmaları belirtir.",
  gp = grid::gpar(fontsize = 8, col = "#666666", fontfamily = "sans"),
  x = 0.01, hjust = 0)

arranged <- gridExtra::arrangeGrob(
  tbl, top = title, bottom = caption, padding = unit(0.6, "line"))

n_rows <- nrow(disp)
save_grid_figure(function() { grid::grid.newpage(); grid::grid.draw(arranged) },
                 "Fig09_summary_table",
                 width = 10, height = max(5, n_rows * 0.46 + 2))

# DOCX export for manuscript (guarded; never stops the script).
tryCatch({
  suppressPackageStartupMessages(library(flextable))
  ft <- flextable(disp) |> fontsize(size = 9, part = "all") |>
    bold(part = "header") |> autofit() |>
    set_caption("Öncelikli 12 aday gen")
  save_as_docx(ft, path = file.path(FIG_DIR, "Fig09_summary_table.docx"))
}, error = function(e) message("DOCX export skipped: ", conditionMessage(e)))

message("Fig09_summary_table saved")
