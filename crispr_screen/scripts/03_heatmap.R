#!/usr/bin/env Rscript
# Fig03 -- Karşılaştırmalı LFC haritası: anlamlı genlerin Log2 kat değişimi.
# Ölçülmeyen hücreler NA olarak işlenir ve BOŞ (çok açık gri) bırakılır; renk
# ölçeğine katkı vermez ve 0.0 olarak gösterilmez. Gerçek sıfır, simetrik ayrışan
# ölçeğin merkezidir. Satırlar 1. karşılaştırmanın KD'sine göre sıralanır.
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

pval_thresh <- 0.05

gene_data <- list()
all_genes <- character()
for (comp in names(COMPARISONS)) {
  df <- load_gene_summary(comp, pval_thresh)
  gene_data[[comp]] <- df
  all_genes <- union(all_genes, df$id[df$sig])
}

mat <- do.call(cbind, lapply(names(COMPARISONS), function(cn) {
  df <- gene_data[[cn]]
  setNames(df$LFC, df$id)[all_genes]
}))
rownames(mat) <- all_genes
colnames(mat) <- c(
  "Dirençli: Tümör / Gün 0",
  "Duyarlı / Dirençli (Tümör)",
  "Duyarlı: Tümör / Gün 0"
)

# Ölçülmeyen hücreler NA olarak KALIR (boş bırakılır); renk ölçeğine katkı vermez.

# Satırları 1. karşılaştırmanın (Dirençli: Tümör / Gün 0) KD'sine göre azalan
# sırada diz: turuncu (artan) üstte, mavi (azalan) altta. NA'lar sona iner.
mat <- mat[order(mat[, 1], decreasing = TRUE, na.last = TRUE), , drop = FALSE]

# Simetrik ayrışan ölçek, 0 merkezli: güçlü azalma (mavi) - nötr 0 (beyaz) -
# güçlü artış (turuncu). NA hücreler çok açık gri (C_NA) ve ölçeğe dahil değil.
max_abs <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs),
                      c(C_DEPLETED, "white", C_ENRICHED))

ht <- Heatmap(
  mat,
  name = "Log2KD",
  col = col_fun,
  na_col = C_NA,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat[i, j]
    if (is.na(v)) return(invisible(NULL))   # ölçülmeyen hücre: boş bırak
    grid.text(sprintf("%.1f", v), x, y,
              gp = gpar(fontsize = 7.5, fontfamily = "sans",
                        col = ifelse(abs(v) > max_abs * 0.5, "white", "#1A1A1A"),
                        fontface = "bold"))
  },
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontfamily = "sans", fontface = "italic"),
  column_names_side = "bottom",
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 9, fontfamily = "sans", fontface = "bold"),
  border = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(
    title = expression(Log[2]~"LFC"),
    title_gp = gpar(fontsize = 9, fontfamily = "sans", fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontfamily = "sans"),
    legend_height = unit(3, "cm")
  )
)

draw_ht <- function() {
  draw(ht,
       gap = unit(c(2, 5), "mm"),
       padding = unit(c(20, 4, 8, 4), "mm"))
  # Açıklayıcı dipnot: boş hücreler ilgili karşılaştırmada ölçüm olmadığını gösterir.
  grid.text("Boş hücre: ölçülmedi",
            x = unit(2, "mm"), y = unit(2.5, "mm"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 8, col = "#666666", fontfamily = "sans"))
}

save_grid_figure(draw_ht, "Fig03_heatmap",
                 width = 7, height = max(5, nrow(mat) * 0.32 + 3.5))
