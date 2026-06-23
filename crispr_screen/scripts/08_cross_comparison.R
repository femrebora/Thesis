#!/usr/bin/env Rscript
# Fig08 -- Compact correlation heatmap: pairwise Pearson r of gene LFCs across comparisons.
source(here::here("crispr_screen", "R", "utils.R"))

pval_thresh <- 0.05
comp_names  <- names(COMPARISONS)
n <- length(comp_names)

gene_data <- lapply(comp_names, load_gene_summary, pval_thresh = pval_thresh)
names(gene_data) <- comp_names
lfc_vecs <- lapply(gene_data, function(df) setNames(df$LFC, df$id))

short <- wrap_labels(COMPARISONS, width = 14)
mat <- matrix(NA, n, n, dimnames = list(short, short))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) { mat[i, j] <- 1; next }
    common <- intersect(names(lfc_vecs[[i]]), names(lfc_vecs[[j]]))
    if (length(common) < 3) next
    mat[i, j] <- cor(lfc_vecs[[i]][common], lfc_vecs[[j]][common], method = "pearson")
  }
}

plot_df <- as.data.frame(as.table(mat)) %>%
  filter(!is.na(Freq)) %>%
  rename(comparison1 = Var1, comparison2 = Var2, pearson_r = Freq)

# Compact tile-based correlation heatmap (geom_tile, not bubble).
# Text colour: white on strong correlations (|r| > 0.5), dark otherwise.
plot_df <- plot_df %>%
  mutate(
    text_col  = ifelse(abs(pearson_r) > 0.5, "white", "#1A1A1A"),
    r_label   = sprintf("%.2f", pearson_r)
  )

p <- ggplot(plot_df, aes(comparison1, comparison2)) +
  geom_tile(aes(fill = pearson_r), colour = "white", linewidth = 1.2) +
  geom_text(aes(label = r_label, colour = text_col),
            size = 4.2, family = BASE_FAMILY, fontface = "bold",
            show.legend = FALSE) +
  scale_fill_gradient2(
    low    = C_DEPLETED,
    mid    = "white",
    high   = C_ENRICHED,
    midpoint = 0,
    limits = c(-1, 1),
    name   = "Pearson r"
  ) +
  scale_colour_identity() +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Karşılaştırmalar arası gen etki büyüklüğü uyumu"
  ) +
  coord_equal() +
  theme_thesis() +
  theme(
    axis.text.x      = element_text(angle = 25, hjust = 1),
    axis.text.y      = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_blank(),
    legend.position  = "right",
    legend.key.height = unit(1.8, "lines"),
    legend.key.width  = unit(0.6, "lines")
  )

save_figure(p, "Fig08_cross_comparison", width = 6.2, height = 5.2)
