#!/usr/bin/env Rscript
# Fig04 -- Top significant genes by |log2 FC| as horizontal bars, per comparison.
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages(library(patchwork))

pval_thresh <- 0.05
top_n       <- 12

make_bar <- function(comp) {
  df <- load_gene_summary(comp, pval_thresh) %>%
    mutate(direction = classify_direction_concordant(LFC, neg_p, pos_p, pval_thresh)) %>%
    filter(direction != "Anlamlı değil") %>%
    mutate(abs_lfc = abs(LFC)) %>%
    arrange(desc(abs_lfc)) %>%
    slice_head(n = top_n) %>%
    # Order: direction first (Artan group on top, Azalan below), then |LFC|
    # descending within each group so the strongest effects sit at the outside.
    arrange(direction, abs_lfc) %>%
    mutate(gene = factor(id, levels = id))

  ggplot(df, aes(LFC, gene, fill = direction)) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey75",
               linewidth = 0.25) +
    geom_col(width = 0.7) +
    facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = DIRECTION_COLORS, name = "Yön",
                      breaks = c("Azalan", "Artan")) +
    labs(title = COMPARISONS[comp],
         x = bquote(Log[2]~"kat değişimi"), y = NULL) +
    theme_thesis() +
    theme(axis.text.y = element_text(face = "italic"),
          panel.grid.major.y = element_blank())
}

panels <- lapply(names(COMPARISONS), make_bar)

combined <- wrap_plots(panels, nrow = 1, guides = "collect") +
  plot_annotation(
    tag_levels = "A"
  ) &
  tag_theme() &
  theme(
    legend.position = "bottom",
    plot.title    = element_text(margin = margin(l = 16, b = 2)),
    plot.subtitle = element_text(margin = margin(l = 16, b = 6))
  )

save_figure(combined, "Fig04_barplot", width = 12, height = 5.2)
