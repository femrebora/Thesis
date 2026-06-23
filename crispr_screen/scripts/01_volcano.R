#!/usr/bin/env Rscript
# Fig01 -- Volcano plots (one panel per comparison), thesis-publication style.
# Custom ggplot2 + ggrepel for full control over colour, labels and layout.
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(ggrepel)
  library(patchwork)
})

pval_thresh <- 0.05
top_n       <- 8

make_volcano <- function(comp) {
  df <- load_gene_summary(comp, pval_thresh)
  # Direction-concordant colouring: a gene is highlighted only when its
  # significant tail matches its LFC sign. Discordant near-threshold points
  # (e.g. TECR/TKTL1/APOC1 in panel B) fall to "Anlamlı değil" -> grey, no label.
  df$direction <- classify_direction_concordant(df$LFC, df$neg_p, df$pos_p, pval_thresh)
  df$hl <- df$direction != "Anlamlı değil"
  n_hl <- sum(df$hl)

  xr <- range(df$LFC, na.rm = TRUE)
  lab_df <- df %>%
    filter(hl) %>%
    arrange(pmin(neg_p, pos_p)) %>%
    slice_head(n = top_n) %>%
    # Target x: Azalan (blue) labels to the left edge, Artan (orange) to the right.
    mutate(target_x = ifelse(LFC < 0, xr[1] - 0.5, xr[2] + 0.5))

  ggplot(df, aes(LFC, log10p)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey80", linewidth = 0.25) +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed",
               colour = "grey80", linewidth = 0.25) +
    # Non-highlighted first (background), then coloured hits on top
    geom_point(data = subset(df, !hl), colour = C_NS,
               size = 1.3, alpha = 0.5) +
    geom_point(data = subset(df, hl), aes(colour = direction),
               size = 2.2, alpha = 0.9) +
    geom_text_repel(data = lab_df, aes(label = id, colour = direction),
                    size = 2.9, fontface = "italic", family = BASE_FAMILY,
                    nudge_x = lab_df$target_x - lab_df$LFC,
                    direction = "y",
                    hjust = ifelse(lab_df$LFC < 0, 1, 0),
                    min.segment.length = 0, max.overlaps = Inf,
                    force = 1, force_pull = 0.15,
                    segment.size = 0.25, segment.colour = "grey60",
                    segment.alpha = 0.7,
                    box.padding = 0.3, point.padding = 0.2,
                    seed = 42, show.legend = FALSE) +
    scale_colour_manual(values = DIRECTION_COLORS, name = "Yön",
                        drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.22)) +
    labs(
      title    = COMPARISONS[comp],
      subtitle = sprintf("%d gen (p < %.2f)", n_hl, pval_thresh),
      x = bquote(Log[2]~"kat değişimi"),
      y = bquote(-Log[10]~"p-değeri")
    ) +
    theme_thesis()
}

panels <- lapply(names(COMPARISONS), make_volcano)

combined <- wrap_plots(panels, nrow = 1, guides = "collect") +
  plot_annotation(
    title   = "Diferansiyel sgRNA temsili",
    tag_levels = "A"
  ) &
  tag_theme() &
  theme(
    legend.position = "bottom",
    plot.title    = element_text(margin = margin(l = 16, b = 2)),
    plot.subtitle = element_text(margin = margin(l = 16, b = 6))
  )

save_figure(combined, "Fig01_volcano", width = 12, height = 4.8)
