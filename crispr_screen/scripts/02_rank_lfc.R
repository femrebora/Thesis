#!/usr/bin/env Rscript
# Fig02 -- Ranked log2 fold-change ("waterfall") plots, one panel per comparison.
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(ggrepel)
  library(patchwork)
})

pval_thresh <- 0.05
top_n       <- 6

make_rank <- function(comp) {
  df <- load_gene_summary(comp, pval_thresh) %>%
    arrange(LFC) %>%
    mutate(rank = row_number(),
           direction = classify_direction_concordant(LFC, neg_p, pos_p, pval_thresh),
           hl = direction != "Anlamlı değil")
  n_genes <- nrow(df)
  yr      <- range(df$LFC, na.rm = TRUE)
  lab_df  <- df %>% filter(hl) %>% group_by(direction) %>%
    arrange(desc(abs(LFC))) %>% slice_head(n = top_n) %>% ungroup() %>%
    # Push labels into empty corners: Artan (orange) top-left, Azalan (blue) bottom-right.
    mutate(
      nx = ifelse(direction == "Artan", -n_genes * 0.16,  n_genes * 0.16),
      ny = ifelse(direction == "Artan",  diff(yr) * 0.12, -diff(yr) * 0.12)
    )

  p <- ggplot(df, aes(rank, LFC)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80",
               linewidth = 0.25) +
    geom_point(data = subset(df, !hl), colour = C_NS, size = 1.3, alpha = 0.5) +
    geom_point(data = subset(df, hl), aes(colour = direction),
               size = 2.2, alpha = 0.9)

  # No gene may clear a strict FDR threshold; report the null result on-panel
  # instead of drawing an unlabelled plot that looks like a rendering failure.
  if (nrow(lab_df) > 0) {
    p <- p + geom_text_repel(data = lab_df, aes(label = id, colour = direction),
                    size = 2.9, fontface = "italic", family = BASE_FAMILY,
                    nudge_x = lab_df$nx, nudge_y = lab_df$ny,
                    min.segment.length = 0, max.overlaps = Inf,
                    force = 2, force_pull = 0.1,
                    segment.size = 0.25, segment.colour = "grey60",
                    segment.alpha = 0.7,
                    box.padding = 0.4, point.padding = 0.2,
                    seed = 42, show.legend = FALSE)
  } else {
    p <- p + annotate("text", x = n_genes / 2, y = Inf, vjust = 1.8,
                      label = sprintf("%s < %.2f kriterini karşılayan gen yok",
                                      SIG_SHORT_TR, pval_thresh),
                      size = 2.9, colour = "#777777", family = BASE_FAMILY)
  }

  p +
    scale_colour_manual(values = DIRECTION_COLORS, name = "Yön",
                        breaks = c("Azalan", "Artan"), drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    labs(title = COMPARISONS[comp],
         # Baseline (p-value) figure is thesis-locked and carries no subtitle;
         # the FDR rebuild states its hit count so the null result is explicit.
         subtitle = if (SIG_METRIC == "fdr")
                      sprintf("%d gen (%s < %.2f)", sum(df$hl), SIG_SHORT_TR, pval_thresh)
                    else NULL,
         x = "Gen",
         y = bquote(Log[2]~"kat değişimi")) +
    theme_thesis()
}

panels <- lapply(names(COMPARISONS), make_rank)

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

save_figure(combined, "Fig02_rank_lfc", width = 12, height = 4.8)
