#!/usr/bin/env Rscript
# Fig07 -- Replicate correlation of sgRNA counts (log10) per condition.
source("utils.R")
suppressPackageStartupMessages({
  library(ggrastr)
  library(patchwork)
})

count_file <- file.path(DATA_DIR, "count_table.txt")
raw <- read_tsv(count_file, show_col_types = FALSE)

pairs <- list(
  c("S0_1", "S0_2"),
  c("Spost_1", "Spost_2"),
  c("R0_1", "R0_2"),
  c("Rpost_1", "Rpost_2")
)

# Shared, equal axis limits across all panels for fair visual comparison.
L <- ceiling(max(vapply(unlist(pairs), function(c) max(log10(raw[[c]] + 1), na.rm = TRUE),
                        numeric(1))))

make_corr <- function(pr) {
  x <- log10(raw[[pr[1]]] + 1)
  y <- log10(raw[[pr[2]]] + 1)
  mask <- x > 0 | y > 0
  rho <- cor(x[mask], y[mask], method = "spearman")
  cond <- sub("_.*", "", pr[1])
  df <- data.frame(x = x, y = y)

  ggplot(df, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey60", linewidth = 0.3) +
    geom_point_rast(size = 0.5, alpha = 0.25, colour = COND_COLORS[cond],
                    raster.dpi = 300) +
    geom_rug(data = subset(df, mask), colour = COND_COLORS[cond],
             alpha = 0.08, length = unit(0.02, "npc")) +
    geom_smooth(method = "lm", se = FALSE, colour = "#1A1A1A", linewidth = 0.6) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.12, vjust = 1.4,
             label = sprintf("rho == %.2f", rho), parse = TRUE,
             size = 3.3, colour = "#1A1A1A", family = BASE_FAMILY) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 3.4,
             label = sprintf("n = %s saptanan", format(sum(mask), big.mark = ".")),
             size = 2.8, colour = "#555555", family = BASE_FAMILY) +
    labs(title = paste0(COND_LABELS[cond], " (tekrar 1 / 2)"),
         x = bquote(Log[10]~"(okuma + 1), tekrar 1"),
         y = bquote(Log[10]~"(okuma + 1), tekrar 2")) +
    coord_equal(xlim = c(0, L), ylim = c(0, L)) +
    theme_thesis()
}

panels <- lapply(pairs, make_corr)

combined <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title    = "Tekrar uyumu",
    subtitle = "Spearman ρ; en az bir tekrarda saptanan sgRNA'lar üzerinden, kesikli çizgi x = y",
    tag_levels = "A"
  ) &
  tag_theme() &
  theme(plot.title    = element_text(margin = margin(l = 16, b = 2)),
        plot.subtitle = element_text(margin = margin(l = 16, b = 6)))

save_figure(combined, "Fig07_replicate_corr", width = 8.5, height = 8.5)
