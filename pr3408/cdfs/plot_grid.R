# Bundle: heat map of the 266-cell grid CSV (needs ggplot2).
# Supports the figure of the #3408 grid claim.
# Heat map of log10 relative error of d/dalpha gamma_lccdf(y | alpha, 1),
# in the style of the plot in the stan-dev/math#3408 comment by mjcarter95.
# Usage: Rscript plot_grid.R IN_CSV OUT_PNG TITLE
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(args[1], stringsAsFactors = FALSE)
title <- args[3]

worst <- max(d$log10_rel_err, na.rm = TRUE)
n_bad <- sum(is.na(d$log10_rel_err))
d$alpha_f <- factor(d$alpha, levels = sort(unique(d$alpha)))
d$ratio_f <- factor(d$ratio, levels = sort(unique(d$ratio)))
d$shown <- pmin(pmax(d$log10_rel_err, -8), 0)

p <- ggplot(d, aes(x = alpha_f, y = ratio_f, fill = shown)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlGn", limits = c(-8, 0),
                       na.value = "black",
                       name = "log10(relative error\nvs high-precision\nreference)") +
  labs(x = "alpha (shape)", y = "y / alpha  (ratio to the mean, beta=1)",
       title = title,
       subtitle = sprintf("worst relative error %.1e; non-finite cells: %d",
                          10^worst, n_bad)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(args[2], p, width = 7.5, height = 6.5, dpi = 150, bg = "white")
cat(sprintf("worst %.2e, non-finite %d\n", 10^worst, n_bad))
