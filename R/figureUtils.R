# Shared settings for the figures of the paper.
# save_fig() writes each figure twice, as figures/Fig<n>.png for the manuscript
# and as figures/Fig<n>.tif for PLOS.
# PLOS asks for TIFF with LZW compression, RGB without alpha channel, 300 to 600 dpi,
# a width of 2.63 to 7.5 in, a height of at most 8.75 in, and Arial, Times or Symbol at 8 to 12 pt.

library(ggplot2)

fig_family <- "Arial"
if (!grepl("arial", systemfonts::font_info(fig_family)$family, ignore.case = TRUE)) {
  warning("Arial is not installed, the figures use ", systemfonts::font_info(fig_family)$family, " instead")
}

# base_size 10 gives 8 pt for axis and legend text and 10 pt for titles
theme_fig <- theme_bw(base_size = 10, base_family = fig_family)

# node colors of the networks, separable under color blindness
fig_colors <- c(blue = "#2a78d6", orange = "#eb6834", green = "#1baf7a")

# plots is a ggplot or a list of ggplots, and ... is passed on to gridExtra::arrangeGrob()
# tags puts the panel labels in a row below a single row of panels
# the panels are arranged while the device is open, so that the text is measured in the right font
save_fig <- function(plots, n, width, height, dpi = 600, tags = NULL, ...) {
  stopifnot(width >= 2.63, width <= 7.5, height <= 8.75, dpi >= 300, dpi <= 600)
  if (inherits(plots, "ggplot")) plots <- list(plots)
  dir.create("figures", showWarnings = FALSE)
  draw <- function() {
    if (is.null(tags)) {
      fig <- gridExtra::arrangeGrob(grobs = plots, ...)
    } else {
      stopifnot(length(tags) == length(plots))
      labels <- lapply(tags, grid::textGrob, gp = grid::gpar(fontfamily = fig_family, fontsize = 10))
      fig <- gridExtra::arrangeGrob(grobs = c(plots, labels), nrow = 2,
                                    heights = grid::unit.c(grid::unit(1, "null"), grid::unit(0.2, "in")), ...)
    }
    grid::grid.draw(fig)
  }
  ragg::agg_png(sprintf("figures/Fig%d.png", n), width = width, height = height,
                units = "in", res = dpi, background = "white")
  draw()
  dev.off()
  ragg::agg_tiff(sprintf("figures/Fig%d.tif", n), width = width, height = height,
                 units = "in", res = dpi, compression = "lzw", background = "white")
  draw()
  dev.off()
  invisible(NULL)
}
