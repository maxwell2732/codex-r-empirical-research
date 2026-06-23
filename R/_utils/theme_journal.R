# ------------------------------------------------------------------------------
# File:     R/_utils/theme_journal.R
# Purpose:  Shared figure aesthetics for the project — a small, opinionated
#           palette + ggplot2 theme that reads as Nature / Science / Lancet rather
#           than "default ggplot blue and red". Every figure in production R/
#           and in explorations/ should source this and use the helpers below.
# Usage:
#           source("R/_utils/theme_journal.R")
#           ggplot(...) + theme_journal() + scale_fill_manual(values = pal_journal[...])
# ------------------------------------------------------------------------------

# Cool, watercolor-leaning palette: predominantly blue with mint-teal and
# lavender-lilac accents. Reserves a deep navy for line strokes / accents and
# a muted slate for neutral text. Names are colour-descriptive so future
# palette swaps stay self-documenting.
pal_journal <- c(
  blue    = "#7FA8D9",   # periwinkle blue — primary group / default histogram fill
  teal    = "#82C8C2",   # soft mint-teal  — secondary group / "Yes" vs "No" partner
  lilac   = "#B49CCF",   # lavender lilac  — third category / negative diverging end
  navy    = "#3B6EA5",   # deep navy       — line strokes, accents
  slate   = "#6B85A0"    # muted slate     — fifth, neutral text accent
)

#' Journal-style ggplot2 theme.
#'
#' theme_classic-derived: serif font, bold axis titles, no gridlines,
#' panel-only axis lines at linewidth 0.4. Pair with `pal_journal` for
#' fills/colours.
#'
#' @param base_size Base font size in pt (default 11; bump to 12-14 for posters).
theme_journal <- function(base_size = 11) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required. Run scripts/setup_r.R.")
  }
  ggplot2::theme_classic(base_size = base_size, base_family = "serif") +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold",
                                               size = base_size + 1, hjust = 0,
                                               margin = ggplot2::margin(b = 4)),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, hjust = 0,
                                               margin = ggplot2::margin(b = 6),
                                               colour = "grey25"),
      plot.caption     = ggplot2::element_text(size = base_size - 3, hjust = 1,
                                               colour = "grey45"),
      axis.title       = ggplot2::element_text(face = "bold"),
      axis.title.x     = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y     = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      axis.text        = ggplot2::element_text(colour = "black"),
      axis.line        = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = ggplot2::element_line(colour = "black", linewidth = 0.4),
      panel.grid       = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key       = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(face = "bold"),
      plot.margin      = ggplot2::margin(8, 8, 8, 8)
    )
}
