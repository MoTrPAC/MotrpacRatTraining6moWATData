#' @title Create Volcano Plot
#'
#' @description Create a volcano plot with optional labels.
#'
#' @param x a `matrix` or `data.frame` containing differential analysis results.
#' @param logFC character; name of a column in `x` containing log2 fold-changes.
#'   Default is "logFC".
#' @param pval character; name of a column in `x` contain (adjusted) p-values.
#'   Default is "adj.P.Val".
#' @param pval_cutoff numeric; cutoff for p-values to be considered significant.
#'   If provided, a horizontal dashed line will be added to the plot.
#' @param label character; name of a column in `x` used to label the points. Any
#'   points that should not be labelled should be `NA`.
#' @param scale numeric; scaling factor used when saving the plot. Default is
#'   2.5.
#'
#' @md
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom scales extended_breaks
#'
#' @export plot_volcano

plot_volcano <- function(x,
                         logFC = "logFC",
                         pval = "adj.P.Val",
                         pval_cutoff,
                         label,
                         scale = 2.5)
{
  x[["log10_pval"]] <- -log10(x[[pval]])
  x[["nudge_x"]] <- 0.15*sum(range(x[[logFC]], na.rm = T))*sign(x[[logFC]])

  p <- ggplot(data = x,
              mapping = aes(x = !!sym(logFC), y = log10_pval)) +
    geom_point(size = 0.9*scale, alpha = 0.35) +
    scale_x_continuous(
      limits = range_extend(x[[logFC]], nearest = 0.5),
      expand = expansion(mult = 0),
      breaks = scales::extended_breaks(n = 5)
    ) +
    # labs(x = TeX("$log_{2}$(Fold-Change)"),
    #      y = TeX("$-log_{10}$(BH-adjusted p-value)")) +
    theme_bw() +
    theme(text = element_text(size = 6.5*scale, color = "black"),
          line = element_line(size = 0.3*scale, color = "black"),
          panel.border = element_rect(size = 0.4*scale,
                                      fill = NULL,
                                      color = "black"),
          axis.line.y.right = element_blank(),
          axis.text = element_text(size = 5*scale),
          strip.background = element_blank(),
          strip.text = element_text(size = 6.5*scale),
          panel.spacing = unit(0.08*scale, "in"),
          plot.title = element_text(size = 7*scale))

  # y-axis limits
  y_lims <- range_extend(x[["log10_pval"]], nearest = 1)
  y_lims[1] <- 0
  y_lims[2] <- max(3, y_lims[2])

  if (!missing(pval_cutoff)) {
    p <- p +
      geom_hline(yintercept = -log10(pval_cutoff),
                 size = 0.3*scale,
                 lty = "dashed") +
      scale_y_continuous(
        breaks = scales::extended_breaks(n = 5),
        limits = y_lims,
        sec.axis = sec_axis(trans = ~ 10^(-.),
                            breaks = pval_cutoff),
        expand = expansion(mult = c(5e-3, 0.1))
      )
  } else {
    p <- p +
      scale_y_continuous(
        breaks = scales::extended_breaks(n = 5),
        limits = y_lims,
        expand = expansion(mult = c(5e-3, 0.1))
      )
  }

  if (!missing(label)) {
    p <- p +
      geom_label_repel(
        mapping = aes(label = !!sym(label)),
        na.rm = TRUE,
        size = 5*scale*0.352778, # convert points to mm
        max.overlaps = Inf,
        nudge_x = x[["nudge_x"]],
        nudge_y = 0.1,
        fill = alpha("white", 0.65),
        force = 10, seed = 0,
        color = "darkred",
        min.segment.length = 0,
        label.padding = 0.15)
  }

  return(p)
}



