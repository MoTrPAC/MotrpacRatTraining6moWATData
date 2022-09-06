#' @title Create Volcano Plot
#'
#' @description Create a volcano plot.
#'
#' @param x a `matrix` or `data.frame` containing differential analysis results.
#' @param pval_cutoff numeric; cutoff for p-values to be considered significant.
#'   Adds a dashed horizontal line to the plot.
#' @param scale numeric; scaling factor used when saving the plot. Default is
#'   2.5.
#' @param colors character; length 3 vector of significantly-negative,
#'   significantly-positive, and non-significant (NS) logFC colors.
#'
#' @returns A `ggplot2` object.
#'
#' @md
#'
#' @import ggplot2
#' @importFrom scales extended_breaks
#' @importFrom data.table setDT `:=` `.N`
#'
#' @export plot_volcano

plot_volcano <- function(x,
                         pval_cutoff = 0.05,
                         scale = 2.5,
                         colors = c("#3366ff", "darkred", "grey"))
{
  setDT(x)

  p <- ggplot(x, aes(x = logFC, y = log10_pval)) +
    geom_point(aes(color = sign_logFC), alpha = 0.5, size = 3) +
    scale_x_continuous(
      limits = range(x[["logFC"]], na.rm = T)*1.1,
      expand = expansion(mult = 0),
      breaks = scales::extended_breaks(n = 5)
    ) +
    theme_bw() +
    theme(text = element_text(size = 6.5*scale, color = "black"),
          line = element_line(size = 0.3*scale, color = "black"),
          axis.line.y.right = element_blank(),
          axis.text = element_text(size = 5*scale, color = "black"),
          axis.title = element_text(size = 6.5*scale, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 6.5*scale,
                                    hjust = 0, #face = "bold",
                                    margin = margin(0, 0, 7*scale, 0)),
          panel.spacing = unit(0.08*scale, "in"),
          plot.title = element_text(size = 7*scale,
                                    color = "black"),
          # plot.margin = unit(c(5, 5, 5, 5)*scale,
          #                    units = "pt"),
          legend.position = "none",
          plot.tag = element_text(size = 10*2.5,
                                  color = "black",
                                  face = "bold"))

  # y-axis limits
  y_lims <- range_extend(x[["log10_pval"]], nearest = 1)
  y_lims[1] <- 0
  y_lims[2] <- max(3, y_lims[2])


  p <- p +
    geom_hline(yintercept = -log10(pval_cutoff),
               size = 0.3*scale,
               lty = "dashed") +
    scale_y_continuous(
      breaks = scales::extended_breaks(n = 5),
      limits = y_lims,
      sec.axis = sec_axis(trans = ~ 10^(-.),
                          breaks = pval_cutoff),
      expand = expansion(mult = c(0.01, 0.1))
    ) +
    scale_color_manual(values = colors,
                       breaks = levels(x[["sign_logFC"]]))

  # add annotations
  p <- p +
    geom_label(data = unique(x[, list(contrast, label)]),
               aes(label = label, x = -Inf, y = Inf),
               size = 5*2.5*0.35, label.size = NA,
               label.padding = unit(4, "pt"),
               fill = alpha("white", 0.5),
               hjust = 0, vjust = 0) +
    coord_cartesian(clip = "off")

  return(p)

  ## Label points (old code)
  # if (!missing(label)) {
  #   label_args <- list(mapping = aes(label = !!sym(label)),
  #                      na.rm = TRUE,
  #                      size = 5*scale*0.352778, # convert points to mm
  #                      max.overlaps = Inf,
  #                      nudge_x = x[["nudge_x"]],
  #                      nudge_y = 0.1,
  #                      fill = alpha("white", 0.65),
  #                      force = 10, seed = 0,
  #                      color = "darkred",
  #                      min.segment.length = 0,
  #                      label.padding = 0.15)
  #   label_args <- modifyList(x = label_args, val = list(...),
  #                            keep.null = TRUE)
  #
  #   p <- p +
  #     do.call(what = geom_label_repel, args = label_args)
  # }
}

