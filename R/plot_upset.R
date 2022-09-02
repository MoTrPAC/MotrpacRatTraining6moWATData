#' @title Create Annotated UpSet Plot
#'
#' @description Create an UpSet plot annotated with intersection and set sizes.
#'   Essentially a modified version of \code{\link[ComplexHeatmap]{UpSet}}.
#'
#' @inheritParams ComplexHeatmap::make_comb_mat
#' @inheritParams enrichmat
#' @param top_n_comb integer; number of largest intersections to display.
#'   Default is 10.
#' @param scale_set_bars logical;
#' @param row_labels character; labels for each set (row of UpSet plot).
#'   Defaults to rownames of \code{\link[ComplexHeatmap]{make_comb_mat}} output.
#' @param row_names_gp graphical parameters specified by
#'   \code{\link[grid]{gpar}}. Modifies `row_labels`.
#' @param annotation_name_gp graphical parameters specified by
#'   \code{\link[grid]{gpar}}. Passed to
#'   \code{\link[ComplexHeatmap]{HeatmapAnnotation}}. Modifies the barplot
#'   titles.
#' @param cell_height `unit` object; the height of each heatmap cell. Default is
#'   `unit(4, "mm")`.
#' @param cell_width same as `cell_height`.
#' @param extend numeric;
#' @param rot numeric; angle in \[0,360\] to rotate set intersection labels.
#' @param hjust numeric; a numeric vector specifying horizontal justification of
#'   intersection size bar labels.
#' @param vjust numeric; A numeric vector specifying vertical justification of
#'   intersection size bar labels.
#' @param bar_label_gp graphical parameters specified by
#'   \code{\link[grid]{gpar}}. Passed to
#'   \code{\link[ComplexHeatmap]{HeatmapAnnotation}}. Modifies the annotation
#'   labels.
#'
#' @md
#'
#' @details See \code{\link[ComplexHeatmap]{UpSet}} for details and sections
#'   mentioned in Arguments.
#'
#' @import ComplexHeatmap
#' @importFrom grid gpar unit grid.rect grid.lines grid.text grid.points
#' @importFrom utils modifyList
#' @importFrom grDevices dev.off
#'
#' @export plot_upset
#'
#' @references Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister,
#'   H. (2014). UpSet: Visualization of Intersecting Sets. *IEEE transactions on
#'   visualization and computer graphics, 20*(12), 1983--1992.
#'   \url{https://doi.org/10.1109/TVCG.2014.2346248}

plot_upset <- function(...,
                       mode = c("distinct", "intersect", "union"),
                       top_n_comb = 10,
                       scale_set_bars = FALSE,
                       row_labels = character(0),
                       row_names_gp = gpar(fontsize = 6.5),
                       annotation_name_gp = gpar(fontsize = 6.5,
                                                 fontface = "bold"),
                       heatmap_args = list(),
                       cell_height = unit(4, "mm"),
                       cell_width = cell_height,
                       extend = 0.2,
                       rot = 45, hjust = 0, vjust = 0,
                       bar_label_gp = gpar(fontsize = 6),
                       filename = character(0),
                       height = 2.5,
                       width = 2.7,
                       save_args = list())
{
  m <- make_comb_mat(..., mode = mode)
  cs <- comb_size(m)
  ss <- set_size(m)
  m <- m[, order(-cs)[1:min(top_n_comb, length(cs))]]
  cs <- comb_size(m)
  bar_ratio <- ifelse(scale_set_bars, max(ss)/max(cs), 1)*0.8
  row_order <- seq_along(ss)
  column_order <- order(-cs)
  if (identical(row_labels, character(0))) {
    row_labels <- rownames(m)
  }
  # I don't remember why I couldn't just use ComplexHeatmap::UpSet,
  # but I'm sure there was some justification...
  heatmap_args <- modifyList(x = list(
    matrix = m,
    row_labels = row_labels,
    row_names_gp = row_names_gp,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_max_width = max_text_width(row_labels),
    show_heatmap_legend = FALSE,
    row_order = row_order,
    column_order = column_order,
    height = length(ss)*cell_height,
    width = length(cs)*cell_width,
    rect_gp = gpar(type = "none"),
    ## layer function from ComplexHeatmap::UpSet
    layer_fun = function(j, i, x, y, w, h, fill) {
      n_row <- round(1/as.numeric(h[1]))
      n_col <- round(1/as.numeric(w[1]))
      subm <- matrix(pindex(m, i, j), nrow = n_row, byrow = FALSE)
      for (k in seq_len(n_row)) {
        if (k%%2) {
          grid.rect(y = k/n_row, height = 1/n_row, just = "top",
                    gp = gpar(fill = "#f0f0f0", col = NA))
        }
      }
      grid.points(x, y, size = unit(min(3, min(as.numeric(cell_height),
                                               as.numeric(cell_width))),
                                    "mm"), pch = 16,
                  gp = gpar(col = ifelse(pindex(m, i, j),
                                         "black", "#cccccc")))
      for (k in seq_len(n_col)) {
        i_range <- range(which(subm[, k] == 1))
        grid.lines(c(k - 0.5, k - 0.5)/n_col,
                   (n_row - i_range + 0.5)/n_row,
                   gp = gpar(col = "black", lwd = 2))
      }
    },
    ##
    top_annotation = HeatmapAnnotation(
      "Intersection Size" = anno_barplot(
        x = cs,
        extend = extend,
        border = FALSE,
        gp = gpar(fill = "black"),
        height = unit(1.2, "in"),
        which = "column",
        axis = FALSE
      ),
      annotation_name_gp = annotation_name_gp,
      annotation_name_side = "left",
      annotation_name_rot = 90),
    ##
    right_annotation = HeatmapAnnotation(
      "N differential features" = anno_barplot(
        x = ss,
        extend = extend,
        border = FALSE,
        gp = gpar(fill = "black"),
        width = unit(1.2, "in")*bar_ratio,
        which = "row",
        axis = FALSE
      ),
      which = "row",
      annotation_name_gp = annotation_name_gp,
      annotation_name_side = "bottom",
      annotation_name_rot = 0)
  ),
  val = heatmap_args, keep.null = TRUE) # update arguments

  # Save plot
  if (!identical(filename, character(0))) {
    save_heatmap(filename = filename, height = height,
                 width = width, save_args = save_args)
    on.exit(expr = dev.off())
  }
  # draw base UpSet plot
  ht <- do.call(what = Heatmap, args = heatmap_args)
  draw(ht)
  column_order <- heatmap_args$column_order
  row_order <- rev(heatmap_args$row_order)
  # Add counts to bars
  decorate_annotation(annotation = "Intersection Size", {
    grid.text(label = cs[column_order],
              x = seq_along(cs),
              y = unit(cs[column_order], "native") + unit(3, "pt"),
              vjust = vjust, hjust = hjust,
              rot = rot, default.units = "native", gp = bar_label_gp)
  })
  decorate_annotation(annotation = "N differential features", {
    grid.text(label = ss[row_order],
              x = unit(ss[row_order], "native") + unit(3, "pt"),
              y = seq(1/(length(ss)*2), 1 - 1/(length(ss)*2),
                      length.out = length(ss)),
              hjust = 0, gp = bar_label_gp)
  })
}


