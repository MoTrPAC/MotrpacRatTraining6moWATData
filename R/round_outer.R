#' @title Round to Outer Limits
#'
#' @description Round limits out to nearest value.
#'
#' @param x a numeric vector.
#' @param nearest numeric; the range of \code{x} will be rounded out to the
#'   value specified by \code{nearest}. Default is 1.
#'
#' @returns A numeric vector of length 2 specifying the range of values after
#'   rounding outward to the nearest value provided by `nearest`.
#'
#' @export round_outer
#'
#' @author Tyler Sagendorf
#'
#' @examples
#' set.seed(0)
#' x <- runif(5, min = -10, max = 10)
#' x
#'
#' round_outer(x) # nearest whole number
#' round_outer(x, nearest = 0.1) # nearest tenth

round_outer <- function(x, nearest = 1) {
  x <- x/nearest
  xmin <- min(x, na.rm = T)
  xmax <- max(x, na.rm = T)
  c(floor(xmin), ceiling(xmax))*nearest
}

