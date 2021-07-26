#' Find the localMaxima in a series
#'
#' @param x Numeric
#'
#' @return A numeric vector with the positions of localMaxima
#' @details This function comes from user Tommy 's answer on SO:
#' https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima


localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}
