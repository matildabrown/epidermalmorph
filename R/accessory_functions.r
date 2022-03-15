#Accessory functions

`%notin%` <- Negate(`%in%`)

rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
tran = function(geo, ang, center) {
  (geo - center) * rot(ang * pi / 180) + center
}

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


