#' Find the mean distance to k nearest neighbours
#'
#' @param cells \code{sf} object.
#' @param k Numeric. The k nearest neighbours to calculate distances to.
#'
#' @return
#' If \code{k} is a single value, returns a vector of same length as
#' \code{cells}, containing the distance (in pixels) between each cell and its
#'kth nearest neighbour. If \code{k} is a vector, returns a matrix with
#'\code{length(k)} variables, with the distance (in pixels) from each
#' cell to its nearest k neighbours.
#' @details This code can be very slow if \code{cells} contains many polygons
#' or \code{k} has more than two values.
#' @export


cell_cell_distances <- function(cells, k){

cellsdist <- sf::st_distance(cells, byid=T)
k.length=length(k)

for (j in 1:nrow(cells)){
  cellsdist[j,j] <- 10^10
}

if (k.length==1) {
  mindist <- Rfast::rowMins(cellsdist, value=T)
}

if (k.length>1) {
mindist=matrix(nrow=nrow(cells), ncol=k.length)
for (i in 1:length(k)){
  k.i <- k[i]

mindist.j <- NULL
for (j in 1:nrow(cellsdist)){
  mindist.j <- c(mindist.j, Rfast::nth(cellsdist[,j], k.i, descending = F))
}
mindist[,i] <- mindist.j
}

}

return(mindist)

}
