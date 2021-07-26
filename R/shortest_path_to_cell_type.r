#' Find the length of the shortest path to any cell of a particular type
#'
#' @param cells SpatialPolygonsDataFrame.
#' @param path.matrix Matrix with path lengths between each pair of cells.
#' Recommended to use output from \code{cell_graph_shortest_paths()}.
#' @param cell.type.value Numeric. The value of stomata.
#'
#' @return A numeric vector with the distance (in cells) to the nearest cell
#' with value \code{cell.type.value}.
#'
#' @details For cells with the same value as \code{cell.type.value}, this
#'  algorithm finds the path length to the nearest other cell of that type.
#'
#' @export


shortest_path_to_cell_type <- function(cells, path.matrix, cell.type.value){

  # Values
  referenceCol <- cells[[1]]

  # Rename spatial matrix
  path.matrix <- as.data.frame(path.matrix)
  colnames(path.matrix) <- paste0(referenceCol)
  path.matrix$id <- rownames(cells@data)

  diag(path.matrix) <- 10000
  dist_from_stom <- Rfast::rowMins(as.matrix(path.matrix[,which(referenceCol==cell.type.value)]), value=T)

  return(dist_from_stom)

}
