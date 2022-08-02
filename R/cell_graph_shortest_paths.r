#' Find the shortest paths between each cell
#'
#' @param cells \code{sf} object.
#' @param snap.tolerance Numeric. Value to buffer each cell to identify neighbours.
#'
#' @return A matrix with the shortest paths between each cell
#' @import methods sf
#'
#' @export


cell_graph_shortest_paths <- function(cells, snap.tolerance){

  boundaries <- sf::st_buffer(cells, snap.tolerance)
  x <- sf::st_intersects(boundaries)
  m = data.frame(row.id=rep(seq_along(x), lengths(x)), col.id=unlist(x))
  g <- igraph::graph_from_data_frame(m, directed = F)
  sp_mat <- igraph::shortest.paths(g)
  return(sp_mat)
}



