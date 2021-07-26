#' Find the shortest paths between each cell
#'
#' @param cells SpatialPolygonsDataFrame.
#'
#' @return A matrix with the shortest paths between each cell
#' @import methods spatialreg
#'


cell_graph_shortest_paths <- function(cells){


boundaries <- cells
nb_q <- spdep::poly2nb(boundaries, snap=5)

nb_B <- spdep::nb2listw(nb_q, style="B", zero.policy=TRUE)
B <- as(nb_B, "symmetricMatrix")
g1 <- igraph::graph.adjacency(B, mode="undirected")
sp_mat <- igraph::shortest.paths(g1)
return(sp_mat)
}


