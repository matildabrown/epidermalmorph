#' Find edges of image, identify which get.edges.from are on the edge and which are
#'
#' @param get.edges.from \code{sf} object. This object will be used to define the edges.
#' @param find.edge.cells \code{sf}object. This is the object that will be evaluated.
#'                        Defaults to \code{get.edges.from} (i.e. finds which cells of
#'                        \code{get.edges.from} lie on the edge of \code{get.edges.from} )
#' @param edge.buffer.width Numeric. The width (in pixels) that defines the edge. Defaults to 2.
#' @param rotated Numeric. The angle (in degrees) that the image has been rotated. Defaults to NULL.
#'
#' @return A list with two components.
#'  \describe{
#' \item{edge}{Numeric. Vector containing the indices of all edge get.edges.from.}
#' \item{lose}{Numeric. Vector containing the indices of all get.edges.from on the
#'             bottom or right-hand edge.}
#' }
#'
#' @details The \code{get.edges.from} parameter can be the output from
#' \code{image_to_poly()}.
#'
#' @import methods
#'
#' @export


find_edge_cells <- function(get.edges.from, find.edge.cells=NULL, edge.buffer.width=2, rotated=NULL){

if(is.null(find.edge.cells)) find.edge.cells <- get.edges.from
if(!is.null(rotated)){
center <- sf::st_centroid(sf::st_union(find.edge.cells))
 sf::st_geometry(get.edges.from) <- tran(sf::st_geometry(get.edges.from), 0-rotated, center)
}

x1 <- 0
x2 <- edge.buffer.width
x3 <- sf::st_bbox(get.edges.from)[3]-edge.buffer.width
x4 <- sf::st_bbox(get.edges.from)[3]

y1 <- 0
y2 <- edge.buffer.width
y3 <- sf::st_bbox(get.edges.from)[4]-edge.buffer.width
y4 <- sf::st_bbox(get.edges.from)[4]


#create edge polygons
top_edge <- sf::st_polygon(list(matrix(c(x1,x1,x4,x4,x1,y3,y4,y4,y3,y3), ncol=2)))#define bottom and right edge to remove partial get.edges.from
bottom_edge <- sf::st_polygon(list(matrix(c(x1,x1,x4,x4,x1,y1,y2,y2,y1,y1), ncol=2)))
left_edge <- sf::st_polygon(list(matrix(c(x1,x1,x2,x2,x1,y1,y4,y4,y1,y1), ncol=2)))
right_edge <- sf::st_polygon(list(matrix(c(x3,x3,x4,x4,x3,y1,y4,y4,y1,y1), ncol=2)))

keep_edges <- sf::st_union(top_edge, left_edge)
lose_edges <- sf::st_union(bottom_edge, right_edge)
edges <- sf::st_union(keep_edges, lose_edges)

if(!is.null(rotated)) edges <- edges * rot(pracma::deg2rad(rotated))

#identify edge get.edges.from (where lose is edge get.edges.from on right or bottom edges)
cells_edge <-which(lengths(sf::st_intersects(find.edge.cells,edges))>0)
cells_lose <-which(lengths(sf::st_intersects(find.edge.cells,lose_edges))>0)
return(list(edge=cells_edge, lose=cells_lose))
  }



