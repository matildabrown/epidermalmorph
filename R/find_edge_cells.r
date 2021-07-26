#' Find edges of image, identify which get.edges.from are on the edge and which are
#'
#' @param get.edges.from SpatialPolygonsDataFrame. This object will be used to define the edges.
#' @param find.edge.cells SpatialPolygonsDataFrame. This is the object that will be evaluated.
#'                        Defaults to \code{get.edges.from} (i.e. finds which cells of
#'                        \code{get.edges.from} lie on the edge of \code{get.edges.from} )
#' @param edge.buffer.width Numeric. The width (in pixels) that defines the edge. Defaults to 2.
#' @param rotated Numeric. The angle that the image has been rotated. Defaults to NULL.
#'
#' @return A list with two components.
#' #' \describe{
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
if(!is.null(rotated)) get.edges.from <- maptools::elide(get.edges.from, rotate=0-rotated)

x1 <- 0
x2 <- edge.buffer.width
x3 <- raster::bbox(get.edges.from)[1,2]-edge.buffer.width
x4 <- raster::bbox(get.edges.from)[1,2]

y1 <- 0
y2 <- edge.buffer.width
y3 <- raster::bbox(get.edges.from)[2,2]-edge.buffer.width
y4 <- raster::bbox(get.edges.from)[2,2]


#create edge polygons
top_edge <- as(raster::extent(x1,x4,y3,y4), 'SpatialPolygons') #define bottom and right edge to remove partial get.edges.from
bottom_edge <- as(raster::extent(x1,x4,y1,y2), 'SpatialPolygons')
left_edge <- as(raster::extent(x1,x2,y1,y4), 'SpatialPolygons')
right_edge <- as(raster::extent(x3,x4,y1,y4), 'SpatialPolygons')
keep_edges <- raster::bind(top_edge, left_edge)
lose_edges <- raster::bind(bottom_edge, right_edge)
edges <- raster::bind(keep_edges, lose_edges)

if(!is.null(rotated)) edges <- maptools::elide(edges, rotate=rotated)

#identify edge get.edges.from (where lose is edge get.edges.from on right or bottom edges)
cells_edge <-which(colSums(rgeos::gIntersects(find.edge.cells, edges, byid=TRUE)==T, na.rm = TRUE) > 0L)
cells_lose <-which(colSums(rgeos::gIntersects(find.edge.cells, lose_edges,byid=TRUE)==T, na.rm = TRUE) > 0L)
return(list(edge=cells_edge, lose=cells_lose))
  }



