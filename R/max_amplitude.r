#' Find the maximum amplitude of cell sinuosity
#'
#' @param cell SpatialPolygon. The cell.
#' @param junction_points Matrix. The x and y coordinates defining the simplified
#' cell, as returned by \code{cell.simplify()}
#'
#' @return Numeric. The maximum distance from the simplified cell edge to the
#' actual cell wall.
#'
#' @import methods
#' @export


max_amplitude <- function(cell, junction_points){

  actuallines = as(cell, 'SpatialLines')

#how many sides does our simplified shape have
n_gon = nrow(junction_points)-1 #number of junction points

# FIND MAX AMPLITUDE OF WIBBLES IE MAX DISTANCE BETWEEN ACTUAL SHAPE AND SIMPLIFIED SHAPE
bestdist=0

for (i in 1:n_gon){

  A = junction_points[i,]
  if (i != n_gon) {
    B = junction_points[i+1,]
  } else {
    B = junction_points[1,]
  }
  slopeAB <-( B[2]- A[2])/(B[1]-A[1])  #rise/run
  perpslope <- -1/slopeAB

  lineAB = raster::spLines(rbind(A,B))


  if(dist(rbind(A,B))/2>1){
    t = suppressWarnings( sp::spsample(lineAB, n = dist(rbind(A,B))/2, type='regular') )

    for (j in 1:length(t)){
      t.j = t@coords[j,]
      length = 1500
      delta.x = (length^0.5*(perpslope^2 +1)^0.5)/(1+perpslope^2)
      delta.y = (length^0.5*perpslope*(perpslope^2 +1)^0.5)/(1+perpslope^2)
      cdistminus = dist(rbind(raster::coordinates(cell), c(t.j[1]-delta.x, t.j[2]-delta.y)))
      cdistplus = dist(rbind(raster::coordinates(cell), c(t.j[1]+delta.x, t.j[2]+delta.y)))

      if(prevR::point.in.SpatialPolygons(t.j[1], t.j[2], cell)==T){
        if (cdistplus>cdistminus){
          lineTJ = raster::spLines(rbind(c(t.j[1], t.j[2]), c(t.j[1]+delta.x, t.j[2]+delta.y)))
        } else {lineTJ = raster::spLines(rbind(c(t.j[1]-delta.x, t.j[2]-delta.y), c(t.j[1], t.j[2])))}
      }else{
        if (cdistminus>cdistplus){
          lineTJ = raster::spLines(rbind(c(t.j[1], t.j[2]), c(t.j[1]+delta.x, t.j[2]+delta.y)))
        } else {
          lineTJ = raster::spLines(rbind(c(t.j[1]-delta.x, t.j[2]-delta.y), c(t.j[1], t.j[2])))
        }
      }

      if(length(raster::intersect(lineTJ, actuallines))>0){
      d = raster::coordinates(raster::intersect(lineTJ, actuallines))
      d2 <-  dist(rbind(t.j, d))
      if(nrow(d)>1){d=d[which.min(d2[1:nrow(d)]),]}
      d2 <-  dist(rbind(t.j, d))
      if(length(d!=0)){
        if(d2 > bestdist) {
          bestdist <- d2
          bestd <- d
          bestt <- t.j
        }
      }
      }
    }
  }
}
return(as.numeric(bestdist))
}
