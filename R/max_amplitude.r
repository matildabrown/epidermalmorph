#' Find the maximum amplitude of cell sinuosity
#'
#' @param cell \code{sf} object; the cell.
#' @param junction_points Matrix. The x and y coordinates defining the simplified
#' cell, as returned by \code{cell.simplify()}
#'
#' @return Numeric. The maximum distance from the simplified cell edge to the
#' actual cell wall.
#'
#' @import methods
#' @export


max_amplitude <- function(cell, junction_points){

  actuallines = suppressWarnings(sf::st_cast(cell, 'LINESTRING'))

#how many sides does our simplified shape have
n_gon = nrow(junction_points)-1 #number of junction points
coords.cent <- suppressWarnings(sf::st_coordinates(sf::st_centroid(cell)))

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

  lineAB = sf::st_linestring(rbind(A,B))


  if(dist(rbind(A,B))/2>1){
    t = suppressWarnings( sf::st_sample(lineAB, size = dist(rbind(A,B))/2, type='regular') )[[1]]

    for (j in 1:nrow(t)){
      t.j = t(t[j,])


      length = 1500
      delta.x = (length^0.5*(perpslope^2 +1)^0.5)/(1+perpslope^2)
      delta.y = (length^0.5*perpslope*(perpslope^2 +1)^0.5)/(1+perpslope^2)
      cdistminus = dist(rbind(coords.cent, c(t.j[1]-delta.x, t.j[2]-delta.y)))
      cdistplus = dist(rbind(coords.cent, c(t.j[1]+delta.x, t.j[2]+delta.y)))
      in.cell <- st_intersects(st_point(t.j),cell)[[1]]
      if(length(in.cell)==1){
        if (cdistplus>cdistminus){
          lineTJ = sf::st_linestring(rbind(c(t.j[1], t.j[2]), c(t.j[1]+delta.x, t.j[2]+delta.y)))
        } else {lineTJ = sf::st_linestring(rbind(c(t.j[1]-delta.x, t.j[2]-delta.y), c(t.j[1], t.j[2])))}
      }else{
        if (cdistminus>cdistplus){
          lineTJ = sf::st_linestring(rbind(c(t.j[1], t.j[2]), c(t.j[1]+delta.x, t.j[2]+delta.y)))
        } else {
          lineTJ = sf::st_linestring(rbind(c(t.j[1]-delta.x, t.j[2]-delta.y), c(t.j[1], t.j[2])))
        }
      }
      t.j <- as.data.frame(t.j)
      colnames(t.j) <- c("X","Y")
      if(length(sf::st_intersects(lineTJ, actuallines)[[1]])>0){
       d <-  sf::st_coordinates(sf::st_intersection(lineTJ,actuallines))
       d <- as.data.frame(d) #why does this work on a separate line
       d <- d[,1:2] #get rid of the extra column

       colnames(d) <- colnames(t.j)
       d2 <-  dist(dplyr::bind_rows(t.j, d))
      if(length(d2)>1){
        d=d[which.min(d2[1:nrow(d)]),]

        d2 <-  dist(dplyr::bind_rows(t.j, d))
        }

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
