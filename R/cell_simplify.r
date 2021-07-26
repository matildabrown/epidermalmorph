#' Get the simplified cell shape from junction points
#'
#' @param cell SpatialPolygon of the original cell
#' @param cell.junctions sp::SpatialPoints object containing all junction points
#' between cells. Can be obtained from the Analyse Skeleton function in ImageJ.
#' @return An ordered matrix containing the x and y coordinates defining the simplified
#' cell.
#' @details The returned matrix has n+1 rows (the first and last rows are the
#' same, so points can be joined up to form a closed polygon).
#' @export


cell.simplify <- function(cell, cell.junctions){

  #get centroid
  centroid <- rgeos::gCentroid(cell) #


  # EXTRACT THE JUNCTION POINTS FROM THE raster::buffer ZONE, MAP ONTO CELL BOUNDARY
  # AND ORDER THEM SO A SENSIBLE SIMPLIFIED POLYGON CAN BE FOUND

  #check that raster::buffer value includes all junctions
  b <- raster::buffer(cell, width=8)

  # get the junctions that fall within the boundary
  if("SpatialPoints" %notin% class(cell.junctions))  cell.junctions <- SpatialPoints(cell.junctions)

  c.j <- sp::over(cell.junctions,b)

  c.j <- cell.junctions[which(is.na(c.j)==F)]

  #this is where we do an angle test to order the points (anticlockwise?)
  lines=NULL;angles=NULL;q=NULL
  for( i in 1:nrow(c.j@coords)){
    # plot(raster::spLines(rbind(centroid@coords, c.j@coords[i,])), add=T)
    lines <- c(lines, raster::spLines(rbind(centroid@coords, c.j@coords[i,])))
    slope <- (c.j@coords[i,2]-centroid@coords[,2])/(c.j@coords[i,1]-centroid@coords[,1])
    angles=c(angles,slope)
    #find which quadrat
    if(c.j@coords[i,2]>centroid@coords[,2]){ #if point below centroid
      if(c.j@coords[i,1]>centroid@coords[,1]) { #if point to right of centroid
        quad <- 1
      } else { #point to left of centroid
        quad <- 2
      }
    } else { #point above centroid
      if(c.j@coords[i,1]>centroid@coords[,1]) { #if point to right of centroid
        quad <- 4
      } else {
        quad <- 3
      }
    }
    q=c(q,quad)
  }
  c.j.ordered=sp::SpatialPoints(c.j@coords[order(q,angles),1:2]) #order points by quadrat and slope (to/from centroid)

  c.j.fixed.list <- list()
  #move the points onto the edge of the polygon (i.e. out of the raster::buffer zone)
  for (i in 1:nrow(c.j.ordered@coords)){
    c.j.fixed.list[[i]] <- rgeos::gNearestPoints(cell, c.j.ordered[i])[1]@coords# finds the nearest point on polygon to wall junction
  }

  c.j.fixed <- do.call("rbind", c.j.fixed.list)
  c.j.fixed <- rbind(c.j.fixed, c.j.fixed[1,]) #add first point to close the line
  rownames(c.j.fixed) <- 1:nrow(c.j.fixed)


  junction_points <- sp::SpatialPoints(c.j.fixed) #1,2,3,...n,1 in order
  junction_points <- junction_points@coords #1,2,3,...n,1 in order

  return(junction_points)
}
