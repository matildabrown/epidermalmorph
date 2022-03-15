#' Get the simplified cell shape from junction points
#'
#' @param cell Obect of class \code{sf}; the original cell
#' @param cell.junctions Dataframe with two columns containing all junction points
#' between cells.
#' @param snap.tolerance Numeric. Value to buffer each cell to identify junction points.
#' @return An ordered matrix containing the x and y coordinates defining the simplified
#' cell.
#' @details The returned matrix has n+1 rows (the first and last rows are the
#' same, so points can be joined up to form a closed polygon).
#' @export


cell.simplify <- function(cell, cell.junctions, snap.tolerance){

  #get centroid
  centroid <- suppressWarnings( sf::st_centroid(cell) )


  # EXTRACT THE JUNCTION POINTS FROM THE raster::buffer ZONE, MAP ONTO CELL BOUNDARY
  # AND ORDER THEM SO A SENSIBLE SIMPLIFIED POLYGON CAN BE FOUND

  #check that raster::buffer value includes all junctions
  b <- sf::st_buffer(cell, dist=snap.tolerance)

  # get the junctions that fall within the boundary
  if("sfc_MULTIPOINT" %notin% class(cell.junctions))  cell.junctions <- sf::st_geometry(sf::st_cast(sf::st_multipoint(cell.junctions),"MULTIPOINT"))

  c.j <- sf::st_intersection(cell.junctions,b)


  #this is where we do an angle test to order the points (anticlockwise?)
  lines=NULL;angles=NULL;q=NULL

  #set up coords
  coords.cent <- sf::st_coordinates(centroid)
  coords.cj <-  sf::st_coordinates(c.j)[,1:2]

  for( i in 1:nrow(c.j[[1]])){
    coords.i <- coords.cj[i,]

    lines <- c(lines, sf::st_linestring(rbind(coords.cent,
                                              coords.i)))
    slope <- (coords.i[2]-coords.cent[,2])/(coords.i[1]-coords.cent[,1])
    angles=c(angles,slope)
    #find which quadrat
    if(coords.cj[2]>coords.cent[,2]){ #if point below centroid
      if(coords.i[1]>coords.cent[,1]) { #if point to right of centroid
        quad <- 1
      } else { #point to left of centroid
        quad <- 2
      }
    } else { #point above centroid
      if(coords.i[1]>coords.cent[,1]) { #if point to right of centroid
        quad <- 4
      } else {
        quad <- 3
      }
    }
    q=c(q,quad)
  }
  c.j.ordered=sf::st_multipoint(coords.cj[order(q,angles),]) #order points by quadrat and slope (to/from centroid)

  c.j.fixed.list <- list()
  #move the points onto the edge of the polygon (i.e. out of the raster::buffer zone)
  for (i in 1:nrow(sf::st_coordinates(c.j.ordered))){
    c.j.fixed.list[[i]] <- sf::st_nearest_points(cell, sf::st_point(c.j.ordered[i,]))[[1]][1,]# finds the nearest point on polygon to wall junction
  }

  junction_points <- do.call("rbind", c.j.fixed.list)
  junction_points <- rbind( junction_points,  junction_points[1,]) #add first point to close the line
  rownames( junction_points) <- 1:nrow( junction_points)

  return(junction_points)
}
