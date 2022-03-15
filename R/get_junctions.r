#' Clean up junction point coordinates
#'
#' @param wj Two-column numeric matrix with x and y coordinates of junction points.
#' @param buffer.width Numeric, the distance within which two cell junctions are
#' merged. Suggested value: 5 pixels (or equivalent if image is scaled)
#'
#' @return \code{sf MULTIPOINT} object of the cleaned junction points
#'
#'@export

get_junctions <- function(wj, buffer.width){

  wjp <- sf::st_multipoint(as.matrix(wj))

  #CLEAN UP JUNCTIONS
  #remove junctions within 5px of each other
  wall.nb.dist <- spdep::dnearneigh(wjp, 0,buffer.width)
 # summary(wall.nb.dist)
  wjpc <- as.data.frame(as.matrix(wjp))
  wjpc$keep <- NA
  for (i in 1:nrow(wjpc)){
    if(wall.nb.dist[[i]][1]==0){
      wjpc[i,"keep"] <- 1
    } else {
      if (min(c(i,wall.nb.dist[[i]]))==i) {
        wjpc[i,"keep"] <- 1
      } else {
        wjpc[i,"keep"] <- 0
      }
    }
  }

  wjp <- sf::st_multipoint(as.matrix(wjpc[which(wjpc$keep==1),1:2]))
  return(wjp)
}
