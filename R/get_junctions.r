#' Clean up junction point coordinates
#'
#' @param wj Coordinates
#'
#' @return Some SpatialPoints
#'

get_junctions <- function(wj){

  wjp <- sp::SpatialPoints(wj)

  #CLEAN UP JUNCTIONS
  #remove junctions within 5px of each other
  wall.nb.dist <- spdep::dnearneigh(wjp, 0,5)
 # summary(wall.nb.dist)
  wjpc <- data.frame(wjp@coords)
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

  wjp <- sp::SpatialPoints(wjpc[which(wjpc$keep==1),1:2])
  return(wjp)
}
