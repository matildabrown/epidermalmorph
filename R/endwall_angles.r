#' Find the endwall (furthest left and right) angles of cells
#'
#' @description
#' For a set of junction points that describe a simplified cell, this function
#' extracts the left-most wall and the right-most wall and calculates their
#' slope and bearing (in radians).
#'
#'@param junction_points This should be a matrix of the junction points defining
#'the polygon. The first and last points are the same (i.e. the line formed by
#'the points is closed).
#'@return A list with the following components:
#'#' \describe{
#' \item{slope1}{The slope of the leftmost endwall. Useful for plotting. }
#' \item{slope2}{The slope of the rightmost endwall.}
#' \item{angle1}{The angle of the leftmost endwall. Useful for analysis. }
#' \item{angle2}{The angle of the rightmost endwall.}
#' \item{points1}{The points defining the leftmost endwall.}
#' \item{points2}{The points defining the rightmost endwall.}
#' }
#' @details The \code{angle} outputs have been transformed to be within -90 and
#'  90, so that the distribution is approximately normal around 0.
#'  @export




endwall_angles <- function(junction_points){

  endpoint1.1 <- which.min(junction_points[,1])
  if(endpoint1.1-1==0) endpoint1.2 <- c(endpoint1.1+1, (nrow(junction_points)-1))[which.min(junction_points[c(endpoint1.1+1, (nrow(junction_points)-1)),1])]
  if(endpoint1.1+1>(nrow(junction_points)-1)) endpoint1.2 <- c(1, endpoint1.1-1)[which.min(junction_points[c(1, endpoint1.1-1),1])]
  if(endpoint1.1-1!=0 & endpoint1.1+1<=(nrow(junction_points)-1)) endpoint1.2 <- c(endpoint1.1+1, endpoint1.1-1)[which.min(junction_points[c(endpoint1.1+1, endpoint1.1-1),1])]

  x1 <- junction_points[endpoint1.2,1]
  x2 <- junction_points[endpoint1.1,1]
  y1 <- junction_points[endpoint1.2,2]
  y2 <- junction_points[endpoint1.1,2]

  s1 <- (y1-y2)/(x1-x2)

  slope1 <- atan(s1)
  if(slope1<0) slope1 <-pi+slope1

  endpoint2.1 <- which.max(junction_points[,1])
  if(endpoint2.1-1==0) endpoint2.2 <- c(endpoint2.1+1, (nrow(junction_points)-1))[which.max(junction_points[c(endpoint2.1+1, (nrow(junction_points)-1)),1])]
  if(endpoint2.1+1>(nrow(junction_points)-1)) endpoint2.2 <- c(1, endpoint2.1-1)[which.max(junction_points[c(1, endpoint2.1-1),1])]
  if(endpoint2.1-1!=0 & endpoint2.1+1<=(nrow(junction_points)-1)) endpoint2.2 <- c(endpoint2.1+1, endpoint2.1-1)[which.max(junction_points[c(endpoint2.1+1, endpoint2.1-1),1])]

  x3 <- junction_points[endpoint2.2,1]
  x4 <- junction_points[endpoint2.1,1]
  y3 <- junction_points[endpoint2.2,2]
  y4 <- junction_points[endpoint2.1,2]

  s2 <- (y3-y4)/(x3-x4)

  slope2 <- atan(s2)
  if(slope2<0) slope2 <-pi+slope2


  return(list(slope1=s1,slope2=s2,angle1=slope1,angle2=slope2,
              points1=matrix(c(x1,x2,y1,y2), nrow=2, dimnames=list(c("1","2"), c("x","y"))),
              points2=matrix(c(x3,x4,y3,y4), nrow=2, dimnames=list(c("1","2"), c("x","y")))))
}
