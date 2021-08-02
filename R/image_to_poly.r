#' Convert image to SpatialPolygons
#'
#' @param dir String. Directory for images.
#' @param file String. File name.
#'
#'#' @return A list with components:
#' \describe{
#' \item{cells}{SpatialPolygonsDataFrame. All cells, with values stored in a 'value' column.}
#' \item{junction_coords}{Dataframe with two columns: x and y coordinates of tri-cell junctions.}
#' }
#'
#' @import methods
#'
#' @export


image_to_poly <- function(dir, file){

  im_st <- stars::read_stars(paste0(dir, "/", file))[,]
  p <-  sf::st_as_sf(im_st, as_points = FALSE, merge = TRUE)
  values <- sf::st_drop_geometry(p)
  p <- as(p, "Spatial")

  colnames(p@data) <- "value"

  p.cells <- p[which(p$value!=0),]
  p.cells <- p.cells[which(rgeos::gArea(p.cells, byid = T)>10),]
  walls <- p[which(values==0),]

  im_st_walls <- im_st[[1]]
  im_st_walls[which(im_st_walls>0)] <- 100
  im_st_walls[which(im_st_walls==0)] <- 255
  im_st_walls[which(im_st_walls==100)] <- 0

  man <- imager::cimg2magick(imager::as.cimg(im_st_walls))
  jcn <- man %>%
   image_flop() %>%
   image_morphology('Thinning','Skeleton') %>%
   image_morphology('Thinning','Corners') %>%
   image_morphology('Thinning','LineEnds', iterations=-1) %>%
   image_morphology('Thinning','Skeleton') %>%
   image_morphology('Thinning','Diagonals') %>%
   image_morphology('Thinning','Corners') %>%
   image_morphology("HMT","LineJunctions")



  jcn_data <- t(matrix(as.integer(jcn[[1]]), nrow=attributes(im_st)$dimensions$y$to))
  #jcn_data[which(skel==0)] <- 0



#plot(imager::as.cimg(jcn_points), interpolate = F)
jcncoords <- setNames(reshape2::melt(jcn_data),c("x","y","val"))
jcncoords <- jcncoords[which(jcncoords$val!=0 & jcncoords$x<=attributes(im_st)$dimensions$x$to),1:2]
jcncoords$y <- attributes(im_st)$dimensions$y$to - jcncoords$y

  return(list(cells=p.cells, junction_coords=jcncoords) )
}



