#' Convert image to SpatialPolygons
#'
#' @param dir String. Directory for images.
#' @param fname String. File name.
#'
#'#' @return A list with components:
#' \describe{
#' \item{pave}{SpatialPolygonsDataFrame. Pavement cells. }
#' \item{stom}{SpatialPolygonsDataFrame. Stomata (pore + guard cells).}
#' \item{subs}{SpatialPolygonsDataFrame. Lateral subsidiary cells. }
#' }
#' @return List of SpatialPolygonsDataFrames: pavement cells, subsidiary cells and stomata.
#' @details This function assumes that the pixels are coloured in the following order
#' (darkest to lightest): wall, stomate, subsidiary, pavement.
#'
#' @import methods
#'
#' @export


image_to_poly <- function(dir, fname){

  im_st <- stars::read_stars(paste0(dir, "/", fname))[,]
  p <-  sf::st_as_sf(im_st, as_points = FALSE, merge = TRUE)
  values <- sf::st_drop_geometry(p)
  table(values)
  t <-  table(values)

  # t <- t[order(t)]

  paveno <- which(values==as.numeric(names(t)[4]))
  stomno <- which(values==as.numeric(names(t)[2]))
  subsno <- which(values==as.numeric(names(t)[3]))
  p1 <- p[paveno,]
  p2 <- p[stomno,]
  p3 <- p[subsno,]
  pave <- as(p1, "Spatial")
  stom <- as(p2, "Spatial")
  subs <- as(p3, "Spatial")
  pave <- pave[which(rgeos::gArea(pave, byid = T)>5),] #remove rogue pixels
  stom <- stom[which(rgeos::gArea(stom, byid = T)>5),] #remove rogue pixels
  subs <- subs[which(rgeos::gArea(subs, byid = T)>5),] #remove rogue pixels

  return(list(pave=pave, stom=stom, subs=subs))
}
