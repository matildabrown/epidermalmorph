#' Find the number of stomatal rows
#'
#' @param cells SpatialPolygonsDataFrame, containing all cell polygons (i.e. not
#'  just the stomatal polygons). The specific colour of the stomata is irrelevant,
#'  but they should be the second-lowest value.
#'
#' @return A list with components:
#' \describe{
#' \item{n.rows}{Numeric. The number of rows of stomata}
#' \item{row.bounds}{Numeric. The y coordinates defining the edge of each row}
#' \item{row.wiggliness}{Numeric. The width of the stomatal band; see Details}
#' \item{row.consistency}{Numeric. The number of stomata per row.}
#' \item{stomrow.density}{Numeric. The number of rows of stomata per unit
#' length (per pixel if unscaled)}
#' }
#'
#' @details This function assumes that the pixels are coloured in the following order
#' (darkest to lightest): wall, stomate, subsidiary, pavement. To calculate
#' \code{row.wiggliness}, principal components analysis was computed for the
#' stomatal pixels in each row. The value \code{row.wiggliness} is the range
#' of PC2; it is recommended that this value be standardised based on stomatal
#'  complex width, to account for stomatal size variability.
#'
#'  @import stats
#'
#' @export

number_rows_stomata <- function(cells) {
#NUMBER OF STOMATAL ROWS
  boundbox <- sf::st_bbox(cells)
r <- raster::raster(xmn=boundbox[1], xmx=boundbox[3], ymn=boundbox[2], ymx=boundbox[4])
imstom <- fasterize::fasterize(cells,r, background=0, field=colnames(cells)[1])

t <- table(raster::values(imstom))

imstom[imstom!=as.numeric(names(t)[2])] <- NA
xy_stom <- data.frame(raster::rasterToPoints(imstom))
stomdens <- density(xy_stom$y, bw=30)
n.rows.stom <- length(localMaxima(stomdens$y))

stomdens2 <- stomdens$y*-1
stom_row_boundaries <- stomdens$x[localMaxima(stomdens2)]
stom_density_peaks <- stomdens$x[localMaxima(stomdens$y)]

stom_centroids <- suppressWarnings(
  sf::st_coordinates(sf::st_centroid(cells[which(cells$value==85),]))
)

stom_centroids_df <- data.frame(stom_centroids)
colnames(stom_centroids_df) <- c("x","y")

row.wiggliness <- NULL
row.consistency <- NULL
stom_centroids_df$rowID <- NA

for (k in 1:n.rows.stom){
  stom_centroids_df[which(stom_centroids_df$y>stom_row_boundaries[k] & stom_centroids_df$y<stom_row_boundaries[k+1]),"rowID"] <- k
  if(nrow(na.omit(stom_centroids_df[stom_centroids_df$rowID==k,1:2]))>1){
    fit <- prcomp(na.omit(stom_centroids_df[stom_centroids_df$rowID==k,1:2]))
    row.wiggliness <- c(row.wiggliness, summary(fit)$importance[2,2]*100)
    row.consistency <- c(row.consistency, nrow(na.omit(stom_centroids_df[stom_centroids_df$rowID==k,1:2])))
  }
}

stomrow.density <- (max(stom_density_peaks)-min(stom_density_peaks))/(length(stom_density_peaks)-1)

return (list(n.rows=n.rows.stom,
             row.bounds=stom_row_boundaries,
             row.wiggliness=row.wiggliness,
             row.consistency=row.consistency,
             dist.between.stom.rows= stomrow.density
))
}

