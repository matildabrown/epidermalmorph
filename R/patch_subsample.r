#' Sample a connected patch of cells within an image

#' @param cells SpatialPolygonsDataFrame of all cells, with edge cells removed.
#' @param k.cells Numeric. Desired number of cells to be sampled in patch
#' @param k.pavement Numeric. Desired number of pavement cells to be sampled in patch
#' @param k.stomata Numeric. Desired number of stomata to be sampled in patch
#' @param plot.patch Logical. Plot the sampling iterations? Defaults to FALSE.
#' @param seed Numeric. For reproducibility, the starting cell can be set using
#' this argument.
#' @param stom.val Numeric vector of length 1. The value of stomatal cells (guard cells + pore)
#' @param subs.val Numeric vector of length 1. The value of subsidiary cells
#' @param pave.val Numeric vector of length 1. The value of pavement cells
#'
#' @return A numeric vector containing the indices of the sampled polygons
#'
#' @details If setting a seed using the \code{seed} argument, cell numbers can be
#' displayed on a plot of the cells using: \code{text(cells, labels = 1:length(cells))}.

#'
#' @examples
#' \dontrun{ #
#' #minimum patch size of 40 cells
#' patch_subsample(cells, k.cells=40)
#' #patch must include at least 10 stomatal complexes and 30 pavement cells
#' patch_subsample(cells,k.pavement=30, k.stomata=10)
#' }
#' @export

patch_subsample <- function(cells,
                            k.cells=NULL,
                            k.pavement=NULL,
                            k.stomata=NULL,
                            plot.patch=FALSE,
                            seed=NULL,
                            stom.val=85,
                            subs.val=170,
                            pave.val=255){

  #make sure there is a stopping value
  if (is.null(k.cells) & is.null(k.pavement) & is.null(k.stomata)){
    stop("At least one k-value must be provided as a stopping threshold.")
  }

  if (plot.patch==TRUE) {
    plot(sf::st_geometry(cells))
    colfunc <- grDevices::colorRampPalette(c("#6e010a", "#ebcccf"))
    if (!is.null(k.cells)) pal <- colfunc(ceiling(log10(k.cells)*3.5))
    if (is.null(k.cells)) {
      if(!is.null(k.stomata)) {pal <- colfunc(k.stomata)
      } else {pal <- colfunc(k.pavement/2) }
    }
  }


  if(is.null(k.cells)) k.cells <- 0
  if(is.null(k.stomata)) k.stomata <- 0
  if(is.null(k.pavement)) k.pavement <- 0

  nb_q <- spdep::poly2nb(cells, snap=5)
  blacklisted_cells <- which(lapply(nb_q, length)==1) #cells disconnected from rest of sample
  row.names(cells) <- cells$ID

  ############ set seed #### note that this is different to set.seed()
  if(is.null(seed))  seed <- sample(1:nrow(cells),1)

  if (plot.patch==TRUE) {
    plot(sf::st_geometry(cells[seed,]), add=T, col=pal[1])
  }

  #if it belongs to a stomatal complex, include the whole thing
  if (cells$value[seed]<200){
    if(cells$value[seed]==85) {
      stomID <- seed
    } else {
      stomID <- nb_q[[seed]][which(cells$value[nb_q[[seed]]]==stom.val)]
      if(length(stomID)>1){
        centroids <- suppressWarnings(sf::st_coordinates(sf::st_centroid(cells[c(seed,stomID),])))
        dists = as.matrix(dist(centroids))[-1,1]
        stomID <- stomID[which(dists==min(dists))]
      }
      if(length(stomID)==0){
        return(NULL)
      }
    }
    c0 <- c(stomID,nb_q[[stomID]][which(cells$value[nb_q[[stomID]]]<200)])
  } else {
    c0 <- seed
  }

  if (plot.patch==TRUE) plot(sf::st_geometry(cells[c0,]), add=T, col=pal[1])

  ######################### 1 away from seed ####

  c1.list <- grow_connected_region(cells,nb_q,seed=c0,blacklisted_cells=blacklisted_cells)
  ci <- c1.list[[1]]
  ci <- ci[which(ci %notin% c0)]


  sampled_cells <- c(c0,ci)

  if (plot.patch==TRUE) plot(sf::st_geometry(cells[ci,]), add=T, col=pal[2])

  #################### 2+ away from seed ####

  i=3
  while(
    length(sampled_cells)<k.cells ||
    length(which(cells$value[sampled_cells]==85))<k.stomata |
    length(which(cells$value[sampled_cells]==255))<k.pavement
  ){
    ci.list <- grow_connected_region(cells,nb_q,seed=ci,existing.patch=sampled_cells,blacklisted_cells=blacklisted_cells)
    ci <- ci.list[[1]]
    sampled_cells <- c(sampled_cells,ci)
    if (plot.patch==TRUE) plot(st_geometry(cells[ci,]), add=T, col=pal[i])
    i=i+1
  }

  return(cells$ID[sampled_cells])
}
