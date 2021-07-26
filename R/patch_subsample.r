#' Sample a connected patch of cells within an image

#' @param cells SpatialPolygonsDataFrame of all cells, with edge cells removed.
#' @param k.cells Numeric. Desired number of cells to be sampled in patch
#' @param k.pavement Numeric. Desired number of pavement cells to be sampled in patch
#' @param k.stomata Numeric. Desired number of stomata to be sampled in patch
#' @param plot.patch Logical. Plot the sampling iterations? Defaults to FALSE.
#' @param seed Numeric. For reproducibility, the starting cell can be set using
#' this argument.
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

patch_subsample <- function(cells, k.cells=NULL, k.pavement=NULL, k.stomata=NULL, plot.patch=FALSE, seed=NULL){

  #make sure there is a stopping value
  if (is.null(k.cells) & is.null(k.pavement) & is.null(k.stomata)){
    stop("At least one k-value must be provided.")
  }

  if (plot.patch==TRUE) {
    sp::plot(cells)
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
  blacklisted_cells <- which(lapply(nb_q, length)==1)
  #subsidiary cells without neighbouring stomatal cells

  #cells disconnected from rest of sample




  ############ set seed #### note that this is diffeent to set.seed()
  if(is.null(seed))  seed <- sample(1:length(cells),1)

  if (plot.patch==TRUE) {
    sp::plot(cells[seed,], add=T, col=pal[1])
  }

  #if it belongs to a stomatal complex, include the whole thing
  if (cells@data[seed,1]<200){
    if(cells@data[seed,1]==85) {
      stomID <- seed
    } else {
      stomID <- nb_q[[seed]][which(cells@data[nb_q[[seed]],1]==85)]
      if(length(stomID)>1){
        centroids <- data.frame(t(sapply(slot(cells[c(seed,stomID),], "polygons"), function(x) slot(x, "labpt"))))
        dists = as.matrix(dist(centroids))[-1,1]
        stomID <- stomID[which(dists==min(dists))]
      }
      if(length(stomID)==0){
        return(NULL)
      }
    }
    c0 <- c(stomID,nb_q[[stomID]][which(cells@data[nb_q[[stomID]],1]<200)])
  } else {
    c0 <- seed
  }

  if (plot.patch==TRUE) sp::plot(cells[c0,], add=T, col=pal[1])

  ######################### 1 away from seed ####

  c1.list <- grow_connected_region(cells,nb_q,seed=c0,blacklisted_cells=blacklisted_cells)
  ci <- c1.list[[1]]


  sampled_cells <- c(c0,ci)

  if (plot.patch==TRUE) sp::plot(cells[ci,], add=T, col=pal[2])

  #################### 2+ away from seed ####

  i=3
  while(
    length(sampled_cells)<k.cells ||
    length(which(cells@data[sampled_cells,1]==85))<k.stomata |
    length(which(cells@data[sampled_cells,1]==255))<k.pavement
  ){
    ci.list <- grow_connected_region(cells,nb_q,seed=ci,existing.patch=sampled_cells,blacklisted_cells=blacklisted_cells)
    ci <- ci.list[[1]]
    sampled_cells <- c(sampled_cells,ci)
    if (plot.patch==TRUE) sp::plot(cells[ci,], add=T, col=pal[i])
    i=i+1
  }

  return(sampled_cells)
}
