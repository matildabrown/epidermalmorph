#' Grow a patch of cells
#'
#' @param cells sf object of all cells, with edge cells removed.
#' @param nb Object of class "nb", created from the 'spdep' function \code{poly2nb()}.
#' @param seed Numeric. Index of first cell in patch (the 'seed' cell). Defaults to NULL.
#' @param existing.patch Numeric. Indices of cells in the patch which is to
#' be 'grown'. Defaults to NULL.
#' @param blacklisted_cells Numeric. Indices of cells to exclude (e.g. the subsidiary
#' cell of a stomatal complex that is cut off by the image edge).
#'
#' @return A numeric vector containing the indices of the sampled polygons
#'
#' @export



grow_connected_region <- function(cells, nb, seed=NULL, existing.patch=NULL,blacklisted_cells){


  if (is.null(seed)){
  seed <- sample(1:length(cells),1)
}

c.i <- NULL
for (m in 1:length(seed)){
  c.i <- c(c.i,nb[[seed[m]]])
}

c.i <- unique(c.i)
c.i <- c.i[which(c.i %notin% existing.patch)]

if (length(which(cells$value[c.i]<200))>0){
  stomcells <- c.i[which(cells$value[c.i]<200)]

  stomID <- NULL
  for (h in 1:length(stomcells)){
    if(cells$value[stomcells[h]]==85) {
      stomIDh <- stomcells[h]
    } else {
      stomIDh <- nb[[stomcells[h]]][which(cells$value[nb[[stomcells[h]]]]==85)]
    }
    stomID <- unique(c(stomID,stomIDh))
  }

  if(length(stomID)==0){
    blacklisted_cells <- c(blacklisted_cells, stomcells)
  } else {
    for (h in 1:length(stomID)){
      c.i <- unique(c(c.i,nb[[stomID[h]]][which(cells$value[nb[[stomID[h]]]]<200)] ))
    }
  }
  c.i <- unique(c(c.i,stomID))
}

c.i <- c.i[which(c.i %notin% blacklisted_cells)]
c.i <- c.i[which(c.i %notin% existing.patch)]

return(list(c.i,blacklisted_cells))
}
