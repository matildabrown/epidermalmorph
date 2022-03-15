#' Extract measurements from sampled_data.frame and polygons
#'
#' @param x List of the form output by \code{extract_epidermal_traits()} (see documentation of that function for details)
#' @param polygons sf object containing the cell polygons that were used to generate sample_index
#' @param sample_index Numeric vector of cells to be measured.
#' @param trait.subset Character. Specific traits to measure. Defaults to NULL (all traits measured).
#' @param stom.val Numeric. Value of stomata. Defaults to 85.
#' @param subs.val Numeric. Value of subsidiary cells. Defaults to 170.
#' @param pave.val Numeric. Value of pavement cells. Defaults to 255.
#'
#' @return A vector with whole-image values.
#'
#' @details This function is designed to be used when the image has already
#' been processed and all cells have already been measured (e.g.
#' when evaluating sampling effort). For unmeasured images, see
#' \code{extract_epidermal_traits()}.
#'
#' @import stats
#'
#' @export



subsampling_trait_extraction <- function(x, polygons, sample_index,
                                         trait.subset=NULL,
                                         stom.val=85,
                                         subs.val=170,
                                         pave.val=255) {

  ID <- stom <- NULL

dat_im <- x$image_data
if(!is.null(trait.subset)){
  meta.cols <- c("image.ID","n","k.cells", "k.pavement","k.stomata",
                 "rotation","n.cells","n.pavement","n.stomata",
                 "n.subsidiary","n.pavezone", "n.stomzone","n.polar")
  dat_im <- dplyr::select(dat_im, dplyr::all_of(colnames(dat_im)[c(which(colnames(dat_im) %in% meta.cols),
                                      which(colnames(dat_im) %in% trait.subset))]))
}
cells <- polygons[polygons$ID %in% sample_index,]
#set up data frames

  cell_df_stom <- x$stomata_individual_data %>% dplyr::filter(ID %in% sample_index)
  cell_df_pavezone <- x$pavezone_individual_data %>% dplyr::filter(ID %in% sample_index)
  cell_df_stomzone <- x$stomzone_individual_data %>% dplyr::filter(ID %in% sample_index)
  cell_df_polar <- x$polar_individual_data %>% dplyr::filter(ID %in% sample_index)

  dat_im[2,"n.cells"] <- length(sample_index)
  dat_im[2,"n.stomata"] <- nrow(cell_df_stom)
  dat_im[2,"n.pavement"] <- nrow( cell_df_pavezone)+nrow(cell_df_stomzone)+nrow(cell_df_polar)
  dat_im[2,"n.subsidiary"] <- nrow(cells[which(cells$value==subs.val),])
  dat_im[2,"n.pavezone"] <- nrow( cell_df_pavezone)
  dat_im[2,"n.stomzone"] <- nrow( cell_df_stomzone)
  dat_im[2,"n.polar"] <- nrow( cell_df_polar)

  if(dat_im[2,"n.stomata"]>0){
  if("stomatal.index" %in% names(dat_im)) dat_im[2,"stomatal.index"] <- nrow(cells[which(cells$value==stom.val),])/nrow(cells)*100
  if("stomatal.density" %in% names(dat_im)) dat_im[2,"stomatal.density"] <- nrow(cells[which(cells$value==stom.val),]
  )/(sf::st_bbox(cells)[3]*sf::st_bbox(cells)[4])

  if(dat_im[2,"n.stomata"]>5){
  if("dist.between.stom.rows" %in% names(dat_im) |
     "row.wiggliness" %in% names(dat_im) |
     "row.consistency" %in% names(dat_im) ) {
    stomatal_rows <- number_rows_stomata(cells)
    if("dist.between.stom.rows" %in% names(dat_im)) dat_im[2,"dist.between.stom.rows"] <- stomatal_rows[[5]]
    if("row.wiggliness" %in% names(dat_im)) dat_im[2,"row.wiggliness"] <- mean(stomatal_rows[[3]])
    if("row.consistency" %in% names(dat_im)) dat_im[2,"row.consistency"] <- sd(stomatal_rows[[4]])
  }
  }

  # mean distance to kth nearest neighbours

  if("stom.distNN.mean" %in% names(dat_im) |
     "stom.distNN.mean" %in% names(dat_im) |
     "stom.distNN.mean" %in% names(dat_im) |
     "stom.dist2NN.sd" %in% names(dat_im)) {

    stomdist <- cell_cell_distances(stom, k=1:2)

    if("stom.distNN.mean" %in% names(dat_im)) dat_im[2,"stom.distNN.mean"] <- mean(stomdist[,1])
    if("stom.distNN.sd" %in% names(dat_im)) dat_im[2,"stom.distNN.sd"]<- sd(stomdist[,1])

    if("stom.dist2NN.mean" %in% names(dat_im)) dat_im[2,"stom.dist2NN.mean"] <- mean(stomdist[,2])
    if("stom.dist2NN.sd" %in% names(dat_im)) dat_im[2,"stom.dist2NN.sd"]<- sd(stomdist[,2])

  }

  if("stom.spacingNN.mean" %in% names(dat_im)) dat_im[2,"stom.spacingNN.mean"] <- mean(cells$dist_from_stom[which(cells$value==stom.val)],na.rm=TRUE)
  if("stom.spacingNN.sd" %in% names(dat_im)) dat_im[2,"stom.spacingNN.sd"] <- sd(cells$dist_from_stom[which(cells$value==stom.val)], na.rm=TRUE)

  if("stom.nsubcells.mean" %in% names(dat_im)) dat_im[2,"stom.nsubcells.mean"] <- nrow(cells[which(cells$value==subs.val),])/nrow(cells[which(cells$value==stom.val),])
  if("stom.subsarea.mean" %in% names(dat_im)) dat_im[2,"stom.subsarea.mean"] <-as.numeric(sum(sf::st_area(cells[which(cells$value==subs.val),])))/nrow(cells[which(cells$value==stom.val),])

  if("stom.gclength.mean" %in% names(dat_im)) dat_im[2,"stom.gclength.mean"] <-  mean(cell_df_stom$gclength, na.rm=T)
  if("stom.gclength.sd" %in% names(dat_im)) dat_im[2,"stom.gclength.sd"] <- sd(cell_df_stom$gclength, na.rm=T)
  if("stom.AR.mean" %in% names(dat_im)) dat_im[2,"stom.AR.mean"] <- mean(cell_df_stom$length/cell_df_stom$width, na.rm=T)
  if("stom.AR.sd" %in% names(dat_im)) dat_im[2,"stom.AR.sd"] <- sd(cell_df_stom$length/cell_df_stom$width, na.rm=T)
  if("stom.butterfly.mean" %in% names(dat_im)) dat_im[2,"stom.butterfly.mean"] <- mean(cell_df_stom$butterfly, na.rm=T)
  if("stom.symmetry.mean" %in% names(dat_im)) dat_im[2,"stom.symmetry.mean"] <- mean(cell_df_stom$symmetry, na.rm=T)
  if("stom.angle.sd" %in% names(dat_im)) dat_im[2,"stom.angle.sd"] <- sd(cell_df_stom$angle)
}


  if("pavezone.area.median" %in% names(dat_im)) dat_im[2,"pavezone.area.median"]  <- median(cell_df_pavezone$area)
  if("pavezone.area.sd" %in% names(dat_im)) dat_im[2,"pavezone.area.sd"]  <- sd(cell_df_pavezone$area)
  if("pavezone.AR.median" %in% names(dat_im))  dat_im[2,"pavezone.AR.median"]  <- median(cell_df_pavezone$AR)
  if("pavezone.AR.sd" %in% names(dat_im)) dat_im[2,"pavezone.AR.sd"]  <- sd(cell_df_pavezone$AR)
  if("pavezone.njunctionpts.mean" %in% names(dat_im))  dat_im[2,"pavezone.njunctionpts.mean"]  <- mean(cell_df_pavezone$n.junction.points)
  if("pavezone.njunctionpts.sd" %in% names(dat_im))  dat_im[2,"pavezone.njunctionpts.sd"]  <- sd(cell_df_pavezone$n.junction.points)
  if("pavezone.complexity.median" %in% names(dat_im))  dat_im[2,"pavezone.complexity.median"]  <- median(cell_df_pavezone$complexity)
  if("pavezone.complexity.sd" %in% names(dat_im))  dat_im[2,"pavezone.complexity.sd"]  <- sd(cell_df_pavezone$complexity)
  if("pavezone.undulation.amp.median" %in% names(dat_im))  dat_im[2,"pavezone.undulation.amp.median"]  <- median(cell_df_pavezone$undulation.amp)
  if("pavezone.undulation.amp.sd" %in% names(dat_im))  dat_im[2,"pavezone.undulation.amp.sd"]  <- sd(cell_df_pavezone$undulation.amp)
  if("pavezone.undulation.freq.mean" %in% names(dat_im))  dat_im[2,"pavezone.undulation.freq.mean"]  <- mean(cell_df_pavezone$undulation.freq)
  if("pavezone.undulation.freq.sd" %in% names(dat_im))  dat_im[2,"pavezone.undulation.freq.sd"]  <- sd(cell_df_pavezone$undulation.freq)
  if("pavezone.endwallangle.mean" %in% names(dat_im))  dat_im[2,"pavezone.endwallangle.mean"]  <- mean(cell_df_pavezone$endwall.angle.mean)
  if("pavezone.endwalldiff.mean" %in% names(dat_im))  dat_im[2,"pavezone.endwalldiff.mean"]  <- mean(cell_df_pavezone$endwall.angle.diff)
  if("pavezone.angle.median" %in% names(dat_im))  dat_im[2,"pavezone.angle.median"]  <- mean(cell_df_pavezone$angle)
  if("pavezone.angle.sd" %in% names(dat_im))  dat_im[2,"pavezone.angle.sd"]  <- sd(cell_df_pavezone$angle)

  if("stomzone.area.median" %in% names(dat_im)) dat_im[2,"stomzone.area.median"] <- median(cell_df_stomzone$`area`)
  if("stomzone.area.sd" %in% names(dat_im)) dat_im[2,"stomzone.area.sd"] <-  sd(cell_df_stomzone$`area`)
  if("stomzone.AR.median" %in% names(dat_im)) dat_im[2,"stomzone.AR.median"] <-  median(cell_df_stomzone$AR)
  if("stomzone.complexity.median" %in% names(dat_im)) dat_im[2,"stomzone.complexity.median"] <- median(cell_df_stomzone$complexity)
  if("stomzone.undulation.amp.median" %in% names(dat_im)) dat_im[2,"stomzone.undulation.amp.median"] <- median(cell_df_stomzone$undulation.amp)
  if("stomzone.undulation.freq.mean" %in% names(dat_im)) dat_im[2,"stomzone.undulation.freq.mean"] <- mean(cell_df_stomzone$undulation.freq)

  if("polar.area.median" %in% names(dat_im)) dat_im[2,"polar.area.median"] <- median(cell_df_polar$area)
  if("polar.area.sd" %in% names(dat_im))  dat_im[2,"polar.area.sd"] <-  sd(cell_df_polar$area)
  if("polar.AR.median" %in% names(dat_im)) dat_im[2,"polar.AR.median"] <-  median(cell_df_polar$AR)
  if("polar.complexity.median" %in% names(dat_im)) dat_im[2,"polar.complexity.median"] <- median(cell_df_polar$complexity)
  if("polar.undulation.amp.median" %in% names(dat_im)) dat_im[2,"polar.undulation.amp.median"] <- median(cell_df_polar$undulation.amp)
  if("polar.undulation.freq.mean" %in% names(dat_im)) dat_im[2,"polar.undulation.freq.mean"] <- mean(cell_df_polar$undulation.freq)


  return(dat_im)
  }

