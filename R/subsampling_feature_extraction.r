#' Extract measurements from sampled_data.frame and polygons
#'
#' @param sampled_polygons SpatialPolygonsDataFrame, containing the cell polygons.
#' @param sampled_data Data.frame, containing the measurements for all cells (some empty
#'  columns because all cell types included in one file)
#'
#' @return A vector with whole-image values.
#'
#' @details This function is designed to be used when the image has already
#' been processed and all cells have already been measured (e.g.
#' when evaluating sampling effort). For unmeasured images, see
#' \code{extract_image_features()}.
#'
#' @import stats
#'
#' @export



subsampling_feature_extraction <- function(sampled_polygons, sampled_data) {


  dat_im <- rep(NA,times=55)
  names(dat_im) <- c("image.ID","n.cells", "n.stom","n.pave","n.pavezone","n.stomzone",
                     "n.polar","stomatal.index",
                     "stomatal.row.density","row.wiggliness","row.consistency",
                     "stom.m.DNN","stom.sd.DNN","stom.m.D2NN","stom.sd.D2NN",
                     "stom.nsubcells.mean","stom.subsarea.mean","pavezone.area.median",
                     "pavezone.area.sd.as.pc.median","pavezone.AR.median","pavezone.AR.sd",
                     "pavezone.njunctionpts.mean","pavezone.njunctionpts.sd",
                     "pavezone.complexity.median","pavezone.complexity.sd",
                     "pavezone.maxamplitude.median","pavezone.maxamplitude.sd",
                     "pavezone.ncrosses.mean","pavezone.ncrosses.sd",
                     "pavezone.endwallangle.mean","pavezone.endwalldiff.mean",
                     "pavezone.angle.median","stomzone.area.median","stomzone.area.sd",
                     "stomzone.AR.median","stomzone.complexity.median",
                     "stomzone.maxamplitude.median","stomzone.ncrosses.mean",
                     "polar.area.median","polar.area.sd","polar.AR.median",
                     "polar.complexity.median","polar.maxamplitude.median",
                     "polar.ncrosses.mean","stom.gclength.mean","stom.gclength.sd",
                     "stom.AR.mean","stom.AR.sd","stom.butterfly.mean",
                     "stom.symmetry.mean","stom.angle.sd","stomatal.density.px2",
                     "m.stom.stom.dist","var.stom.stom.dist",
                     "pavezone.angle.sd")

  dat_stom <- sampled_data[which(sampled_data$type=="stom"),]
  dat_pave <- sampled_data[which(sampled_data$type=="pavezone"),]
  dat_stz <- sampled_data[which(sampled_data$type=="stomzone"),]
  dat_polar <- sampled_data[which(sampled_data$type=="polar"),]
  dat_pave.all <- sampled_data[which(sampled_data$type != "stom"),]

  dat_im["n.cells"] <- length(sampled_polygons)
  dat_im["n.stom"] <- nrow(dat_stom)
  dat_im["n.pave"] <- nrow(dat_pave.all)
  dat_im["n.pavezone"] <- nrow(dat_pave)
  dat_im["n.stomzone"] <- nrow(dat_stz)
  dat_im["n.polar"] <- nrow(dat_polar)

  dat_im["stomatal.index"] <- nrow(dat_stom)/(nrow(dat_pave)+nrow(dat_stz)+nrow(dat_polar))
  dat_im["stomatal.density.px2"] <- nrow(dat_stom)/rgeos::gArea(sampled_polygons, byid = FALSE)

  # There is a bug in the rasterize() function called by number_rows_stomata() so I am skipping them for now
  stomatal_rows <- number_rows_stomata(sampled_polygons)
  dat_im["stomatal.row.density"] <- stomatal_rows[[5]]
  dat_im["row.wiggliness"] <- mean(stomatal_rows[[3]])
  dat_im["row.consistency"] <- sd(stomatal_rows[[4]])

  # dat_im["stomatal.row.density"] <- NA
  # dat_im["row.wiggliness"] <- NA
  # dat_im["row.consistency"] <- NA

  # mean distance to kth nearest neighbours
  stomdist <- cell_cell_distances(sampled_polygons[which(sampled_polygons@data[,1]=="85"),], k=1:2)

  dat_im["stom.m.DNN"] <- mean(stomdist[,1])
  dat_im["stom.sd.DNN"]<- sd(stomdist[,1])/mean(stomdist[,1])

  dat_im["stom.m.D2NN"] <- mean(stomdist[,2])
  dat_im["stom.sd.D2NN"]<- sd(stomdist[,2])/mean(stomdist[,2])

  dat_im["m.stom.stom.dist"] <- mean(sampled_polygons@data$dist_from_stom[which(sampled_polygons@data[,1]==85)])
  dat_im["var.stom.stom.dist"] <- var(sampled_polygons@data$dist_from_stom[which(sampled_polygons@data[,1]==85)])

  dat_im["stom.nsubcells.mean"] <- nrow(sampled_polygons@data[which(sampled_polygons$lose==0 & sampled_polygons@data[,1]==170),])/nrow(sampled_polygons@data[which(sampled_polygons$lose==0 & sampled_polygons@data[,1]==85),])
  dat_im["stom.subsarea.mean"] <- rgeos::gArea(sampled_polygons[which(sampled_polygons@data[,1]==170),])/length(sampled_polygons[which(sampled_polygons@data[,1]==85),])

  dat_im["stom.gclength.mean"] <-  mean(dat_stom$gclength, na.rm=T)
  dat_im["stom.gclength.sd"] <- sd(dat_stom$gclength, na.rm=T)
  dat_im["stom.AR.mean"] <- mean(dat_stom$length/dat_stom$width, na.rm=T)
  dat_im["stom.AR.sd"] <- sd(dat_stom$length/dat_stom$width, na.rm=T)
  dat_im["stom.butterfly.mean"] <- mean(dat_stom$butterfly, na.rm=T)
  dat_im["stom.symmetry.mean"] <- mean(dat_stom$symmetry, na.rm=T)
  dat_im["stom.angle.sd"] <- sd(dat_stom$angle, na.rm=T)

  # Pavement cell measurements ####
  cell_df_pavezone <- sampled_data[which(sampled_data$type=="pavezone"),]
  if(nrow(cell_df_pavezone)>0){


    dat_im["pavezone.area.median"]  <- median(cell_df_pavezone$`area.px.2`)
    dat_im["pavezone.area.sd.as.pc.median"]  <- sd(cell_df_pavezone$`area.px.2`)/median(cell_df_pavezone$`area.px.2`)
    dat_im["pavezone.AR.median"]  <- median(cell_df_pavezone$AR)
    dat_im["pavezone.AR.sd"]  <- sd(cell_df_pavezone$AR)
    dat_im["pavezone.njunctionpts.mean"]  <- mean(cell_df_pavezone$n.junction.points)
    dat_im["pavezone.njunctionpts.sd"]  <- sd(cell_df_pavezone$n.junction.points)
    dat_im["pavezone.complexity.median"]  <- median(cell_df_pavezone$complexity)
    dat_im["pavezone.complexity.sd"]  <- sd(cell_df_pavezone$complexity)
    dat_im["pavezone.maxamplitude.median"]  <- median(cell_df_pavezone$max.dist)
    dat_im["pavezone.maxamplitude.sd"]  <- sd(cell_df_pavezone$max.dist)
    dat_im["pavezone.ncrosses.mean"]  <- mean(cell_df_pavezone$n.crosses)
    dat_im["pavezone.ncrosses.sd"]  <- sd(cell_df_pavezone$n.crosses)
    dat_im["pavezone.endwallangle.mean"]  <- mean(cell_df_pavezone$endwall.angle.mean)
    dat_im["pavezone.endwalldiff.mean"]  <- mean(cell_df_pavezone$endwall.angle.diff)
    dat_im["pavezone.angle.median"]  <- mean(cell_df_pavezone$angle)
    dat_im["pavezone.angle.sd"]  <- sd(cell_df_pavezone$angle)

  } else {
    dat_im["pavezone.area.median"]  <- NA
    dat_im["pavezone.area.sd.as.pc.median"]  <- NA
    dat_im["pavezone.AR.median"]  <- NA
    dat_im["pavezone.AR.sd"]  <- NA
    dat_im["pavezone.njunctionpts.mean"]  <- NA
    dat_im["pavezone.njunctionpts.sd"]  <- NA
    dat_im["pavezone.complexity.median"]  <- NA
    dat_im["pavezone.complexity.sd"]  <- NA
    dat_im["pavezone.maxamplitude.median"]  <- NA
    dat_im["pavezone.maxamplitude.sd"]  <- NA
    dat_im["pavezone.ncrosses.mean"]  <- NA
    dat_im["pavezone.ncrosses.sd"]  <- NA
    dat_im["pavezone.endwallangle.mean"]  <- NA
    dat_im["pavezone.endwalldiff.mean"]  <- NA
    dat_im["pavezone.angle.median"]  <- NA
    dat_im["pavezone.angle.sd"]  <- NA
  }
  # stomzone ####
  cell_df_stomzone <- sampled_data[which(sampled_data$type=="stomzone"),]
  dat_im["stomzone.area.median"] <- median(cell_df_stomzone$`area.px.2`)
  dat_im["stomzone.area.sd"] <-  sd(cell_df_stomzone$`area.px.2`)
  dat_im["stomzone.AR.median"] <-  median(cell_df_stomzone$AR)
  dat_im["stomzone.complexity.median"] <- median(cell_df_stomzone$complexity)
  dat_im["stomzone.maxamplitude.median"] <- median(cell_df_stomzone$max.dist)
  dat_im["stomzone.ncrosses.mean"] <- mean(cell_df_stomzone$n.crosses)

  # polar ####
  cell_df_polar <- sampled_data[which(sampled_data$type=="polar"),]
  dat_im["polar.area.median"] <- median(cell_df_polar$`area.px.2`)
  dat_im["polar.area.sd"] <-  sd(cell_df_polar$`area.px.2`)
  dat_im["polar.AR.median"] <-  median(cell_df_polar$AR)
  dat_im["polar.complexity.median"] <- median(cell_df_polar$complexity)
  dat_im["polar.maxamplitude.median"] <- median(cell_df_polar$max.dist)
  dat_im["polar.ncrosses.mean"] <- mean(cell_df_polar$n.crosses)

  return(dat_im)
  }

