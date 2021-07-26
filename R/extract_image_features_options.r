#' Extract cell measurements from an image
#'
#' @param image.file Character. File name of image.
#' @param wall.skeleton.file Character. File name of corresponding cell wall
#' skeleton image (see Details).
#' @param dir Character. Directory for images. Defaults to NULL.
#' @param stomata.present Logical. Does the image contain stomata? Defaults
#' to TRUE, see details.
#' @param stomatal.shape Logical. Include measurements of stomatal shape?
#' Defaults to TRUE, see details.
#' @param stomatal.arrangement Logical. Include measurements of stomatal
#'  arrangement? Defaults to TRUE, see details.
#' @param pavement.cells Logical. Measure pavement cells? Defaults to TRUE,
#'  see details.
#' @param sd.measures Logical. Return standard deviation of variables?
#'  Defaults to TRUE, see details.
#' @param specific.inclusions Character. Cell measurements to include;
#'  see details for options.
#' @param specific.exclusions Character. Cell measurements to exclude;
#'  see details for options.
#' @param verbose Logical. Should progress output be printed in the console?
#' Defaults to FALSE.
#'
#' @return A list of length three, containing a vector with whole-image values, a
#' data.frame of pavement cell values and a data.frame of stomatal values.
#'
#' @details This function contains the full pipeline from image to measurements.
#' The labelled image is required, along with the cell wall skeleton (with
#' junction points labelled), such as that produced by ImageJ's 'AnalyzeSkeleton'
#' function. If only a subset of variables is desired, this can can be set using the
#' arguments \code{ stomatal.shape, stomatal.arrangement, pavement.cells, sd.measures,
#' specific.inclusions,  specific.exclusions}. For a full list of variables and
#' graphical descriptions of what they measure, see the vignette here. ((add link))
#'
#' @importFrom graphics par
#' @import methods
#' @import stats
#'
#' @export



extract_image_features_options <- function(image.file, wall.skeleton.file, dir=".",
                                   stomata.present = TRUE,
                                   stomatal.shape = TRUE,
                                   stomatal.arrangement = TRUE,
                                   pavement.cells = TRUE,
                                   sd.measures = TRUE,
                                   specific.inclusions = NULL,
                                   specific.exclusions = NULL,
                                   verbose=FALSE) {


basic_names <- c(
  "image.ID")
stomatal_basic_names <- c(
  "stomatal.north",
  "stomatal.index",
  "stomatal.density.px2")
stom_arr_names <- c(
  "dist.between.stom.rows",
  "row.wiggliness",
  "row.consistency",
  "stom.m.DNN",
  "stom.m.D2NN",
  "m.stom.stom.dist",
  "stom.angle.sd")
stom_shape_names <- c(
  "stom.nsubcells.mean",
  "stom.subsarea.mean",
  "stom.gclength.mean",
  "stom.AR.mean",
  "stom.butterfly.mean",
  "stom.symmetry.mean")
pavement_names <- c(
  "pavezone.area.median",
  "pavezone.AR.median",
  "pavezone.njunctionpts.mean",
  "pavezone.complexity.median",
  "pavezone.maxamplitude.median",
  "pavezone.ncrosses.mean",
  "pavezone.endwallangle.mean",
  "pavezone.endwalldiff.mean",
  "pavezone.angle.median",
  "stomzone.area.median",
  "stomzone.AR.median",
  "stomzone.maxamplitude.median",
  "stomzone.ncrosses.mean",
  "stomzone.complexity.median",
  "pavezone.complexity.sd",
  "pavezone.maxamplitude.sd",
  "pavezone.ncrosses.sd",
  "polar.area.median",
  "polar.AR.median",
  "polar.complexity.median",
  "polar.maxamplitude.median",
  "polar.ncrosses.mean")

stom_arr_sd_names <- c(
  "stom.sd.DNN",
  "stom.sd.D2NN",
  "var.stom.stom.dist")
stom_shape_sd_names <- c(
  "stom.gclength.sd",
  "stom.AR.sd")
pavement_sd_names <- c(
  "pavezone.area.sd.as.pc.median",
  "pavezone.AR.sd",
  "pavezone.njunctionpts.sd",
    "pavezone.angle.sd",
  "stomzone.area.sd",
  "polar.area.sd")

#set up the variables
selected_names <- basic_names

if(stomata.present==TRUE) {
  selected_names <- c(selected_names, stomatal_basic_names)
} else {
  stomatal.arrangement <- FALSE
  stomatal.shape <- FALSE
}

if(stomatal.arrangement==T) {
  selected_names <- c(selected_names, stom_arr_names)
  if(sd.measures==T) selected_names <- c(selected_names, stom_arr_sd_names)
}

if(stomatal.shape==T) {
  selected_names <- c(selected_names, stom_shape_names)
  if(sd.measures==T) selected_names <- c(selected_names, stom_shape_sd_names)
}

if(pavement.cells==T) {
  if(stomata.present==T) {
  selected_names <- c(selected_names, pavement_names)
  } else {
    selected_names <- c(selected_names, pavement_names[1:9])
  }
  if(sd.measures==T) {
    if(stomata.present==T) {
      selected_names <- c(selected_names, pavement_sd_names)
  } else {
    selected_names <- c(selected_names, pavement_sd_names[1:5])
  }
  }
}

if(!is.null(specific.inclusions)) selected_names <- c(selected_names, specific.inclusions)

if(!is.null(specific.exclusions)) selected_names <- selected_names[which(selected_names %notin% specific.exclusions)]



dat_im <- rep(NA,times=length(selected_names))
names(dat_im) <- selected_names

dat_im["image.ID"] <- image.file
  #Read in image, separate cell types ####
  cells <- epidermalmorph::image_to_poly(dir, image.file)
  pave <- cells[[1]]
  stom <- cells[[2]]
  subs <- cells[[3]]
  allcells <- raster::bind(pave, stom, subs)

  #wall junctions
  wj <- raster::raster(paste0("../ImageJ_results/Skeletons/", image.file))
  wjp <- get_junctions(wj)

  #Identify edge cells ####
  edges <- epidermalmorph::find_edge_cells(dir, image.file, allcells)

  allcells@data$edge <- 0
  allcells@data$edge[edges[[1]]] <- 1
  allcells@data$lose <- 0
  allcells@data$lose[edges[[2]]] <- 1

  # Find stomatal north, rotate sp::Polygons accordingly ####
  stom_edge <- find_edge_cells(dir, image.file, stom)
  stom_angles = NULL
  for (j in 1:length(stom)){
    if (j %notin% stom_edge[[1]]){
      c.smooth <- smoothr::smooth(stom[j,], method = "ksmooth", smoothness = 3)
      ####ellipse bit
      coords_for_ellipse <- c.smooth@polygons[[1]]@Polygons[[1]]@coords
      ellipDirect <- conicfit::EllipseDirectFit(coords_for_ellipse)
      ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG

      angle <- ellipDirectG[5]
      if(angle>(pi/2)) angle=0-pi+angle
      stom_angles <- c(stom_angles, angle)
    }
  }
  stomatal_north <- pracma::rad2deg(mean(stom_angles))

  # align to stomatal north
  pave <- maptools::elide(pave, rotate=stomatal_north, center=rowMeans(sp::bbox(allcells)))
  stom <- maptools::elide(stom, rotate=stomatal_north, center=rowMeans(sp::bbox(allcells)))
  subs <- maptools::elide(subs, rotate=stomatal_north, center=rowMeans(sp::bbox(allcells)))
  allcells <- maptools::elide(allcells, rotate=stomatal_north, center=rowMeans(sp::bbox(allcells)))
  wjp <- maptools::elide(wjp, rotate=stomatal_north, center=rowMeans(sp::bbox(allcells)))

  # Sort out pave zone/stom zone/polar cells ####

  if(stomata.present==TRUE) {
  paths <- cell_graph_shortest_paths(allcells)
  allcells@data$dist_from_stom <- shortest_path_to_cell_type(allcells,paths,85) #85 is stomatal value

  stom_keep <- allcells[which(allcells@data[,1]==85 & allcells@data$edge==0),]
  pavezone_pave <- allcells[which(allcells@data[,1]==255 & allcells@data$edge==0 & allcells$dist_from_stom>2),]
  stomzone_pave <- allcells[which(allcells@data[,1]==255 & allcells@data$edge==0 & allcells$dist_from_stom==2),]
  polar_pave <- allcells[which(allcells@data[,1]==255 & allcells@data$edge==0 & allcells$dist_from_stom==1),]
  subs_keep <- allcells[which(allcells@data[,1]==170 & allcells@data$edge==0),]
  pave_keep <- allcells[which(allcells@data[,1]==255 & allcells@data$edge==0),]

  } else {
    pave_keep <- allcells[which(allcells@data[,1]==255 & allcells@data$edge==0),]
    pavezone_pave <- pave_keep
  }

  ##################################################################################
  # MEASUREMENTS ----
  if(verbose==TRUE) print("Extracting cell measurements...")
  # Stomatal measurements ####

  if("stomatal.north" %in% names(dat_im)) dat_im["stomatal.north"] <- stomatal_north

  if("stomatal.index" %in% names(dat_im)) dat_im["stomatal.index"] <- nrow(allcells@data[which(allcells$lose==0 & allcells@data[,1]==85),])/nrow(allcells@data[which(allcells$lose==0 & allcells@data[,1]!=85),])*100
  if("stomatal.density.px2" %in% names(dat_im)) dat_im["stomatal.density.px2"] <- nrow(allcells@data[which(allcells$lose==0 & allcells@data[,1]==85),])/(sp::bbox(allcells)[1,2]*sp::bbox(allcells)[2,2])

  if("dist.between.stom.rows" %in% names(dat_im) |
     "row.wiggliness" %in% names(dat_im) |
     "row.consistency" %in% names(dat_im) ) {
  stomatal_rows <- number_rows_stomata(allcells[which(allcells@data$edge==0),])
  if("dist.between.stom.rows" %in% names(dat_im)) dat_im["dist.between.stom.rows"] <- stomatal_rows[[5]]
  if("row.wiggliness" %in% names(dat_im)) dat_im["row.wiggliness"] <- mean(stomatal_rows[[3]])
  if("row.consistency" %in% names(dat_im)) dat_im["row.consistency"] <- sd(stomatal_rows[[4]])
  }

  # mean distance to kth nearest neighbours


  if("stom.m.DNN" %in% names(dat_im) |
     "stom.sd.DNN" %in% names(dat_im) |
     "stom.m.D2NN" %in% names(dat_im) |
     "stom.sd.D2NN" %in% names(dat_im)) {

  stomdist <- cell_cell_distances(stom, k=1:2)

  if("stom.m.DNN" %in% names(dat_im)) dat_im["stom.m.DNN"] <- mean(stomdist[,1])
  if("stom.sd.DNN" %in% names(dat_im)) dat_im["stom.sd.DNN"]<- sd(stomdist[,1])/mean(stomdist[,1])

  if("stom.m.D2NN" %in% names(dat_im)) dat_im["stom.m.D2NN"] <- mean(stomdist[,2])
  if("stom.sd.D2NN" %in% names(dat_im)) dat_im["stom.sd.D2NN"]<- sd(stomdist[,2])/mean(stomdist[,2])

  }

  if("m.stom.stom.dist" %in% names(dat_im)) dat_im["m.stom.stom.dist"] <- mean(allcells@data$dist_from_stom[which(allcells@data[,1]==85)])
  if("var.stom.stom.dist" %in% names(dat_im)) dat_im["var.stom.stom.dist"] <- var(allcells@data$dist_from_stom[which(allcells@data[,1]==85)])



  if("stom.nsubcells.mean" %in% names(dat_im)) dat_im["stom.nsubcells.mean"] <- nrow(allcells@data[which(allcells$lose==0 & allcells@data[,1]==170),])/nrow(allcells@data[which(allcells$lose==0 & allcells@data[,1]==85),])
  if("stom.subsarea.mean" %in% names(dat_im)) dat_im["stom.subsarea.mean"] <- rgeos::gArea(subs_keep)/length(stom_keep)


  # stomatal complex measurements ####

  if(stomata.present==T){

  if(verbose==TRUE) print("Calculating stomatal features...")

  bstom  <- raster::buffer(stom_keep, width=5, dissolve=F)
  cell_df_stom <- data.frame(ID=as.numeric(sapply(stom_keep@polygons, function(x) slot(x, name = "ID"))))
  cell_df_stom$ID <- as.numeric(as.character(cell_df_stom$ID))
  cell_df_stom[,2:3] <- rgeos::gCentroid(stom_keep, byid=T)@coords; colnames(cell_df_stom)[2:3] <- c("centroid.x", "centroid.y")

  stom_comp_key <- rgeos::gIntersects(bstom,subs,byid=TRUE)


  for (j in 1:nrow(cell_df_stom)){
    if(verbose==TRUE) print(paste0(image.file,": stomatal complex shape : ",j))

    stom_j <- stom_keep[j,]
    bstom_j <- bstom[j,]

    subs_j <- subs[which(stom_comp_key[,j]==TRUE),]
    stom_comp_j <- raster::bind(stom_j, subs_j)

    if(length(find_edge_cells(dir, image.file, stom_comp_j, rotated=stomatal_north)$edge)==0){

      if(length(subs_j)>2){
        bsubs <- raster::buffer(subs_j, width=2, dissolve=F)
        subs_key <-  rgeos::gIntersects(bsubs, byid=T)
        subs1 <- which(subs_key[1,]==TRUE)
        subs2 <- which(1:length(subs_j) %notin% subs1)
        if(length(subs2)==0){
          subskey2 <- subs_key
          for (m in 1:(nrow(subskey2)-1)){
            for(n in (m+1):nrow(subskey2)){
              if(subskey2[m,n]==T)  subskey2[m,n] <- rgeos::gArea(rgeos::gIntersection(bsubs[m,], bsubs[n,]))
              if(subskey2[m,n]==F)  subskey2[m,n] <- 0
            }
          }
          subs1 <- as.numeric(which(subskey2==max(subskey2), arr.ind = T))
          subs2 <- as.numeric(which(1:length(subs_j) %notin% subs1))
        }
      }


      coords_for_ellipse <- stom_j@polygons[[1]]@Polygons[[1]]@coords
      ellipDirect <- conicfit::EllipseDirectFit(coords_for_ellipse)
      ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG

      angle <- ellipDirectG[5]
      if(angle>(pi/2)) angle=0-pi+angle

      stom_comp_j <- maptools::elide(stom_comp_j, rotate=angle-180)

      cell_df_stom[j,"gclength"] <- stom_j@bbox[1,2]-stom_comp_j@bbox[1,1]

      cell_df_stom[j,"length"] <- stom_comp_j@bbox[1,2]-stom_comp_j@bbox[1,1]
      cell_df_stom[j,"width"] <- stom_comp_j@bbox[2,2]-stom_comp_j@bbox[2,1]
      cell_df_stom[j,"butterfly"] <- (stom_j@bbox[1,2]-stom_j@bbox[1,1])/(stom_comp_j@bbox[1,2]-stom_comp_j@bbox[1,1])
      if(length(subs_j)==2) cell_df_stom[j,"symmetry"] <- rgeos::gArea(subs_j[1,])/rgeos::gArea(subs_j[2,])
      if(length(subs_j)>2) cell_df_stom[j,"symmetry"] <- rgeos::gArea(subs_j[subs1,])/rgeos::gArea(subs_j[subs2,])
      if(length(subs_j)<2) {cell_df_stom[j,"symmetry"] <- NA
      }else{if(cell_df_stom[j,"symmetry"] >1) cell_df_stom[j,"symmetry"]<- 1/cell_df_stom[j,"symmetry"]}


    } else {
      if(verbose==TRUE) print(paste("stomatal complex on the edge - not measured"))
      cell_df_stom[j,] <- NA
    }

  }


  if("stom.gclength.mean" %in% names(dat_im)) dat_im["stom.gclength.mean"] <-  mean(cell_df_stom$gclength, na.rm=T)
  if("stom.gclength.sd" %in% names(dat_im)) dat_im["stom.gclength.sd"] <- sd(cell_df_stom$gclength, na.rm=T)
  if("stom.AR.mean" %in% names(dat_im)) dat_im["stom.AR.mean"] <- mean(cell_df_stom$length/cell_df_stom$width, na.rm=T)
  if("stom.AR.sd" %in% names(dat_im)) dat_im["stom.AR.sd"] <- sd(cell_df_stom$length/cell_df_stom$width, na.rm=T)
  if("stom.butterfly.mean" %in% names(dat_im)) dat_im["stom.butterfly.mean"] <- mean(cell_df_stom$butterfly, na.rm=T)
  if("stom.symmetry.mean" %in% names(dat_im)) dat_im["stom.symmetry.mean"] <- mean(cell_df_stom$symmetry, na.rm=T)
  if("stom.angle.sd" %in% names(dat_im)) dat_im["stom.angle.sd"] <- sd(stom_angles)

  }


  if(pavement.cells==TRUE)
 { # Pavement cell measurements ####
  if(length(pavezone_pave)>0){
    # pavezone ####

    cell_df_pavezone <- data.frame(ID=sapply(pavezone_pave@polygons, function(x) methods::slot(x, name = "ID")))
    cell_df_pavezone$ID <- as.numeric(as.character(cell_df_pavezone$ID))

    #vectorised measurements
    cell_df_pavezone[,2:3] <- rgeos::gCentroid(pavezone_pave, byid=T)@coords; colnames(cell_df_pavezone)[2:3] <- c("centroid.x", "centroid.y")
    cell_df_pavezone[,"area.px^2"] <- rgeos::gArea(pavezone_pave, byid=T)
    cell_df_pavezone[,"perimeter.px"] <- rgeos::gLength(rgeos::gBoundary(pavezone_pave, byid=T), byid=T)

    for (j in 1:nrow(cell_df_pavezone)){

      if(verbose==TRUE) print(paste0(image.file," : pavezone : ",j))

      cell_j <- pavezone_pave[j,]

      cell_df_pavezone[j,"AR"] <- (cell_j@bbox[1,2]-cell_j@bbox[1,1])/(cell_j@bbox[2,2]-cell_j@bbox[2,1])

      c.smooth<- smoothr::smooth(cell_j, method = "ksmooth", smoothness = 3)

      #GET THE JUNCTION POINTS THAT DEFINE THE SIMPLIFIED CELL
      junction_points <- cell.simplify(cell=c.smooth, cell.junctions=wjp)

      #make sure no double ups
      if(nrow(unique(junction_points))<(nrow(junction_points)-1)){
        n <- nrow(junction_points)
        junction_points <- rbind(junction_points[1,], unique(junction_points[2:(n-1),]), junction_points[n,])
      }
      cell_df_pavezone[j,"n.junction.points"] <- nrow(junction_points)-1



      p = sp::Polygon(junction_points)
      ps = sp::Polygons(list(p),1)
      straightpoly = sp::SpatialPolygons(list(ps))
      straightlines = raster::spLines(junction_points)
      actuallines = as(c.smooth, 'SpatialLines')


      #GET difference between straight perimeter and actual perimeter
      cell_df_pavezone[j,"complexity"] <- rgeos::gLength(actuallines)/rgeos::gLength(straightlines)

      #GET MAX DISTANCE
      cell_df_pavezone[j,"max.dist"] <- max_amplitude(c.smooth, junction_points)

      #where does the actual boundary cross the simplified boundary?
      crosses = rgeos::gIntersection(straightlines, actuallines)
      cell_df_pavezone[j,"n.crosses"] <- length(crosses)

      #ANGLE OF ENDPOINT WALLS

      endwalls <- endwall_angles(junction_points)

      cell_df_pavezone[j,"endwall.angle.mean"] <- mean(c(endwalls[[3]],endwalls[[4]]))
      cell_df_pavezone[j,"endwall.angle.diff"] <- abs(endwalls[[3]]-endwalls[[4]])

      #angle of cell
      coords_for_ellipse <- c.smooth@polygons[[1]]@Polygons[[1]]@coords
      ellipDirect <- conicfit::EllipseDirectFit(coords_for_ellipse)
      ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG

      angle <- ellipDirectG[5]
      if(angle>(pi/2)) angle=0-pi+angle
      cell_df_pavezone[j,"angle"] <- angle



    }


    if("pavezone.area.median" %in% names(dat_im)) dat_im["pavezone.area.median"]  <- median(cell_df_pavezone$`area.px^2`)
    if("pavezone.area.sd.as.pc.median" %in% names(dat_im)) dat_im["pavezone.area.sd.as.pc.median"]  <- sd(cell_df_pavezone$`area.px^2`)/median(cell_df_pavezone$`area.px^2`)
    if("pavezone.AR.median" %in% names(dat_im))  dat_im["pavezone.AR.median"]  <- median(cell_df_pavezone$AR)
    if("pavezone.AR.sd" %in% names(dat_im)) dat_im["pavezone.AR.sd"]  <- sd(cell_df_pavezone$AR)
    if("pavezone.njunctionpts.mean" %in% names(dat_im))  dat_im["pavezone.njunctionpts.mean"]  <- mean(cell_df_pavezone$n.junction.points)
    if("pavezone.njunctionpts.sd" %in% names(dat_im))  dat_im["pavezone.njunctionpts.sd"]  <- sd(cell_df_pavezone$n.junction.points)
    if("pavezone.complexity.median" %in% names(dat_im))  dat_im["pavezone.complexity.median"]  <- median(cell_df_pavezone$complexity)
    if("pavezone.complexity.sd" %in% names(dat_im))  dat_im["pavezone.complexity.sd"]  <- sd(cell_df_pavezone$complexity)
    if("pavezone.maxamplitude.median" %in% names(dat_im))  dat_im["pavezone.maxamplitude.median"]  <- median(cell_df_pavezone$max.dist)
    if("pavezone.maxamplitude.sd" %in% names(dat_im))  dat_im["pavezone.maxamplitude.sd"]  <- sd(cell_df_pavezone$max.dist)
    if("pavezone.ncrosses.mean" %in% names(dat_im))  dat_im["pavezone.ncrosses.mean"]  <- mean(cell_df_pavezone$n.crosses)
    if("pavezone.ncrosses.sd" %in% names(dat_im))  dat_im["pavezone.ncrosses.sd"]  <- sd(cell_df_pavezone$n.crosses)
    if("pavezone.endwallangle.mean" %in% names(dat_im))  dat_im["pavezone.endwallangle.mean"]  <- mean(cell_df_pavezone$endwall.angle.mean)
    if("pavezone.endwalldiff.mean" %in% names(dat_im))  dat_im["pavezone.endwalldiff.mean"]  <- mean(cell_df_pavezone$endwall.angle.diff)
    if("pavezone.angle.median" %in% names(dat_im))  dat_im["pavezone.angle.median"]  <- mean(cell_df_pavezone$angle)
    if("pavezone.angle.sd" %in% names(dat_im))  dat_im["pavezone.angle.sd"]  <- sd(cell_df_pavezone$angle)

  } else {
    if("pavezone.area.median" %in% names(dat_im))  dat_im["pavezone.area.median"]  <- NA
    if("pavezone.area.sd.as.pc.median" %in% names(dat_im))  dat_im["pavezone.area.sd.as.pc.median"]  <- NA
    if("pavezone.AR.median" %in% names(dat_im))  dat_im["pavezone.AR.median"]  <- NA
    if("pavezone.AR.sd" %in% names(dat_im))  dat_im["pavezone.AR.sd"]  <- NA
    if("pavezone.njunctionpts.mean" %in% names(dat_im))  dat_im["pavezone.njunctionpts.mean"]  <- NA
    if("pavezone.njunctionpts.sd" %in% names(dat_im))  dat_im["pavezone.njunctionpts.sd"]  <- NA
    if("pavezone.complexity.median" %in% names(dat_im))  dat_im["pavezone.complexity.median"]  <- NA
    if("pavezone.complexity.sd" %in% names(dat_im))  dat_im["pavezone.complexity.sd"]  <- NA
    if("pavezone.maxamplitude.median" %in% names(dat_im))  dat_im["pavezone.maxamplitude.median"]  <- NA
    if("pavezone.maxamplitude.sd" %in% names(dat_im))  dat_im["pavezone.maxamplitude.sd"]  <- NA
    if("pavezone.ncrosses.mean" %in% names(dat_im))  dat_im["pavezone.ncrosses.mean"]  <- NA
    if("pavezone.ncrosses.sd" %in% names(dat_im))  dat_im["pavezone.ncrosses.sd"]  <- NA
    if("pavezone.endwallangle.mean" %in% names(dat_im))  dat_im["pavezone.endwallangle.mean"]  <- NA
    if("pavezone.endwalldiff.mean" %in% names(dat_im))  dat_im["pavezone.endwalldiff.mean"]  <- NA
    if("pavezone.angle.median" %in% names(dat_im))  dat_im["pavezone.angle.median"]  <- NA
    if("pavezone.angle.sd" %in% names(dat_im))  dat_im["pavezone.angle.sd"]  <- NA
  }

  if(stomata.present==TRUE) {
  # stomzone ####
  cell_df_stomzone <- data.frame(ID=sapply(stomzone_pave@polygons, function(x) slot(x, name = "ID")))
  cell_df_stomzone$ID <- as.numeric(as.character(cell_df_stomzone$ID))

  #vectorised measurements
  cell_df_stomzone[,2:3] <- rgeos::gCentroid(stomzone_pave, byid=T)@coords; colnames(cell_df_stomzone)[2:3] <- c("centroid.x", "centroid.y")
  cell_df_stomzone[,"area.px^2"] <- rgeos::gArea(stomzone_pave, byid=T)
  cell_df_stomzone[,"perimeter.px"] <- rgeos::gLength(rgeos::gBoundary(stomzone_pave, byid=T), byid=T)

  for (j in 1:nrow(cell_df_stomzone)){

    if(verbose==TRUE) print(paste0(image.file," : stomzone : ",j))

    cell_j <- stomzone_pave[j,]

    cell_df_stomzone[j,"AR"] <- (cell_j@bbox[1,2]-cell_j@bbox[1,1])/(cell_j@bbox[2,2]-cell_j@bbox[2,1])

    c.smooth<- smoothr::smooth(cell_j, method = "ksmooth", smoothness = 3)

    #GET THE JUNCTION POINTS THAT DEFINE THE SIMPLIFIED CELL
    junction_points <- cell.simplify(cell=c.smooth, cell.junctions=wjp)

    #make sure no double ups
    if(nrow(unique(junction_points))<(nrow(junction_points)-1)){
      n <- nrow(junction_points)
      junction_points <- rbind(junction_points[1,], unique(junction_points[2:(n-1),]), junction_points[n,])
    }


      p = sp::Polygon(junction_points)
      ps = sp::Polygons(list(p),1)
      straightpoly = sp::SpatialPolygons(list(ps))
      straightlines = raster::spLines(junction_points)
      actuallines = as(c.smooth, 'SpatialLines')



    #GET difference between straight perimeter and actual perimeter
    cell_df_stomzone[j,"complexity"] <- rgeos::gLength(actuallines)/rgeos::gLength(straightlines)

    #GET MAX DISTANCE
    cell_df_stomzone[j,"max.dist"] <- max_amplitude(c.smooth, junction_points)

    #where does the actual boundary cross the simplified boundary?
    crosses = rgeos::gIntersection(straightlines, actuallines)
    cell_df_stomzone[j,"n.crosses"] <- length(crosses)
  }

  if("stomzone.area.median" %in% names(dat_im)) dat_im["stomzone.area.median"] <- median(cell_df_stomzone$`area.px^2`)
  if("stomzone.area.sd" %in% names(dat_im)) dat_im["stomzone.area.sd"] <-  sd(cell_df_stomzone$`area.px^2`)
  if("stomzone.AR.median" %in% names(dat_im)) dat_im["stomzone.AR.median"] <-  median(cell_df_stomzone$AR)
  if("stomzone.complexity.median" %in% names(dat_im)) dat_im["stomzone.complexity.median"] <- median(cell_df_stomzone$complexity)
  if("stomzone.maxamplitude.median" %in% names(dat_im)) dat_im["stomzone.maxamplitude.median"] <- median(cell_df_stomzone$max.dist)
  if("stomzone.ncrosses.mean" %in% names(dat_im)) dat_im["stomzone.ncrosses.mean"] <- mean(cell_df_stomzone$n.crosses)

  # polar ####

  cell_df_polar <- data.frame(ID=sapply(polar_pave@polygons, function(x) slot(x, name = "ID")))
  cell_df_polar$ID <- as.numeric(as.character(cell_df_polar$ID))

  #vectorised measurements
  cell_df_polar[,2:3] <- rgeos::gCentroid(polar_pave, byid=T)@coords; colnames(cell_df_polar)[2:3] <- c("centroid.x", "centroid.y")
  cell_df_polar[,"area.px^2"] <- rgeos::gArea(polar_pave, byid=T)
  cell_df_polar[,"perimeter.px"] <- rgeos::gLength(rgeos::gBoundary(polar_pave, byid=T), byid=T)


  for (j in 1:nrow(cell_df_polar)){

    if(verbose==TRUE) print(paste0(image.file," : polar : ",j))

    cell_j <- polar_pave[j,]

    cell_df_polar[j,"AR"] <- (cell_j@bbox[1,2]-cell_j@bbox[1,1])/(cell_j@bbox[2,2]-cell_j@bbox[2,1])

    c.smooth<- smoothr::smooth(cell_j, method = "ksmooth", smoothness = 3)

    #GET THE JUNCTION POINTS THAT DEFINE THE SIMPLIFIED CELL
    junction_points <- cell.simplify(cell=c.smooth, cell.junctions=wjp)

    #make sure no double ups
    if(nrow(unique(junction_points))<(nrow(junction_points)-1)){
      n <- nrow(junction_points)
      junction_points <- rbind(junction_points[1,], unique(junction_points[2:(n-1),]), junction_points[n,])
    }


      p = sp::Polygon(junction_points)
      ps = sp::Polygons(list(p),1)
      straightpoly = sp::SpatialPolygons(list(ps))
      straightlines = raster::spLines(junction_points)
      actuallines = as(c.smooth, 'SpatialLines')


    #GET difference between straight perimeter and actual perimeter
    cell_df_polar[j,"complexity"] <- rgeos::gLength(actuallines)/rgeos::gLength(straightlines)

    #GET MAX DISTANCE
    cell_df_polar[j,"max.dist"] <- max_amplitude(c.smooth, junction_points)

    #where does the actual boundary cross the simplified boundary?
    crosses = rgeos::gIntersection(straightlines, actuallines)
    cell_df_polar[j,"n.crosses"] <- length(crosses)
  }

  if("polar.area.median" %in% names(dat_im)) dat_im["polar.area.median"] <- median(cell_df_polar$`area.px^2`)
  if("polar.area.sd" %in% names(dat_im))  dat_im["polar.area.sd"] <-  sd(cell_df_polar$`area.px^2`)
  if("polar.AR.median" %in% names(dat_im)) dat_im["polar.AR.median"] <-  median(cell_df_polar$AR)
  if("polar.complexity.median" %in% names(dat_im)) dat_im["polar.complexity.median"] <- median(cell_df_polar$complexity)
  if("polar.maxamplitude.median" %in% names(dat_im)) dat_im["polar.maxamplitude.median"] <- median(cell_df_polar$max.dist)
  if("polar.ncrosses.mean" %in% names(dat_im)) dat_im["polar.ncrosses.mean"] <- mean(cell_df_polar$n.crosses)

  } else {
    cell_df_stomzone <- NULL
    cell_df_polar <- NULL
    cell_df_stom <- NULL
}
}
  return(list(dat_im, cell_df_pavezone, cell_df_stomzone, cell_df_polar, cell_df_stom, allcells))
  }
