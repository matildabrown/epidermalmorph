#' Extract cell measurements from an image
#'
#' @param cells SpatialPolygonsDataFrame. The image to be measured, converted to polygons (e.g. using \code{image_to_poly()})
#' @param junction.points Either a numeric matrix or SpatialPoints with the coordinates of the cell wall junctions (e.g. as output by \code{image_to_poly()})
#' @param cells.present Character. A vector containing the names of all the types of cells present. Defaults to \code{c("pavement", "stomata", "subsidiary")}
#' @param cell.values A vector containing the values of all the types of cells present (usually numeric), see details.
#' @param NA.value A vector of length 1, with the value for NA regions of the image (e.g. damaged or un-traced regions)
#' @param stomatal.shape Logical. Include measurements of stomatal shape?
#' Defaults to TRUE, see details.
#' @param stomatal.arrangement Logical. Include measurements of stomatal
#'  arrangement? Defaults to TRUE, see details.
#' @param sd.measures Logical. Return standard deviation of variables?
#'  Defaults to TRUE, see details.
#' @param specific.inclusions Character. Cell measurements to include;
#'  see details for options.
#' @param specific.exclusions Character. Cell measurements to exclude;
#'  see details for options.
#' @param sd.as.percent.of.mean Logical. Should variability measures be presented as % of the mean? Defaults to FALSE.
#' @param verbose Logical. Should progress output be printed in the console?
#' Defaults to FALSE.
#'
#' @return A list of length three, containing a vector with whole-image values, a
#' data.frame of pavement cell values and a data.frame of stomatal values.
#'
#' @details This function contains the script for all measurements. Each type
#' of cell (e.g. pavement, stomate) should have a name and corresponding value
#' given in \code{cells.present} and \code{cell.values}, respectively. If
#' there are unnamed values in the image, these will be measured as
#' 'othertype1', 'othertype2' and so forth. Only limited numbers of
#' measurements will be taken for these cell types (index, density, area,
#' spacing and distance to nearest neighbour). If only a subset of variables is desired,
#' this can can be set using the arguments \code{ stomatal.shape, stomatal.arrangement, pavement.cells, sd.measures,
#' specific.inclusions,  specific.exclusions}. For a full list of variables and
#' graphical descriptions of what they measure, see the vignette here. ((add link)).
#'
#' @importFrom graphics par
#' @import methods
#' @import stats
#'
#' @export


extract_epidermal_traits <- function(cell.polygons, junction.points,
                                           cells.present = c("pavement", "stomata","subsidiary"),
                                           cell.values,
                                           NA.value=NULL,
                                           stomatal.shape = TRUE,
                                           stomatal.arrangement = TRUE,
                                           sd.measures = TRUE,
                                           specific.inclusions = NULL,
                                           specific.exclusions = NULL,
                                           sd.as.percent.of.mean = FALSE,
                                           verbose=FALSE) {
# SETUP ----
#Pre-empt some common input errors ####
  if (length(cells.present) != length(cell.values)) stop("Need to supply the same number of cell values as cell types.")
  if (length(unique(cell.values))!=length(cell.values)) stop("Some cell values repeated.")
  if (length(unique(cells.present))!=length(cells.present)) stop("Some cell types repeated.")
  if(!is.null(NA.value)) if(NA.value %in% cell.values) stop("NA.value muct be different from defined cell type values.")

  if (length(setdiff(unique(cell.polygons$value),cell.values))>0) warning("Not all values in cells have been allocated a type. These cells will be measured as 'other.type.1', 'other.type.2', etc.")

  if(!is.null(NA.value)) cells=cell.polygons[which(cell.polygons$value!=NA.value),] else cells=cell.polygons

#Identify edge cells ####
edges <- find_edge_cells(cell.polygons,cells)

cells@data$edge <- 0
cells@data$edge[edges[[1]]] <- 1
cells@data$edge.notcounted <- 0
cells@data$edge.notcounted[edges[[2]]] <- 1

# Clean up junction points ####
wjp <- get_junctions(junction.points)

# Find stomatal north, rotate sp::Polygons accordingly ####
if ("stomata" %in% cells.present){
  stom.val <- cell.values[which(cells.present=="stomata")]
  stom <- cells[which(cells$value==stom.val),]

stom_angles = NULL
for (j in 1:length(stom)){
  if (stom$edge[j]==0){
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
cells <- maptools::elide(cells, rotate=stomatal_north, center=rowMeans(sp::bbox(cells)))
wjp <- maptools::elide(wjp, rotate=stomatal_north, center=rowMeans(sp::bbox(cells)))@coords
}

# Split up cells into cell types ####

#set up stomata (note: stomata re-extracted from rotated polygons)
  if ("stomata" %in% cells.present){
    stom.val <- cell.values[which(cells.present=="stomata")]
    stom <- cells[which(cells$value==stom.val),]
    stom_keep  <- stom[which(stom$edge==0),]
    # Sort out pave zone/stom zone/polar cells ####
      paths <- cell_graph_shortest_paths(cells)
      cells@data$dist_from_stom <- shortest_path_to_cell_type(cells, paths, stom.val)

  }

#set up subsidiary cells
  if ("subsidiary" %in% cells.present){
    subs.val <- cell.values[which(cells.present=="subsidiary")]
    subs <- cells[which(cells$value==subs.val),]
    subs_keep <- subs[which(subs$edge==0),]
  }

#set up pavement cells, split up polygons by distance fom stomata
if ("pavement" %in% cells.present){
  pave.val <- cell.values[which(cells.present=="pavement")]
  pave <- cells[which(cells$value==pave.val),]
  pave_keep <- pave[which(pave$edge==0),]

  #Split up pavement cells based on distance to stomata
  if ("stomata" %in% cells.present){
  pavezone <- pave[which(pave@data$edge==0 & pave$dist_from_stom>2),]
  stomzone <- pave[which(pave@data$edge==0 & pave$dist_from_stom==2),]
  polar <- pave[which(pave@data$edge==0 & pave$dist_from_stom==1),]

} else {
  pavezone <- pave_keep #if there are no stomata, all pavement cells are 'pavezone'
}
}

  if (length(setdiff(cells.present, c("pavement", "stomata","subsidiary")))>0){
    othertypes <- data.frame(celltype=setdiff(cells.present, c("pavement", "stomata","subsidiary")), val=NA)
    for (i in 1:nrow(othertypes)){
      othertypes[i,2] <- cell.values[which(cells.present==othertypes[i,1])]
    }
  } else { othertypes <- NULL }

alltypes <- unique(rbind(data.frame(celltype=cells.present, val=cell.values), othertypes))

if (length(setdiff(unique(cells$value),cell.values)>0)){
  warning("Not all values in cells have been allocated a type. These cells will be
          measured as 'othertype1', 'othertype2', etc.")
  othertypes <- data.frame(celltype=NA, val=setdiff(unique(cells$value),cell.values))
  for (i in 1:nrow(othertypes)){
    othertypes[i,1] <- paste0("othertype",i)
  }
}

# Now set up empty data.frame ####
basic_names <- c(
  "image.ID",
  "rotation",
  paste0("n.",alltypes[,1])
  )

stomatal_basic_names <- c(
  "stomatal.density.px2",
  "stomatal.index"
    )

stom_arrangement_names <- c(
  "dist.between.stom.rows",
  "row.consistency",
  "row.wiggliness",
  "stom.angle.sd",
  "stom.distNN.mean",
  "stom.dist2NN.mean",
  "stom.spacingNN.mean"
  )
stom_shape_names <- c(
  "stom.nsubcells.mean",
  "stom.subsarea.mean",
  "stom.gclength.mean",
  "stom.AR.mean",
  "stom.butterfly.mean",
  "stom.symmetry.mean")

pavement_names <- c(
    "pavezone.angle.median",
    "pavezone.AR.median",
    "pavezone.area.median",
    "pavezone.complexity.median",
    "pavezone.endwallangle.mean",
    "pavezone.endwalldiff.mean",
    "pavezone.njunctionpts.mean",
    "pavezone.undulation.amp.median",
    "pavezone.undulation.freq.mean",

    "stomzone.AR.median",
    "stomzone.area.median",
    "stomzone.complexity.median",
    "stomzone.undulation.amp.median",
    "stomzone.undulation.freq.mean",

    "polar.AR.median",
    "polar.area.median",
    "polar.complexity.median",
    "polar.undulation.amp.median",
    "polar.undulation.freq.mean")

other_cell_type_names <- c(
  "area.median",
  "density",
  "distNN.mean", "dist2NN.mean",
  "index",
  "spacingNN.mean"
  )

stom_arr_sd_names <- c(
  "stom.distNN.sd",
  "stom.dist2NN.sd",
  "stom_spacingNN.sd")
stom_shape_sd_names <- c(
  "stom.gclength.sd",
  "stom.AR.sd")
pavement_sd_names <- c(
  "pavezone.area.sd",
  "pavezone.AR.sd",
  "pavezone.njunctionpts.sd",
  "pavezone.angle.sd",
  "stomzone.area.sd",
  "polar.area.sd",
  "pavezone.complexity.sd",
  "pavezone.undulation.amp.sd",
  "pavezone.undulation.freq.sd")
other_cell_type_sd_names <- c(
  "area.sd",
  "distNN.sd", "dist2NN.sd",
  "spacingNN.sd"
)

#set up the variables
selected_names <- basic_names

if ("stomata" %in% cells.present) {
  selected_names <- c(selected_names, stomatal_basic_names)
} else {
  stomatal.arrangement <- FALSE
  stomatal.shape <- FALSE
}

if(stomatal.arrangement==T) {
  selected_names <- c(selected_names, stom_arrangement_names)
  if(sd.measures==T) selected_names <- c(selected_names, stom_arr_sd_names)
}

if(stomatal.shape==T) {
  selected_names <- c(selected_names, stom_shape_names)
  if(sd.measures==T) selected_names <- c(selected_names, stom_shape_sd_names)
}

if ("pavement" %in% cells.present)  {
  if ("stomata" %in% cells.present)  {
    selected_names <- c(selected_names, pavement_names)
  } else {
    selected_names <- c(selected_names, pavement_names[1:9])
  }
  if(sd.measures==T) {
    if ("stomata" %in% cells.present)  {
      selected_names <- c(selected_names, pavement_sd_names)
    } else {
      selected_names <- c(selected_names, pavement_sd_names[1:5])
    }
  }
}

if(!is.null(othertypes)){
  for (i in 1:nrow(othertypes)){
    selected_names <- c(selected_names, paste0(othertypes[i,1],".", other_cell_type_names))
  }
}

if(!is.null(specific.inclusions)) selected_names <- c(selected_names, specific.inclusions)

if(!is.null(specific.exclusions)) selected_names <- selected_names[which(selected_names %notin% specific.exclusions)]


dat_im <- rep(NA,times=length(selected_names))
names(dat_im) <- selected_names

####


# MEASUREMENTS ----
if(verbose==TRUE) print("Extracting cell measurements...")
# Stomatal arrangement measurements ####

#stomatal north
if("rotation" %in% names(dat_im)) dat_im["rotation"] <- stomatal_north

#number of each cell type (NOT ON THE EDGE)
n.cell.types <- table(cells$value[which(cells$edge==0)])
for (i in 1:nrow(alltypes)){
  dat_im[paste0("n.",alltypes[i,1])] <- as.numeric(n.cell.types[which(names(n.cell.types)==alltypes[i,2])])
}

if("stomatal.index" %in% names(dat_im)) dat_im["stomatal.index"] <- nrow(cells@data[which(cells$edge.notcounted==0 &
                                                                                            cells@data[,1]==stom.val),]
                                                                         )/nrow(cells@data[which(cells$edge.notcounted==0),])*100
if("stomatal.density.px2" %in% names(dat_im)) dat_im["stomatal.density.px2"] <- nrow(cells@data[which(cells$edge.notcounted==0 &
                                                                                                        cells@data[,1]==stom.val),]
                                                                                     )/(sp::bbox(cells)[1,2]*sp::bbox(cells)[2,2])

if("dist.between.stom.rows" %in% names(dat_im) |
   "row.wiggliness" %in% names(dat_im) |
   "row.consistency" %in% names(dat_im) ) {
  stomatal_rows <- number_rows_stomata(cells[which(cells@data$edge==0),])
  if("dist.between.stom.rows" %in% names(dat_im)) dat_im["dist.between.stom.rows"] <- stomatal_rows[[5]]
  if("row.wiggliness" %in% names(dat_im)) dat_im["row.wiggliness"] <- mean(stomatal_rows[[3]])
  if("row.consistency" %in% names(dat_im)) dat_im["row.consistency"] <- sd(stomatal_rows[[4]])
}

# mean distance to kth nearest neighbours

if("stom.distNN.mean" %in% names(dat_im) |
   "stom.distNN.mean" %in% names(dat_im) |
   "stom.distNN.mean" %in% names(dat_im) |
   "stom.dist2NN.sd" %in% names(dat_im)) {

  stomdist <- cell_cell_distances(stom, k=1:2)

  if("stom.distNN.mean" %in% names(dat_im)) dat_im["stom.distNN.mean"] <- mean(stomdist[,1])
  if("stom.distNN.sd" %in% names(dat_im)) dat_im["stom.distNN.sd"]<- sd(stomdist[,1])

  if("stom.dist2NN.mean" %in% names(dat_im)) dat_im["stom.dist2NN.mean"] <- mean(stomdist[,2])
  if("stom.dist2NN.sd" %in% names(dat_im)) dat_im["stom.dist2NN.sd"]<- sd(stomdist[,2])

}

if("stom.spacingNN.mean" %in% names(dat_im)) dat_im["stom.spacingNN.mean"] <- mean(cells@data$dist_from_stom[which(cells@data[,1]==stom.val)])
if("stom.spacingNN.sd" %in% names(dat_im)) dat_im["stom.spacingNN.sd"] <- sd(cells@data$dist_from_stom[which(cells@data[,1]==stom.val)])

if("stom.nsubcells.mean" %in% names(dat_im)) dat_im["stom.nsubcells.mean"] <- nrow(cells@data[which(cells$edge.notcounted ==0 & cells@data[,1]==subs.val),])/nrow(cells@data[which(cells$edge.notcounted==0 & cells@data[,1]==stom.val),])
if("stom.subsarea.mean" %in% names(dat_im)) dat_im["stom.subsarea.mean"] <- rgeos::gArea(subs_keep)/length(stom_keep)


# Stomatal shape measurements ####
if("stomata" %in% cells.present){

  if(verbose==TRUE) print("Calculating stomatal features...")

  bstom  <- raster::buffer(stom_keep, width=5, dissolve=F)
  cell_df_stom <- data.frame(ID=as.numeric(sapply(stom_keep@polygons, function(x) slot(x, name = "ID"))))
  cell_df_stom$ID <- as.numeric(as.character(cell_df_stom$ID))
  cell_df_stom[,2:3] <- rgeos::gCentroid(stom_keep, byid=T)@coords; colnames(cell_df_stom)[2:3] <- c("centroid.x", "centroid.y")

  stom_comp_key <- rgeos::gIntersects(bstom,subs,byid=TRUE)


  for (j in 1:nrow(cell_df_stom)){
    if(verbose==TRUE & stomatal.shape==T & (j == round(j/10)*10| j==length(stom))) print(paste0("Measuring stomatal complex shape : ",j,"/",length(stom)))

    stom_j <- stom_keep[j,]
    bstom_j <- bstom[j,]

    subs_j <- subs[which(stom_comp_key[,j]==TRUE),]
    stom_comp_j <- raster::bind(stom_j, subs_j)

    if(length(find_edge_cells(cell.polygons, stom_comp_j, rotated=stomatal_north)$edge)==0){

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

} else cell_df_stom <- NULL

if(verbose==TRUE) print("Calculating pavement cell features...")
# Pavezone (>2 cells away from stomata) ####
if("pavement" %in% cells.present) {

  if(verbose==TRUE) print("Pavezone cells...")
  if(length(pavezone)>0){

    cell_df_pavezone <- data.frame(ID=sapply(pavezone@polygons, function(x) methods::slot(x, name = "ID")))
    cell_df_pavezone$ID <- as.numeric(as.character(cell_df_pavezone$ID))

    #vectorised measurements
    cell_df_pavezone[,2:3] <- rgeos::gCentroid(pavezone, byid=T)@coords; colnames(cell_df_pavezone)[2:3] <- c("centroid.x", "centroid.y")
    cell_df_pavezone[,"area.px2"] <- rgeos::gArea(pavezone, byid=T)
    cell_df_pavezone[,"perimeter.px"] <- rgeos::gLength(rgeos::gBoundary(pavezone, byid=T), byid=T)

    for (j in 1:nrow(cell_df_pavezone)){

      if(verbose==TRUE & (j == round(j/20)*20| j==length(pavezone))) print(paste0("Measuring pavezone cells: ",j,"/",length(pavezone)))

      cell_j <- pavezone[j,]
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
      if("pavezone.complexity.median" %in% names(dat_im) |
         "pavezone.complexity.sd" %in% names(dat_im)
         ) cell_df_pavezone[j,"complexity"] <- rgeos::gLength(actuallines)/rgeos::gLength(straightlines)

      #GET MAX UNDULATION AMPLITUDE
      if("pavezone.undulation.amp.median" %in% names(dat_im) |
         "pavezone.undulation.amp.sd" %in% names(dat_im)
         ) cell_df_pavezone[j,"undulation.amp"] <- max_amplitude(c.smooth, junction_points)

      #GET UNDULATION FREQUENCY
      if("pavezone.undulation.freq.mean" %in% names(dat_im) |
         "pavezone.undulation.freq.sd" %in% names(dat_im)){

      crosses = rgeos::gIntersection(straightlines, actuallines)
      cell_df_pavezone[j,"undulation.freq"] <- length(crosses)/rgeos::gLength(straightlines)
      }

      #ANGLE OF ENDPOINT WALLS
      if("pavezone.endwallangle.mean" %in% names(dat_im) |
         "pavezone.endwallangle.sd" %in% names(dat_im) |
         "pavezone.endwalldiff.mean" %in% names(dat_im) |
         "pavezone.endwalldiff.sd" %in% names(dat_im)
         ) endwalls <- endwall_angles(junction_points)

      if("pavezone.endwallangle.mean" %in% names(dat_im) |
         "pavezone.endwallangle.sd" %in% names(dat_im)
         )  cell_df_pavezone[j,"endwall.angle.mean"] <- mean(c(endwalls[[3]],endwalls[[4]]))


      if("pavezone.endwalldiff.mean" %in% names(dat_im) |
         "pavezone.endwalldiff.sd" %in% names(dat_im)
         ) cell_df_pavezone[j,"endwall.angle.diff"] <- abs(endwalls[[3]]-endwalls[[4]])

      #angle of cell
      if("pavezone.angle.median" %in% names(dat_im)|
         "pavezone.angle.sd" %in% names(dat_im))  {
      coords_for_ellipse <- c.smooth@polygons[[1]]@Polygons[[1]]@coords
      ellipDirect <- conicfit::EllipseDirectFit(coords_for_ellipse)
      ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG

      angle <- ellipDirectG[5]
      if(angle>(pi/2)) angle=0-pi+angle
      cell_df_pavezone[j,"angle"] <- angle
      }
    }


    if("pavezone.area.median" %in% names(dat_im)) dat_im["pavezone.area.median"]  <- median(cell_df_pavezone$area.px2)
    if("pavezone.area.sd" %in% names(dat_im)) dat_im["pavezone.area.sd"]  <- sd(cell_df_pavezone$area.px2)
    if("pavezone.AR.median" %in% names(dat_im))  dat_im["pavezone.AR.median"]  <- median(cell_df_pavezone$AR)
    if("pavezone.AR.sd" %in% names(dat_im)) dat_im["pavezone.AR.sd"]  <- sd(cell_df_pavezone$AR)
    if("pavezone.njunctionpts.mean" %in% names(dat_im))  dat_im["pavezone.njunctionpts.mean"]  <- mean(cell_df_pavezone$n.junction.points)
    if("pavezone.njunctionpts.sd" %in% names(dat_im))  dat_im["pavezone.njunctionpts.sd"]  <- sd(cell_df_pavezone$n.junction.points)
    if("pavezone.complexity.median" %in% names(dat_im))  dat_im["pavezone.complexity.median"]  <- median(cell_df_pavezone$complexity)
    if("pavezone.complexity.sd" %in% names(dat_im))  dat_im["pavezone.complexity.sd"]  <- sd(cell_df_pavezone$complexity)
    if("pavezone.undulation.amp.median" %in% names(dat_im))  dat_im["pavezone.undulation.amp.median"]  <- median(cell_df_pavezone$undulation.amp)
    if("pavezone.undulation.amp.sd" %in% names(dat_im))  dat_im["pavezone.undulation.amp.sd"]  <- sd(cell_df_pavezone$undulation.amp)
    if("pavezone.undulation.freq.mean" %in% names(dat_im))  dat_im["pavezone.undulation.freq.mean"]  <- mean(cell_df_pavezone$undulation.freq)
    if("pavezone.undulation.freq.sd" %in% names(dat_im))  dat_im["pavezone.undulation.freq.sd"]  <- sd(cell_df_pavezone$undulation.freq)
    if("pavezone.endwallangle.mean" %in% names(dat_im))  dat_im["pavezone.endwallangle.mean"]  <- mean(cell_df_pavezone$endwall.angle.mean)
    if("pavezone.endwalldiff.mean" %in% names(dat_im))  dat_im["pavezone.endwalldiff.mean"]  <- mean(cell_df_pavezone$endwall.angle.diff)
    if("pavezone.angle.median" %in% names(dat_im))  dat_im["pavezone.angle.median"]  <- mean(cell_df_pavezone$angle)
    if("pavezone.angle.sd" %in% names(dat_im))  dat_im["pavezone.angle.sd"]  <- sd(cell_df_pavezone$angle)

  } else {
    if("pavezone.area.median" %in% names(dat_im))  dat_im["pavezone.area.median"]  <- NA
    if("pavezone.area.sd" %in% names(dat_im))  dat_im["pavezone.area.sd"]  <- NA
    if("pavezone.AR.median" %in% names(dat_im))  dat_im["pavezone.AR.median"]  <- NA
    if("pavezone.AR.sd" %in% names(dat_im))  dat_im["pavezone.AR.sd"]  <- NA
    if("pavezone.njunctionpts.mean" %in% names(dat_im))  dat_im["pavezone.njunctionpts.mean"]  <- NA
    if("pavezone.njunctionpts.sd" %in% names(dat_im))  dat_im["pavezone.njunctionpts.sd"]  <- NA
    if("pavezone.complexity.median" %in% names(dat_im))  dat_im["pavezone.complexity.median"]  <- NA
    if("pavezone.complexity.sd" %in% names(dat_im))  dat_im["pavezone.complexity.sd"]  <- NA
    if("pavezone.undulation.amp.median" %in% names(dat_im))  dat_im["pavezone.undulation.amp.median"]  <- NA
    if("pavezone.undulation.amp.sd" %in% names(dat_im))  dat_im["pavezone.undulation.amp.sd"]  <- NA
    if("pavezone.undulation.freq.mean" %in% names(dat_im))  dat_im["pavezone.undulation.freq.mean"]  <- NA
    if("pavezone.undulation.freq.sd" %in% names(dat_im))  dat_im["pavezone.undulation.freq.sd"]  <- NA
    if("pavezone.endwallangle.mean" %in% names(dat_im))  dat_im["pavezone.endwallangle.mean"]  <- NA
    if("pavezone.endwalldiff.mean" %in% names(dat_im))  dat_im["pavezone.endwalldiff.mean"]  <- NA
    if("pavezone.angle.median" %in% names(dat_im))  dat_im["pavezone.angle.median"]  <- NA
    if("pavezone.angle.sd" %in% names(dat_im))  dat_im["pavezone.angle.sd"]  <- NA
  }
}
# Stomzone (pavement cells exactly 2 cells from stomata) ####
if(exists("stomzone")) {
  if(verbose==TRUE) print("Stomatal zone cells...")
  cell_df_stomzone <- data.frame(ID=sapply(stomzone@polygons, function(x) slot(x, name = "ID")))
  cell_df_stomzone$ID <- as.numeric(as.character(cell_df_stomzone$ID))

  #vectorised measurements
  cell_df_stomzone[,2:3] <- rgeos::gCentroid(stomzone, byid=T)@coords; colnames(cell_df_stomzone)[2:3] <- c("centroid.x", "centroid.y")
  cell_df_stomzone[,"area.px2"] <- rgeos::gArea(stomzone, byid=T)
  cell_df_stomzone[,"perimeter.px"] <- rgeos::gLength(rgeos::gBoundary(stomzone, byid=T), byid=T)

  for (j in 1:nrow(cell_df_stomzone)){

    if(verbose==TRUE& (j == round(j/20)*20 | j==length(stomzone))) print(paste0("Measuring stomatal zone cells : ",j,"/",length(stomzone)))

    cell_j <- stomzone[j,]

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
    if("stomzone.undulation.amp.median" %in% names(dat_im) |
       "stomzone.undulation.amp.sd" %in% names(dat_im)
    ) cell_df_stomzone[j,"undulation.amp"] <- max_amplitude(c.smooth, junction_points)

    #where does the actual boundary cross the simplified boundary?
    if("stomzone.undulation.freq.mean" %in% names(dat_im) |
       "stomzone.undulation.freq.sd" %in% names(dat_im)){
      crosses = rgeos::gIntersection(straightlines, actuallines)
    cell_df_stomzone[j,"undulation.freq"] <- length(crosses)/rgeos::gLength(straightlines)
    }
    }

  if("stomzone.area.median" %in% names(dat_im)) dat_im["stomzone.area.median"] <- median(cell_df_stomzone$`area.px2`)
  if("stomzone.area.sd" %in% names(dat_im)) dat_im["stomzone.area.sd"] <-  sd(cell_df_stomzone$`area.px2`)
  if("stomzone.AR.median" %in% names(dat_im)) dat_im["stomzone.AR.median"] <-  median(cell_df_stomzone$AR)
  if("stomzone.complexity.median" %in% names(dat_im)) dat_im["stomzone.complexity.median"] <- median(cell_df_stomzone$complexity)
  if("stomzone.undulation.amp.median" %in% names(dat_im)) dat_im["stomzone.undulation.amp.median"] <- median(cell_df_stomzone$undulation.amp)
  if("stomzone.undulation.freq.mean" %in% names(dat_im)) dat_im["stomzone.undulation.freq.mean"] <- mean(cell_df_stomzone$undulation.freq)
} else {
  cell_df_stomzone <- NULL
}
# Polar cells (cells touching stomata - that are not subsidiary) ####
if (exists("polar")){
  if(verbose==TRUE) print("Polar cells...")
  cell_df_polar <- data.frame(ID=sapply(polar@polygons, function(x) slot(x, name = "ID")))
  cell_df_polar$ID <- as.numeric(as.character(cell_df_polar$ID))

  #vectorised measurements
  cell_df_polar[,2:3] <- rgeos::gCentroid(polar, byid=T)@coords; colnames(cell_df_polar)[2:3] <- c("centroid.x", "centroid.y")
  cell_df_polar[,"area.px2"] <- rgeos::gArea(polar, byid=T)
  cell_df_polar[,"perimeter.px"] <- rgeos::gLength(rgeos::gBoundary(polar, byid=T), byid=T)


  for (j in 1:nrow(cell_df_polar)){

    if(verbose==TRUE & (j == round(j/10)*10 | j==length(polar))) print(paste0("Measuring polar cells: ",j,"/",length(polar)))

    cell_j <- polar[j,]

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
    cell_df_polar[j,"undulation.amp"] <- max_amplitude(c.smooth, junction_points)

    #where does the actual boundary cross the simplified boundary?
    if("polar.undulation.freq.mean" %in% names(dat_im) |
       "polar.undulation.freq.sd" %in% names(dat_im)){
    crosses = rgeos::gIntersection(straightlines, actuallines)
    cell_df_polar[j,"undulation.freq"] <- length(crosses)
    }
  }

  if("polar.area.median" %in% names(dat_im)) dat_im["polar.area.median"] <- median(cell_df_polar$area.px2)
  if("polar.area.sd" %in% names(dat_im))  dat_im["polar.area.sd"] <-  sd(cell_df_polar$area.px2)
  if("polar.AR.median" %in% names(dat_im)) dat_im["polar.AR.median"] <-  median(cell_df_polar$AR)
  if("polar.complexity.median" %in% names(dat_im)) dat_im["polar.complexity.median"] <- median(cell_df_polar$complexity)
  if("polar.undulation.amp.median" %in% names(dat_im)) dat_im["polar.undulation.amp.median"] <- median(cell_df_polar$undulation.amp)
  if("polar.undulation.freq.mean" %in% names(dat_im)) dat_im["polar.undulation.freq.mean"] <- mean(cell_df_polar$undulation.freq)
} else { cell_df_polar <- NULL }
# Other cell types ####
if(!is.null(othertypes)) {
  for (i in 1:nrow(othertypes)){
othertypesi_all = cells[which(cells$value==othertypes[i,2] & cells$edge.notcounted==0),]
othertypesi = cells[which(cells$value==othertypes[i,2] & cells$edge==0),]

if(paste0(othertypes[i,1],".area.median") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".area.median")]  <- median(rgeos::gArea(othertypesi, byid = T))
if(paste0(othertypes[i,1],".area.sd") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".area.sd")] <- sd(rgeos::gArea(othertypesi, byid = T))
if(paste0(othertypes[i,1],".density") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".density")]  <- length(othertypesi_all)/(sp::bbox(cells)[1,2]*sp::bbox(cells)[2,2])
if(paste0(othertypes[i,1],".index") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".index")]   <- length(othertypesi_all)/length(cells[which(cells$edge.notcounted==0),])*100

otherdist <- cell_cell_distances(othertypesi, k=1:2)
if(paste0(othertypes[i,1],".distNN.mean") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".distNN.mean")]  <- mean(otherdist[,1])
if(paste0(othertypes[i,1],".dist2NN.mean") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".dist2NN.mean")] <- mean(otherdist[,2])
if(paste0(othertypes[i,1],".distNN.sd") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".distNN.sd")]  <- sd(otherdist[,1])
if(paste0(othertypes[i,1],".dist2NN.sd") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".dist2NN.sd")] <- sd(otherdist[,2])

if(paste0(othertypes[i,1],".spacingNN.mean") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".spacingNN.mean")] <- mean(othertypesi$dist_from_stom)
if(paste0(othertypes[i,1],".spacingNN.sd") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".spacingNN.sd")] <- sd(othertypesi$dist_from_stom)

  }
}

# Calculate SD as percent of mean ####
if(sd.as.percent.of.mean==TRUE){
  message("Stomatal angle sd not transformed")
  sds = c(
  "stom.distNN.sd"  ,
  "stom.dist2NNsd" ,
  "stom.spacingNN.sd",
  "stom.gclength.sd" ,
  "stom.AR.sd",
  "pavezone.area.sd.as.pc.median",
  "pavezone.AR.sd",
  "pavezone.njunctionpts.sd",
  "pavezone.angle.sd",
  "stomzone.area.sd",
  "polar.area.sd",
  "pavezone.complexity.sd" ,
  "pavezone.undulation.amp.sd",
  "pavezone.undulation.freq.sd")

  means = c(
    "stom.distNN.mean"  ,
    "stom.dist2NN.mean" ,
    "stom.spacingNN.mean",
    "stom.gclength.mean" ,
    "stom.AR.mean",
    "pavezone.area.median",
    "pavezone.AR.median",
    "pavezone.njunctionpts.mean",
    "pavezone.angle.median",
    "stomzone.area.median",
    "polar.area.median",
    "pavezone.complexity.median" ,
    "pavezone.undulation.amp.median",
    "pavezone.undulation.freq.mean")

  measured_sds <- sds[which(sds %in% selected_names)]
  measured_means <- means[which(sds %in% selected_names)]

  if (length(which(measured_means %in% selected_names)) != length(measured_means)) {
    warning(paste0("Some sd% unable to be calculated - means not measured. To calculate sd%, please include mean/median measurements for :
                   ",
            paste(measured_means[which(measured_means %notin% selected_names)], collapse=", ")))

    measured_sds <- measured_sds[which(measured_means %in% selected_names)]
    measured_means <- measured_means[which(measured_means %in% selected_names)]

  }

  dat_im[measured_sds] <- dat_im[measured_sds]/dat_im[measured_means]

}


# RETURN ----
return(list(image_data=dat_im,
            pavezone_individual_data=cell_df_pavezone,
            stomzone_individual_data=cell_df_stomzone,
            polar_individual_data=cell_df_polar,
            stomata_individual_data=cell_df_stom,
            polygons=cells))



}
