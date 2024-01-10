#' Extract epidermal traits from an image
#' @description This function extracts various epidermal traits from an image that has been converted to \code{sf} polygons
#'
#' @param image.ID Character. Identifier for the cells (e.g. accession number,
#' species). Defaults to NA.
#' @param cell.polygons Object of class \code{sf}. The image to be measured, converted
#' to polygons (e.g. using \code{image_to_poly()})
#' @param junction.points Either a numeric matrix or \code{sf MULTIPOINTS} object with the
#' coordinates of the cell wall junctions (e.g. as output by
#' \code{image_to_poly()})
#' @param cells.present Character. A vector containing the names of all the
#' types of cells present. Defaults to \code{c("pavement", "stomata",
#' "subsidiary")}
#' @param cell.values A vector containing the values of all the types of cells
#' present (usually numeric), see details.
#' @param NA.value A vector of length 1, with the value for NA regions of the
#' image (e.g. damaged or un-traced regions).
#' @param image.scale Numeric. The number of pixels per micron (or any other unit -
#' measurements will be returned in this unit). Defaults to NULL (measurements
#' in pixels).
#' @param image.alignment Character. The (general) axis of the image on the leaf.
#' If "horizontal" (default), stomatal north is calculated around 0. If the image
#' is vertically aligned, use "vertical" so that stomata are measured around the 90
#' degree axis. If stomata are randomly aligned and/or the leaf axis is unknown,
#' this can be set to "none".
#' @param stomatal.north Numeric.  The (precise) angle of rotation of the image on the leaf.
#' If NULL (the default), stomatal north is calculated from the angles of all stomata.
#' If stomata are randomly aligned, set to 0.
#' If rotation is known, can be set to the clockwise angle (in degrees) that image is to be rotated.
#' In some cases, it might be useful to set this to 'pavement north' (mean angle
#' of pavement cells) and traits re-extracted.
#' @param wall.width Numeric. How many pixels wide are the cell walls? Only
#' used for identifying neighbours; defaults to 3. Note: width is before scaling;
#' can be increased for a more generous buffer (i.e. more inclusive in finding
#' neighbours).
#' @param paired.guard.cells Logical. Are guard cells annotated individually (TRUE), or as a single unit (FALSE)?
#' Defaults to FALSE.
#' @param stomatal.shape Logical. Include measurements of stomatal shape
#' (`nsubcells,subsarea,gclength,AR,butterfly,symmetry`); defaults to TRUE.
#' Defaults to TRUE, see details.
#' @param stomatal.arrangement Logical. Include measurements of stomatal
#'  arrangement (`dist.between.stom.rows, row.consistency, row.wiggliness,stom.angle.sd, stom.distNN, stom.dist2NN, stom.spacingNN`);
#'  defaults to TRUE, see Details.
#' @param sd.measures Logical. Return standard deviation of variables? This
#' argument controls the measurement of variability of traits within an
#' image; it is only recommended to use this when a) the sample size is
#' high, or b) there is strong indication that these traits are useful.
#'  Defaults to FALSE.
#' @param specific.inclusions Character. Cell measurements to include;
#'  see Details for options.
#' @param specific.exclusions Character. Cell measurements to exclude;
#'  can be passed as key words (e.g. "AR", "undulation", "pavezone").
#' @param sd.as.percent.of.mean Logical. Should variability measures be
#' presented as percentage of the mean? Defaults to FALSE.
#' @param verbose Logical. Should progress output be printed in the console?
#' Defaults to TRUE.
#'
#' @return A list of length three, containing a vector with whole-image values,
#' a data.frame of pavement cell values and a data.frame of stomatal values.
#' Variables are coded as 'celltype.trait.statistic' (e.g. `pavezone.area.sd`).
#'
#' @details This function contains the script for all measurements. Each type
#' of cell (e.g. pavement, stomate) should have a name and corresponding value
#' given in \code{cells.present} and \code{cell.values}, respectively. If
#' there are unnamed values in the image, these will be measured as
#' 'othertype1', 'othertype2' and so forth. Only limited numbers of
#' measurements will be taken for these cell types (index, density, area,
#' spacing and distance to nearest neighbour). If only a subset of variables is
#' desired, this can can be set using the arguments \code{ stomatal.shape,
#' stomatal.arrangement, pavement.cells, sd.measures, specific.inclusions,
#' specific.exclusions}. Pavement cells are designated as `polar` (if touching
#' the guard cells; we acknowledge that this is not the most universal term for
#' possible subsidiary cells but lack an alternative), `stomzone` (for cells that
#' are exactly two neightbours from the guard cells and are likely to be distorted)
#' and `pavezone` (cells that are >2 neighbours from the guard cells and are
#' thus unlikely to be distorted).
#'
#' For a full list of variables and descriptions of what they measure,
#' see the data `trait_key`, the corresponding paper (for graphical descriptions),
#' and the vignettes.
#'
#' @importFrom graphics par
#' @import methods
#' @import stats cli
#'
#' @export


extract_epidermal_traits <- function( image.ID=NA, cell.polygons,
                                      junction.points,
                                      cells.present = c("pavement",
                                                         "stomata",
                                                         "subsidiary"),
                                      cell.values,
                                      paired.guard.cells=FALSE,
                                      NA.value=NULL,
                                      image.scale=NULL,
                                      image.alignment="horizontal",
                                      stomatal.north=NULL,
                                      wall.width=3,
                                      stomatal.shape = TRUE,
                                      stomatal.arrangement = TRUE,
                                      sd.measures = TRUE,
                                      specific.inclusions = NULL,
                                      specific.exclusions = NULL,
                                      sd.as.percent.of.mean = FALSE,
                                      verbose=TRUE) {

  . <- rows <- angle.corrected <- NULL

# SETUP ----

  if (!is.na(image.ID)) {
    if(verbose==TRUE) cli_h1(paste0("Extracting epidermal traits from ", image.ID))
  } else {
    if(verbose==TRUE) cli_h1(paste0("Extracting epidermal traits from image"))
  }
#Pre-empt some common input errors ####
  if (length(cells.present) != length(cell.values)) cli_abort("Need to supply the same number of cell values as cell types.")
  if (length(unique(cell.values))!=length(cell.values)) cli_abort("Some cell values repeated.")
  if (length(unique(cells.present))!=length(cells.present)) cli_abort("Some cell types repeated.")
  if(!is.null(NA.value)) if(NA.value %in% cell.values) cli_abort("NA.value muct be different from defined cell type values.")

  if(!is.null(image.scale)){
    sf::st_geometry(cell.polygons) <- sf::st_geometry(cell.polygons)/image.scale
    junction.points <- junction.points/image.scale
    buffervalue <- wall.width*1.5/image.scale
  } else {
    cli_alert_info("No scale set; all measurements are in pixels")
    buffervalue=wall.width*1.5
  }

  if(!is.null(NA.value)) cells=cell.polygons[which(cell.polygons$value!=NA.value),] else cells=cell.polygons


sf::st_agr(cells) = "constant"


#Identify edge cells ####

edges <- find_edge_cells(cell.polygons,cells)

cells$edge <- 0
cells$edge[edges[[1]]] <- 1
cells$edge.notcounted <- 0
cells$edge.notcounted[edges[[2]]] <- 1


# Clean up junction points ####
wjp <- get_junctions(junction.points, buffervalue*0.5)

# Find stomatal north, rotate cells accordingly ####

if (image.alignment=="vertical"){
  center <- suppressWarnings(sf::st_centroid(sf::st_union(cells)))
  sf::st_geometry(cells) <- tran(sf::st_geometry(cells), 90, center)
  wjp <- tran(sf::st_geometry(wjp), 90, center)
}

if ("stomata" %in% cells.present){
  if("subsidiary" %in% cells.present){
    subs.val <- cell.values[which(cells.present=="subsidiary")]
  } else {subs.val <- NULL}

stom_df <- identify_stomatal_complexes(cells=cells,
                                       paired=paired.guard.cells,
                                       buffervalue=buffervalue*2,
                                       stom.val=cell.values[which(cells.present=="stomata")],
                                       subs.val=subs.val)



if(is.null(stomatal.north)){
stomatal_north <- mean(stom_df$angle, na.rm=TRUE)
} else {
  stomatal_north <- stomatal.north
}

# align to stomatal north
center <- suppressWarnings(sf::st_centroid(sf::st_union(cells)))

sf::st_geometry(cells) <- tran(sf::st_geometry(cells), stomatal_north, center)
sf::st_geometry(stom_df) <- tran(sf::st_geometry(stom_df), stomatal_north, center)

wjp <- tran(sf::st_geometry(wjp), stomatal_north, center)
}

# Split up cells into cell types ####

#set up stomata (note: stomata re-extracted from rotated polygons)
  if ("stomata" %in% cells.present){
    stom.val <- cell.values[which(cells.present=="stomata")]

    cells <- cells[cells$value!=stom.val,]
    cells <- rbind(cells, stom_df[,colnames(cells)])
    stom <- cells[which(cells$value==stom.val),]
    stom_keep  <- stom[which(stom$edge==0),]

    # Sort out pave zone/stom zone/polar cells ####
      paths <- cell_graph_shortest_paths(cells, buffervalue)
      cells$dist_from_stom <- shortest_path_to_cell_type(cells, paths, stom.val)

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
    if("subsidiary" %in% cells.present){
      pavezone <- pave[which(pave$edge==0 & pave$dist_from_stom>2),]
      stomzone <- pave[which(pave$edge==0 & pave$dist_from_stom==2),]
      polar <- pave[which(pave$edge==0 & pave$dist_from_stom==1),]
    } else {
      pavezone <- pave[which(pave$edge==0 & pave$dist_from_stom>1),]
      stomzone <- pave[which(pave$edge==0 & pave$dist_from_stom==1),]
      polar <- NULL
    }


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
  cli_alert_info("Not all values in cells have been allocated a type. These cells will be named 'othertype1', 'othertype2', etc.")
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

trait_key <- epidermalmorph::trait_key

stomatal_basic_names <- trait_key$trait[which(trait_key$category=="stomata basic")]

stom_arrangement_names <- trait_key$trait[
  which(
    trait_key$category=="stomatal arrangement" &
    trait_key$sd == 0
    )
  ]

stom_shape_names <- trait_key$trait[
  which(
    trait_key$category=="stomatal shape" &
      trait_key$sd == 0
  )
]

pavement_names <- trait_key$trait[
  which(
    trait_key$category=="pavement" &
      trait_key$sd == 0
  )
]

other_cell_type_names <- trait_key$trait[
  which(
    trait_key$category=="other" &
      trait_key$sd == 0
  )
]

stom_arr_sd_names <- trait_key$trait[
  which(
    trait_key$category=="stomatal arrangement" &
      trait_key$sd == 1
  )
]
stom_shape_sd_names <- trait_key$trait[
  which(
    trait_key$category=="stomatal shape" &
      trait_key$sd == 1
  )
]
pavement_sd_names <- trait_key$trait[
  which(
    trait_key$category=="pavement" &
      trait_key$sd == 1
  )
]
other_cell_type_sd_names <- trait_key$trait[
  which(
    trait_key$category=="other" &
      trait_key$sd == 1
  )
]

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

if(!is.null(specific.exclusions)) {
  exclusions <- NULL
  for(i in 1:length(specific.exclusions)){
    sei <- specific.exclusions[i]
    if(!sei %in% selected_names) {
      sei <- trait_key$trait[which(grepl(sei, trait_key$trait)==TRUE)]
      exclusions <- c(exclusions, sei)
    }
  }
  if(length(exclusions)==0) {
    cli::cli_alert_warning("Could not match selected.exclusions input to traits; no traits have been excluded")
  } else {
  exclusions <- unique(exclusions[order(exclusions)])
    cli::cli_alert_info(paste0("The following traits match 'specific.exclusions' and will not be measured: \n- ",paste(exclusions, collapse="\n- ")))
    selected_names <- selected_names[which(selected_names %notin% exclusions)]
  }
}


dat_im <- rep(NA,times=length(selected_names))
names(dat_im) <- selected_names

if(!is.null(stomatal.north)) if(stomatal.north == 0 & any(grepl("endwall", selected_names))==TRUE) cli::cli_alert_warning("Endwalls may be unreliable if image is not aligned with leaf axis - use specific.exclusions = 'endwall' to disable.")

####


# MEASUREMENTS ----

# Stomatal arrangement measurements ####

#stomatal north
if("rotation" %in% names(dat_im)) dat_im["rotation"] <- stomatal_north

#number of each cell type (NOT ON THE EDGE)
n.cell.types <- table(cells$value[which(cells$edge==0)])
for (i in 1:nrow(alltypes)){
  dat_im[paste0("n.",alltypes[i,1])] <- as.numeric(n.cell.types[which(names(n.cell.types)==alltypes[i,2])])
}

if("stomatal.index" %in% names(dat_im)) dat_im["stomatal.index"] <- nrow(cells[which(cells$edge.notcounted==0 &
                                                                                            cells$value==stom.val),]
                                                                         )/nrow(cells[which(cells$edge.notcounted==0),])*100
if("stomatal.density" %in% names(dat_im)) dat_im["stomatal.density"] <- nrow(cells[which(cells$edge.notcounted==0 &
                                                                                                        cells$value==stom.val),]
                                                                                     )/(sf::st_bbox(cells)[3]*sf::st_bbox(cells)[4])

if("dist.between.stom.rows" %in% names(dat_im) |
   "row.wiggliness" %in% names(dat_im) |
   "row.consistency" %in% names(dat_im) ) {
  stomatal_rows <- number_rows_stomata(cells[which(cells$edge==0),])
  if("dist.between.stom.rows" %in% names(dat_im)) dat_im["dist.between.stom.rows"] <- stomatal_rows[[5]]
  if("row.wiggliness" %in% names(dat_im)) dat_im["row.wiggliness"] <- mean(stomatal_rows[[3]])
  if("row.consistency" %in% names(dat_im)) dat_im["row.structure"] <- sd(stomatal_rows[[4]])
}

# mean distance to kth nearest neighbours

if("stom.distNN.mean" %in% names(dat_im) |
   "stom.distNN.mean" %in% names(dat_im) |
   "stom.distNN.mean" %in% names(dat_im) |
   "stom.dist2NN.sd" %in% names(dat_im)) {

  #updated to deal with paired stom
  stomdist <- cell_cell_distances(stom_df, k=1:2)

  if("stom.distNN.mean" %in% names(dat_im)) dat_im["stom.distNN.mean"] <- mean(stomdist[,1])
  if("stom.distNN.sd" %in% names(dat_im)) dat_im["stom.distNN.sd"]<- sd(stomdist[,1])

  if("stom.dist2NN.mean" %in% names(dat_im)) dat_im["stom.dist2NN.mean"] <- mean(stomdist[,2])
  if("stom.dist2NN.sd" %in% names(dat_im)) dat_im["stom.dist2NN.sd"]<- sd(stomdist[,2])

}

if("stom.spacingNN.mean" %in% names(dat_im)) dat_im["stom.spacingNN.mean"] <- mean(cells$dist_from_stom[which(cells$value==stom.val)],na.rm=TRUE)
if("stom.spacingNN.sd" %in% names(dat_im)) dat_im["stom.spacingNN.sd"] <- sd(cells$dist_from_stom[which(cells$value==stom.val)], na.rm=TRUE)

if("stom.nsubcells.mean" %in% names(dat_im)) {
  if("subsidiary" %in% cells.present) {
    dat_im["stom.nsubcells.mean"] <- nrow(cells[which(cells$edge.notcounted ==0 & cells$value==subs.val),])/nrow(cells[which(cells$edge.notcounted==0 & cells$value==stom.val),])
  } else {
    dat_im["stom.nsubcells.mean"] <- 0
  }
  }

if("stom.subsarea.mean" %in% names(dat_im)) {
  if("subsidiary" %in% cells.present) {
    dat_im["stom.subsarea.mean"] <-as.numeric(sum(sf::st_area(subs_keep)))/nrow(stom_keep)
  } else {
    dat_im["stom.subsarea.mean"] <- 0
  }
}


# Stomatal shape measurements ####
if("stomata" %in% cells.present){

  cell_df_stom <- data.frame(ID=as.numeric(row.names(stom_keep)))

  suppressWarnings(cell_df_stom[,2:3] <- sf::st_coordinates(sf::st_centroid(stom_keep)))
  colnames(cell_df_stom)[2:3] <- c("centroid.x", "centroid.y")


  cell_df_stom$angle <-  stom_df %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(rows=rownames(.)) %>%
    dplyr::filter(rows %in% cell_df_stom$ID) %>%
    dplyr::pull(angle.corrected)

  SUBS <- stom_df$SUBS

  if(verbose==TRUE) cli_progress_bar("Calculating stomatal features", total = nrow(cell_df_stom))
  for (j in 1:nrow(cell_df_stom)){

    if(verbose==TRUE)  cli_progress_update()
    stom_j <- stom_keep[j,]

    if("subsidiary" %in% cells.present) {
    list.element =  which(rownames(stom_df) == as.character(cell_df_stom$ID)[j])
    subs_j <- subs[as.character(SUBS[[list.element]]),]
      stom_comp_j <- rbind(stom_j, subs_j[,-5])
    } else {
      stom_comp_j <- stom_j
    }

    if(length(find_edge_cells(cell.polygons, stom_comp_j, rotated=stomatal_north)$edge)==0){

      if("subsidiary" %in% cells.present) {
      if(nrow(subs_j)>2){
        bsubs <- sf::st_buffer(subs_j, dist=buffervalue/5)
        subs_key <-  as.matrix(sf::st_intersects(bsubs))
        subs1 <- which(subs_key[1,]==TRUE)
        subs2 <- which(1:nrow(subs_j) %notin% subs1)
        if(length(subs2)==0){
          subskey2 <- subs_key
          for (m in 1:(nrow(subskey2)-1)){
            for(n in (m+1):nrow(subskey2)){
              if(subskey2[m,n]==T)  subskey2[m,n] <- sf::st_area(sf::st_intersection(bsubs[m,], bsubs[n,]))
              if(subskey2[m,n]==F)  subskey2[m,n] <- 0
            }
          }
          subs1 <- as.numeric(which(subskey2==max(subskey2), arr.ind = T))
          subs2 <- as.numeric(which(1:length(subs_j) %notin% subs1))
        }
      }
      }


      cell_df_stom[j,"gclength"] <- sf::st_bbox(stom_j)[3]-sf::st_bbox(stom_j)[1]

      cell_df_stom[j,"length"] <- sf::st_bbox(stom_comp_j)[3]-sf::st_bbox(stom_comp_j)[1]
      cell_df_stom[j,"width"] <- sf::st_bbox(stom_comp_j)[4]-sf::st_bbox(stom_comp_j)[2]
      if("stom.butterfly.mean" %in% names(dat_im)) cell_df_stom[j,"butterfly"] <- cell_df_stom[j,"gclength"]/cell_df_stom[j,"length"]
      if("subsidiary" %in% cells.present & "stom.symmetry.mean" %in% names(dat_im)){

      if(nrow(subs_j)==2) cell_df_stom[j,"symmetry"] <- sf::st_area(subs_j[1,])/sf::st_area(subs_j[2,])
      if(nrow(subs_j)>2) cell_df_stom[j,"symmetry"] <- sf::st_area(sf::st_union(subs_j[subs1,]))/sf::st_area(sf::st_union(subs_j[subs2,]))
      if(nrow(subs_j)<2) {cell_df_stom[j,"symmetry"] <- NA

      } else {
        if(cell_df_stom[j,"symmetry"] >1) cell_df_stom[j,"symmetry"]<- 1/cell_df_stom[j,"symmetry"]}
      }

    } else {
      cell_df_stom[j,] <- NA
    }

  }


  if("stom.gclength.mean" %in% names(dat_im)) dat_im["stom.gclength.mean"] <-  mean(cell_df_stom$gclength, na.rm=T)
  if("stom.gclength.sd" %in% names(dat_im)) dat_im["stom.gclength.sd"] <- sd(cell_df_stom$gclength, na.rm=T)
  if("stom.AR.mean" %in% names(dat_im)) dat_im["stom.AR.mean"] <- mean(cell_df_stom$length/cell_df_stom$width, na.rm=T)
  if("stom.AR.sd" %in% names(dat_im)) dat_im["stom.AR.sd"] <- sd(cell_df_stom$length/cell_df_stom$width, na.rm=T)
  if("stom.butterfly.mean" %in% names(dat_im)) dat_im["stom.butterfly.mean"] <- mean(cell_df_stom$butterfly, na.rm=T)
  if("stom.symmetry.mean" %in% names(dat_im)) {
    if("subsidiary" %in% cells.present){
      dat_im["stom.symmetry.mean"] <- mean(cell_df_stom$symmetry, na.rm=T)
    } else {dat_im["stom.symmetry.mean"] <- NA }
  }

  if("stom.angle.sd" %in% names(dat_im)) dat_im["stom.angle.sd"] <- sd(stom_df$angle.corrected, na.rm=TRUE)

} else cell_df_stom <- NULL

if(verbose==TRUE) cli_alert_success("Stomatal features measured.")

# Pavezone (>2 cells away from stomata) ####
if("pavement" %in% cells.present) {


  if(nrow(pavezone)>0){

    cell_df_pavezone <- data.frame(ID=as.numeric(row.names(pavezone)))
    if(verbose==TRUE)  cli_progress_bar("Pavezone cells", total = nrow(cell_df_pavezone))

    #vectorised measurements
    cell_df_pavezone[,2:3] <- suppressWarnings(sf::st_coordinates(sf::st_centroid(pavezone)))
    colnames(cell_df_pavezone)[2:3] <- c("centroid.x", "centroid.y")
    cell_df_pavezone[,"area"] <- sf::st_area(pavezone)
    cell_df_pavezone[,"perimeter"] <- sf::st_length(sf::st_boundary(pavezone))


    for (j in 1:nrow(cell_df_pavezone)){
      if(verbose==TRUE)  cli_progress_update()


      cell_j <- pavezone[j,]
      cell_j_bbox <- sf::st_bbox(cell_j)
      cell_df_pavezone[j,"AR"] <- (cell_j_bbox[3]-cell_j_bbox[1])/(cell_j_bbox[4]-cell_j_bbox[2])
      c.smooth<- smoothr::smooth(cell_j, method = "ksmooth", smoothness = 3)



      #GET THE JUNCTION POINTS THAT DEFINE THE SIMPLIFIED CELL
      junction_points <- cell_simplify(cell=c.smooth, cell.junctions=wjp, snap.tolerance = buffervalue*2)

      #make sure no double ups
      if(nrow(unique(junction_points))<(nrow(junction_points)-1)){
        n <- nrow(junction_points)
        junction_points <- rbind(junction_points[1,], unique(junction_points[2:(n-1),]), junction_points[n,])
      }
      cell_df_pavezone[j,"n.junction.points"] <- nrow(junction_points)-1

      straightpoly = sf::st_polygon(list(junction_points))
      straightlines = sf::st_linestring(junction_points)
      actuallines = suppressWarnings(sf::st_cast(c.smooth, 'LINESTRING'))
      if(nrow(actuallines)>1) actuallines <- actuallines[which(sf::st_length(actuallines) == max(sf::st_length(actuallines))),]

      #GET difference between straight perimeter and actual perimeter
      if("pavezone.complexity.median" %in% names(dat_im) |
         "pavezone.complexity.sd" %in% names(dat_im)
      ) cell_df_pavezone[j,"complexity"] <- sf::st_length(actuallines)/sf::st_length(straightlines)


      #GET MAX UNDULATION AMPLITUDE
      if("pavezone.undulation.amp.median" %in% names(dat_im) |
         "pavezone.undulation.amp.sd" %in% names(dat_im)
         ) cell_df_pavezone[j,"undulation.amp"] <- max_amplitude(c.smooth, junction_points)

      #GET UNDULATION FREQUENCY
      if("pavezone.undulation.freq.mean" %in% names(dat_im) |
         "pavezone.undulation.freq.sd" %in% names(dat_im)){

      crosses = sf::st_intersection(straightlines, actuallines)
      cell_df_pavezone[j,"undulation.freq"] <- length(crosses)/sf::st_length(straightlines)
      }

      #ANGLE OF ENDPOINT WALLS
      if("pavezone.endwallangle.mean" %in% names(dat_im) |
         "pavezone.endwallangle.sd" %in% names(dat_im) |
         "pavezone.endwalldiff.mean" %in% names(dat_im) |
         "pavezone.endwalldiff.sd" %in% names(dat_im)
         ) {
       endwalls <- endwall_angles(junction_points)
      }

      if("pavezone.endwallangle.mean" %in% names(dat_im) |
         "pavezone.endwallangle.sd" %in% names(dat_im)
         )  cell_df_pavezone[j,"endwall.angle.mean"] <- mean(c(endwalls[[3]],endwalls[[4]]))


      if("pavezone.endwalldiff.mean" %in% names(dat_im) |
         "pavezone.endwalldiff.sd" %in% names(dat_im)
         ) cell_df_pavezone[j,"endwall.angle.diff"] <- abs(endwalls[[3]]-endwalls[[4]])

      #angle of cell
      if("pavezone.angle.median" %in% names(dat_im)|
         "pavezone.angle.sd" %in% names(dat_im))  {
      coords_for_ellipse <- sf::st_coordinates(c.smooth)
      ellipDirect <- conicfit::EllipseDirectFit(coords_for_ellipse)
      ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG

      angle <- ellipDirectG[5]
      if(angle>(pi/2)) angle=0-pi+angle
      cell_df_pavezone[j,"angle"] <- angle
      }
    }


    if("pavezone.area.median" %in% names(dat_im)) dat_im["pavezone.area.median"]  <- median(cell_df_pavezone$area)
    if("pavezone.area.sd" %in% names(dat_im)) dat_im["pavezone.area.sd"]  <- sd(cell_df_pavezone$area)
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

  cell_df_stomzone <- data.frame(ID=as.numeric(row.names(stomzone)))

  #vectorised measurements
  cell_df_stomzone[,2:3] <- suppressWarnings(sf::st_coordinates(sf::st_centroid(stomzone)))
  colnames(cell_df_stomzone)[2:3] <- c("centroid.x", "centroid.y")
  cell_df_stomzone[,"area"] <- sf::st_area(stomzone)
  cell_df_stomzone[,"perimeter"] <- sf::st_length(sf::st_boundary(stomzone))

  if(verbose==TRUE) cli_progress_bar("Stomatal zone cells", total = nrow(cell_df_stomzone))
 for (j in 1:nrow(cell_df_stomzone)){
   if(verbose==TRUE) cli_progress_update()

    cell_j <- stomzone[j,]
    cell_j_bbox <- sf::st_bbox(cell_j)
    cell_df_stomzone[j,"AR"] <- (cell_j_bbox[3]-cell_j_bbox[1])/(cell_j_bbox[4]-cell_j_bbox[2])
    c.smooth<- smoothr::smooth(cell_j, method = "ksmooth", smoothness = 3)

    #GET THE JUNCTION POINTS THAT DEFINE THE SIMPLIFIED CELL
    junction_points <- cell_simplify(cell=c.smooth, cell.junctions=wjp,snap.tolerance = buffervalue*2)

    #make sure no double ups
    if(nrow(unique(junction_points))<(nrow(junction_points)-1)){
      n <- nrow(junction_points)
      junction_points <- rbind(junction_points[1,], unique(junction_points[2:(n-1),]), junction_points[n,])
    }
    cell_df_stomzone[j,"n.junction.points"] <- nrow(junction_points)-1

    straightpoly = sf::st_polygon(list(junction_points))
    straightlines = sf::st_linestring(junction_points)
    actuallines = suppressWarnings(sf::st_cast(c.smooth, 'LINESTRING'))
    if(nrow(actuallines)>1) actuallines <- actuallines[which.max(sf::st_length(actuallines)),]

    #GET difference between straight perimeter and actual perimeter
    if("stomzone.complexity.median" %in% names(dat_im) |
       "stomzone.complexity.sd" %in% names(dat_im)
    ) cell_df_stomzone[j,"complexity"] <- sf::st_length(actuallines)/sf::st_length(straightlines)

    #GET MAX UNDULATION AMPLITUDE
    if("stomzone.undulation.amp.median" %in% names(dat_im) |
       "stomzone.undulation.amp.sd" %in% names(dat_im)
    ) cell_df_stomzone[j,"undulation.amp"] <- max_amplitude(c.smooth, junction_points)

    #GET UNDULATION FREQUENCY
    if("stomzone.undulation.freq.mean" %in% names(dat_im) |
       "stomzone.undulation.freq.sd" %in% names(dat_im)){

      crosses = sf::st_intersection(straightlines, actuallines)
      cell_df_stomzone[j,"undulation.freq"] <- length(crosses)/sf::st_length(straightlines)
    }
  }

  if("stomzone.area.median" %in% names(dat_im)) dat_im["stomzone.area.median"] <- median(cell_df_stomzone$`area`)
  if("stomzone.area.sd" %in% names(dat_im)) dat_im["stomzone.area.sd"] <-  sd(cell_df_stomzone$`area`)
  if("stomzone.AR.median" %in% names(dat_im)) dat_im["stomzone.AR.median"] <-  median(cell_df_stomzone$AR)
  if("stomzone.complexity.median" %in% names(dat_im)) dat_im["stomzone.complexity.median"] <- median(cell_df_stomzone$complexity)
  if("stomzone.undulation.amp.median" %in% names(dat_im)) dat_im["stomzone.undulation.amp.median"] <- median(cell_df_stomzone$undulation.amp)
  if("stomzone.undulation.freq.mean" %in% names(dat_im)) dat_im["stomzone.undulation.freq.mean"] <- mean(cell_df_stomzone$undulation.freq)
} else {
  cell_df_stomzone <- NULL
}
# Polar cells (cells touching stomata - that are not subsidiary) ####
if(exists("polar")& !is.null(polar)){


  cell_df_polar <- data.frame(ID=as.numeric(row.names(polar)))

  #vectorised measurements
  cell_df_polar[,2:3] <- suppressWarnings(sf::st_coordinates(sf::st_centroid(polar)))
  colnames(cell_df_polar)[2:3] <- c("centroid.x", "centroid.y")
  cell_df_polar[,"area"] <- sf::st_area(polar)
  cell_df_polar[,"perimeter"] <- sf::st_length(sf::st_boundary(polar))

  if(verbose==TRUE)  cli_progress_bar("Polar cells", total = nrow(cell_df_polar))
  for (j in 1:nrow(cell_df_polar)){
    if(verbose==TRUE) cli_progress_update()

    cell_j_bbox <- sf::st_bbox(cell_j)
    cell_df_polar[j,"AR"] <- (cell_j_bbox[3]-cell_j_bbox[1])/(cell_j_bbox[4]-cell_j_bbox[2])
    c.smooth<- smoothr::smooth(cell_j, method = "ksmooth", smoothness = 3)


    #GET THE JUNCTION POINTS THAT DEFINE THE SIMPLIFIED CELL
    junction_points <- cell_simplify(cell=c.smooth, cell.junctions=wjp,snap.tolerance = buffervalue*2)

    #make sure no double ups
    if(nrow(unique(junction_points))<(nrow(junction_points)-1)){
      n <- nrow(junction_points)
      junction_points <- rbind(junction_points[1,], unique(junction_points[2:(n-1),]), junction_points[n,])
    }
    cell_df_polar[j,"n.junction.points"] <- nrow(junction_points)-1

    straightpoly = sf::st_polygon(list(junction_points))
    straightlines = sf::st_linestring(junction_points)
    actuallines = suppressWarnings(sf::st_cast(c.smooth, 'LINESTRING'))
    if(nrow(actuallines)>1) actuallines <- actuallines[which(sf::st_length(actuallines) == max(sf::st_length(actuallines))),]

    #GET difference between straight perimeter and actual perimeter
    if("polar.complexity.median" %in% names(dat_im) |
       "polar.complexity.sd" %in% names(dat_im)
    ) cell_df_polar[j,"complexity"] <- sf::st_length(actuallines)/sf::st_length(straightlines)

    #GET MAX UNDULATION AMPLITUDE
    if("polar.undulation.amp.median" %in% names(dat_im) |
       "polar.undulation.amp.sd" %in% names(dat_im)
    ) cell_df_polar[j,"undulation.amp"] <- max_amplitude(c.smooth, junction_points)

    #GET UNDULATION FREQUENCY
    if("polar.undulation.freq.mean" %in% names(dat_im) |
       "polar.undulation.freq.sd" %in% names(dat_im)){

      crosses = sf::st_intersection(straightlines, actuallines)
      cell_df_polar[j,"undulation.freq"] <- length(crosses)/sf::st_length(straightlines)
    }
  }

  if("polar.area.median" %in% names(dat_im)) dat_im["polar.area.median"] <- median(cell_df_polar$area)
  if("polar.area.sd" %in% names(dat_im))  dat_im["polar.area.sd"] <-  sd(cell_df_polar$area)
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

paths <- cell_graph_shortest_paths(cells, buffervalue)[which(cells$value==othertypes[i,2] & cells$edge==0),which(cells$value==othertypes[i,2] & cells$edge==0)]
grouping <- other_cluster_size(cells,othertypes[i,2], buffervalue)

spacing <- apply(paths,2,kthsmallest,k=grouping+1)

otherdist <- cell_cell_distances(othertypesi, k=1:grouping)[,grouping]

if(paste0(othertypes[i,1],".area.median") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".area.median")]  <- median(sf::st_area(othertypesi))
if(paste0(othertypes[i,1],".area.sd") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".area.sd")] <- sd(sf::st_area(othertypesi))
if(paste0(othertypes[i,1],".density") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".density")]  <- length(othertypesi_all)/(sf::st_bbox(cells)[3]*sf::st_bbox(cells)[4])
if(paste0(othertypes[i,1],".index") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".index")]   <- length(othertypesi_all[which(othertypesi_all$edge.notcounted==0),])/length(cells[which(cells$edge.notcounted==0),])*100
if(paste0(othertypes[i,1],".grouping") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".grouping")]   <- grouping

if(paste0(othertypes[i,1],".distNN.mean") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".distNN.mean")]  <- mean(otherdist)
if(paste0(othertypes[i,1],".distNN.sd") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".distNN.sd")]  <- sd(otherdist)

if(paste0(othertypes[i,1],".spacingNN.mean") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".spacingNN.mean")] <- mean(spacing)
if(paste0(othertypes[i,1],".spacingNN.sd") %in% names(dat_im)) dat_im[paste0(othertypes[i,1],".spacingNN.sd")] <- sd(spacing)

  }
}


if(verbose==TRUE) cli_alert_success("Pavement cell features measured.")
# Calculate SD as percent of mean ####
if(sd.as.percent.of.mean==TRUE){
  cli_inform("Stomatal angle sd not transformed")
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
    cli_alert_warning(paste0("Some sd% unable to be calculated - means not measured. To calculate sd%, please include mean/median measurements for :
                   ",
            paste(measured_means[which(measured_means %notin% selected_names)], collapse=", ")))

    measured_sds <- measured_sds[which(measured_means %in% selected_names)]
    measured_means <- measured_means[which(measured_means %in% selected_names)]

  }

  dat_im[measured_sds] <- dat_im[measured_sds]/dat_im[measured_means]

}

dat_im <- data.frame(t(dat_im))
dat_im$image.ID <- image.ID
cells$ID <- as.numeric(row.names(cells))

# RETURN ----
if(!exists("cell_df_pavezone")) cell_df_pavezone <- NULL
if(!exists("cell_df_stomzone")) cell_df_stomzone <- NULL
if(!exists("cell_df_polar")) cell_df_polar <- NULL
if(!exists("cell_df_stom")) cell_df_stom <- NULL

return(list(image_data=dat_im,
            pavezone_individual_data=cell_df_pavezone,
            stomzone_individual_data=cell_df_stomzone,
            polar_individual_data=cell_df_polar,
            stomata_individual_data=cell_df_stom,
            polygons=cells))



}
