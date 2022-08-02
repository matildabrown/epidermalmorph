#' Identify which cells belong to particular stomatal complexes
#'
#' @param cells Object of class \code{sf}. The image to be measured, converted
#' to polygons (e.g. using \code{image_to_poly()})
#' @param paired Logical. Are guard cells annotated separately (TRUE), or as one
#' cell (FALSE? Defaults to FALSE.
#' @param buffervalue Amount to buffer cells by. Should be similar to wall width (which varies with image scaling).
#' @param stom.val Numeric. Value of stomata in inpute image.
#' @param subs.val Numeric. Value of subsidiary cells in image (optional).
#' Defaults to NULL.
#'
#' @return sf object with the IDs of all cells belonging to the stomatal complex.
#' @export
#'
#'
identify_stomatal_complexes <- function(cells, paired=F, buffervalue, stom.val, subs.val=NULL){

value.1 <- edge.1 <- edge.notcounted.1 <- angle.corrected <- ID <- NULL

  #get stomata from cells
   stom_keep <- cells[which(cells$value==stom.val),]


#individual stomata
  if (paired==F){
    #set up dataframe
    stom_key <- stom_keep
    stom_key$GC1 <- rownames(stom_keep)

    #get angle, using ellipse method
    stom_key$angle <- NA
    for (i in 1:nrow(stom_key)){
      ####ellipse bit
      c.smooth <-   smoothr::smooth(stom_key[i,], method = "ksmooth", smoothness = 3)
      coords_for_ellipse <- sf::st_coordinates(c.smooth)
      ellipDirect <- conicfit::EllipseDirectFit(coords_for_ellipse)
      ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG

      angle <- ellipDirectG[5]
      #correct so angles fall +/- 0
      if(angle>(pi/2)) angle=0-pi+angle

      stom_key[i,"angle"] <- pracma::rad2deg(angle)

    }

    #else for paired guard cell polygons
  } else {
    #identify pairs of guard cells making up a stomate using intersection of buffers
    bstom  <- sf::st_buffer(stom_keep, dist=buffervalue/2)
    pairs <- unique((sf::st_intersects(bstom,bstom)))

    #check for pairs
    if(max(lengths(pairs))!=2) cli::cli_abort("Cannot identify corresponding guard cell pairs. If stomata are closely packed, consider merging guard cell pairs during image segmentation or setting a smaller `wall.width` (see Details; ?identify_stomatal_complexes)")
    if(abs(length(pairs)-nrow(stom_keep)/2)>3) cli::cli_abort("Cannot identify corresponding guard cell pairs. Consider merging guard cell pairs during image segmentation or setting a smaller `wall.width` (see Details; ?identify_stomatal_complexes)")


    complexes <- data.frame(do.call(rbind, pairs))

    #get ids
    complexes$GC1 <- as.numeric(row.names(stom_keep))[complexes$X1]
    complexes$GC2 <- as.numeric(row.names(stom_keep))[complexes$X2]
    newstom = list()
    for(i in 1:nrow(complexes)) {
      suppressWarnings(newstom[[i]] <-  sf::st_union(stom_keep[complexes$X1[i],],
                                                 stom_keep[complexes$X2[i], ])
                    )
      newstom[[i]]$edge <- max(c(newstom[[i]]$edge, newstom[[i]]$edge.1))
      newstom[[i]]$edge.notcounted <- max(c(newstom[[i]]$edge.notcounted, newstom[[i]]$edge.notcounted.1))
      newstom[[i]] <- newstom[[i]] %>% dplyr::select(-c(value.1, edge.1, edge.notcounted.1))

      if(newstom[[i]]$edge==1) {
        newstom[[i]]$angle <- NA
      } else{
      #get angle, using GC intercection method
      suppressWarnings(pts <- sf::st_intersection(sf::st_boundary(bstom[complexes$X1[i],]),
                          sf::st_boundary(bstom[complexes$X2[i],])) %>% sf::st_coordinates()
                       )
      if(nrow(pts)>2){
        d <- as.matrix(dist(pts))
        d[upper.tri(d, diag = T)] <- NA
        pts <- pts[which(d==max(d, na.rm = T), arr.ind = T),]
      }
      angle <- atan((pts[2,2]-pts[1,2])/(pts[2,1]-pts[1,1]))
      if(angle>(pi/2)) angle=0-pi+angle
      newstom[[i]]$angle <- pracma::rad2deg(angle)
      }
    }
    stom_key = do.call(rbind, newstom)

  }
  #add corrected angle
  stom_key <- stom_key %>%
    dplyr::mutate(angle.corrected = angle/mean(angle, na.rm=TRUE)) %>%
    dplyr::relocate(angle.corrected, .after=angle)

  #add centroid
  suppressWarnings(stom_key[,c("centroid.x", "centroid.y")] <- sf::st_coordinates(sf::st_centroid(stom_key)))

  if(!is.null(subs.val)){
    bstom  <- sf::st_buffer(stom_key, dist=buffervalue)
    subs <- cells[which(cells$value==subs.val),]
    stom_subs <-  sf::st_intersects(bstom,subs)
    subskey <- data.frame(ID=as.numeric(rownames(subs)), n=1:nrow(subs))
    for (i in 1:length(stom_subs)){
      stom_subs[[i]] <- subskey %>%
        dplyr::inner_join(stom_subs[[i]] %>%
                            as.data.frame(), by=c("n"=".")) %>%
        dplyr::pull(ID)
    }
    stom_key$SUBS <- stom_subs
  }
return(stom_key)

}


#' Identify clustering in cells of type "other"
#'
#' @param cells Object of class \code{sf}. The image to be measured, converted
#' to polygons (e.g. using \code{image_to_poly()})
#' @param val Numeric. Value of cell typ eof interest.
#' @param buffervalue Amount to buffer cells by. Should be similar to wall width (which varies with image scaling).
#'
#' @return Numeric with size of median cluster (i.e. these cells are usually found in groups of X).
#' @export
#'
#'

other_cluster_size <- function(cells, val, buffervalue){
  other <- cells[which(cells$value==val),]
  b <- sf::st_buffer(other, dist=buffervalue*3) #wider buffer for these
  cl <- unique((sf::st_intersects(b,b)))
  return(median(lengths(cl)))
}
