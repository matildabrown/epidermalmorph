#' Measure epidermal traits from subsets of cells
#'
#' @param x x List of the form output by \code{extract_epidermal_traits()} (see documentation of that function for details).
#' @param n.runs Numeric. Number of times to sample each value of k.
#' @param k.cells Numeric. Stopping threshold for number of cells.
#' @param k.pavement Numeric. Stopping threshold for number of cells.
#' @param k.stomata Numeric. Stopping threshold for number of cells.
#' @param stom.val Numeric. Value of stomata (guard cells + pore). Defaults to 85.
#' @param subs.val Numeric. Value of subsidiary cells. Defaults to 170.
#' @param pave.val Numeric. Value of pavement cells. Defaults to 255
#' @param verbose Logical. Show progress of sampling? Defaults to TRUE.
#' @param scaling Numeric. Global standard deviation values to standardize the variation of each trait.
#' @param trait.subset Character. Specific traits to measure. Defaults to NULL (all traits measured).
#'
#' @return A data.frame with the differences in trait values (measured in standard
#' deviations if scaling is supplied) between each sampled cell patch and the whole image.
#' @export
#'
subsample_epidermal_traits <- function(x,
                                       n.runs=100,
                                       k.cells=c(50),
                                       k.pavement=NULL,
                                       k.stomata=NULL,
                                       scaling = NULL,
                                       trait.subset=NULL,
                                       stom.val=85,
                                       subs.val=170,
                                       pave.val=255,
                                       verbose=TRUE){





  meta.cols <- c("image.ID","n","k.cells", "k.pavement","k.stomata","rotation","n.cells","n.pavement","n.stomata","n.subsidiary","n.pavezone", "n.stomzone","n.polar")
  trait.cols <- c("stomatal.density", "stomatal.index",
                  "dist.between.stom.rows", "row.consistency",
                  "row.wiggliness", "stom.angle.sd", "stom.distNN.mean",
                  "stom.dist2NN.mean", "stom.spacingNN.mean",
                  "stom.distNN.sd", "stom.dist2NN.sd", "stom.spacingNN.sd",
                  "stom.nsubcells.mean", "stom.subsarea.mean", "stom.gclength.mean",
                  "stom.AR.mean", "stom.butterfly.mean", "stom.symmetry.mean",
                  "stom.gclength.sd", "stom.AR.sd", "pavezone.angle.median",
                  "pavezone.AR.median", "pavezone.area.median",
                  "pavezone.complexity.median", "pavezone.endwallangle.mean",
                  "pavezone.endwalldiff.mean", "pavezone.njunctionpts.mean",
                  "pavezone.undulation.amp.median",
                  "pavezone.undulation.freq.mean", "stomzone.AR.median",
                  "stomzone.area.median", "stomzone.complexity.median",
                  "stomzone.undulation.amp.median",
                  "stomzone.undulation.freq.mean",
                  "polar.AR.median", "polar.area.median",
                  "polar.complexity.median", "polar.undulation.amp.median",
                  "polar.undulation.freq.mean", "pavezone.area.sd",
                  "pavezone.AR.sd", "pavezone.njunctionpts.sd",
                  "pavezone.angle.sd", "stomzone.area.sd", "polar.area.sd",
                  "pavezone.complexity.sd", "pavezone.undulation.amp.sd",
                  "pavezone.undulation.freq.sd")

  cells <- x$polygons[which(x$polygons$edge==0),c("ID","value","geometry","dist_from_stom")]
  if(!is.null(trait.subset))

if(!is.null(k.cells)){
  kcells <- list()

  if(length(which(k.cells>nrow(cells)))!=0) {
    message("Some values of k.cells are larger than the total number of cells - removing these values")
    k.cells <- k.cells[which(k.cells<nrow(cells))]
  }

  df <- NULL
  if(verbose==TRUE) pb <- progress::progress_bar$new(total = (length(k.cells)*n.runs), current = "|")

  for (k in 1:length(k.cells)){


samples <- list()
for (i in 1:n.runs){
  if(verbose==TRUE)  pb$tick()
 sampleindexi <- patch_subsample(cells, k.cells=k.cells[[k]],
                                  k.pavement= k.pavement,
                                  k.stomata = k.stomata,
                                  stom.val =  stom.val,
                                  subs.val = subs.val,
                                  pave.val = pave.val)

   df_sample <- suppressWarnings(subsampling_trait_extraction(x, cells, sampleindexi, trait.subset=trait.subset))
   df_sample$k.cells[2] <- k
   df_sample$n[2] <- i

   df <- dplyr::bind_rows(df, df_sample)

}
}
}

#rearrange cols to have meta at the start, get rid of duplicated reference rows
df <- df %>%
  unique() %>%
  dplyr::select(dplyr::all_of(colnames(df)[c(which(colnames(df) %in% meta.cols),
                                      which(colnames(df) %in% trait.cols))]))

#convert to delta-values (scaled difference from the whole image)
for (i in 1:ncol(df)){
  if(colnames(df)[i] %in% trait.cols) {
    df[,i] <- (df[,i]-df[1,i])
  }
}

#scale to global mean/sd
if(!is.null(scaling)){
  for (i in 1:ncol(df)){
    if(colnames(df)[i] %in% names(scaling)) {
      scalesi <- which(names(scaling)==colnames(df)[i])
      df[,i] <- df[,i]/scaling[scalesi]
    }
  }
}
df$image.ID <- df$image.ID[1]
df$rotation<- df$rotation[1]
df <- df[-1,]
return(df)
}
