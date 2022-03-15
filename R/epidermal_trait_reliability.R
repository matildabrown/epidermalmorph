#' Estimate trait reliability
#'
#' @param x A data.frame with the trait measurements and grouping variable.
#' Names of traits should be the same as those produced by \code{extract_epidermal_traits()}
#' @param grouping.variable Character. Name of the column used to group images (e.g.
#' individual, species)
#' @param heatmap Logical. Plot a heatmap of trait reliability? Defaults to FALSE.
#' @param thresholdval Numeric. Used for heatmap plotting only - what is the
#' threshold for an 'unreliable' trait? Defaults to 0.25 (i.e. where the
#' variation in one group exceeds 25 percent of the variation in the entire imageset)
#'
#' @return Either a list or data.frame, depending on the value of \code{heatmap}. If
#' \code{heatmap = TRUE}, then a list of length 2 is returned, with a data.frame and
#' the ggplot heatmap object. If \code{heatmap = FALSE}, only the data.frame is
#' returned
#' @export
#'
#' @details This function generates estimates of intra-species or intra-
#' individual variability, expressed as a proportion of the variation in the
#' entire imageset.
#'
#' @import ggplot2
#'
#' @examples
#' df <- epidermalmorph::podocarps
#' epidermal_trait_reliability(df, grouping.variable="individual")
#'
epidermal_trait_reliability <- function (x, grouping.variable, heatmap=FALSE, thresholdval=0.25){

  value <- variation <- trait <-  NULL
  names_to_pivot <- c("stomatal.density.px2","stomatal.index",
                      "dist.between.stom.rows","row.consistency",
                      "row.wiggliness","stom.angle.sd","stom.distNN.mean",
                      "stom.dist2NN.mean","stom.spacingNN.mean","stom.distNN.sd",
                      "stom.dist2NN.sd","stom_spacingNN.sd","stom.nsubcells.mean",
                      "stom.subsarea.mean","stom.gclength.mean","stom.AR.mean",
                      "stom.butterfly.mean","stom.symmetry.mean",
                      "stom.gclength.sd","stom.AR.sd","pavezone.angle.median",
                      "pavezone.AR.median","pavezone.area.median",
                      "pavezone.complexity.median","pavezone.endwallangle.mean",
                      "pavezone.endwalldiff.mean","pavezone.njunctionpts.mean",
                      "pavezone.undulation.amp.median",
                      "pavezone.undulation.freq.mean","stomzone.AR.median",
                      "stomzone.area.median","stomzone.complexity.median",
                      "stomzone.undulation.amp.median",
                      "stomzone.undulation.freq.mean","polar.AR.median",
                      "polar.area.median","polar.complexity.median",
                      "polar.undulation.amp.median","polar.undulation.freq.mean",
                      "pavezone.area.sd","pavezone.AR.sd",
                      "pavezone.njunctionpts.sd","pavezone.angle.sd",
                      "stomzone.area.sd","polar.area.sd",
                      "pavezone.complexity.sd","pavezone.undulation.amp.sd",
                      "pavezone.undulation.freq.sd")
  trait_vars <- which(colnames(x) %in% names_to_pivot)
  if(sum(abs(round(colMeans(x[,trait_vars]))), na.rm=TRUE)>0) stop("Variables should be scaled before assessing their reliability.")

  df_im2 <- x %>%
    tidyr::pivot_longer(cols = dplyr::all_of(trait_vars), names_to = "trait")

  df_summary <- df_im2 %>%
    dplyr::group_by_at(dplyr::vars(c(grouping.variable,"trait"))) %>%
    dplyr::summarise(mean=mean(value), variation=sd(value))

  by_individual <- df_summary %>%
    dplyr::group_by_at(dplyr::vars(grouping.variable)) %>%
    dplyr::summarise(mean=mean(variation, na.rm=TRUE))

  by_trait <- df_summary %>%
    dplyr::group_by(trait) %>%
    dplyr::summarise(mean=mean(variation))

  df_summary$groupvar <- dplyr::pull(df_summary[,grouping.variable])
  df_summary$groupvar <- factor(df_summary$groupvar, levels = dplyr::pull(by_individual[order(by_individual$mean),grouping.variable]))

  df_summary$trait <- factor(df_summary$trait, levels = by_trait$trait[order(by_trait$mean,decreasing = TRUE)])



  if(heatmap==TRUE){

    p <- ggplot(df_summary, aes_string( x="groupvar", y="trait", fill="variation")) +
      geom_tile()+
      geom_text(aes(label = round(variation, 2)), size=3, colour="gray30") +
      scale_fill_gradientn(
        colors=c("#CFE2f1","white","#FFCBCB"),
        values=scales::rescale(c(0,thresholdval,max(df_summary$variation, na.rm=TRUE))),
        limits=c(0,max(df_summary$variation)),
        name=paste0("Intra-", grouping.variable,"\nvariation (sd)")
      ) +
      ggpubr::theme_pubclean()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
      xlab(grouping.variable)

    return (list(df_summary, p))
  }else{
    return(df_summary)
  }

}
