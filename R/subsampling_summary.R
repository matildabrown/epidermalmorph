#' Summarise and plot sampling effort for epidermal traits
#'
#' @param df A data.frame of the type output by \code{subsample_epidermal_traits}
#' @param traits Character vector of the traits to summarise.
#' @param thresholdvalue Numeric. Desired threshold for sample variation. Defaults to 0.2.
#' @param number.of.bins Numeric. How many bins should be used to summarise samples? Recommended to use the same number of bins as values of k in \code{subsample_epidermal_traits}.
#'
#' @return A list of length 2: the first element (a list) contains the ggplot objects and the second (also a list) contains summarised data.
#' @export
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' #This example takes quite a while to run - see vignette for more details
#'
#'sampled <-  subsample_epidermal_traits(traits1.1, n.runs = 100,
#'k.cells = c(100,150,200,250,300,350,400),
#'scaling =globalSDs,
#'trait.subset = reliable_traits)
#'
#'summ <- subsampling_summary(sampled,
#'                            traits=c("stomzone.area.median" ,
#'                                     "stomatal.index",
#'                                     "stomzone.undulation.amp.median",
#'                                      "pavezone.undulation.amp.median"))
#' }
#'
subsampling_summary <- function(df, traits,
                                     thresholdvalue=0.2, number.of.bins=6){

  #binding global variables
  n <- n.cells <- n.cells.binned <-  n.pavement <-  n.pavement.binned <- NULL
  n.pavezone <- n.pavezone.binned <-  n.polar <- n.polar.binned <- NULL
  n.stomata <- n.stomata.binned <- n.stomzone <- n.stomzone.binned <- NULL
  Trait <- value <- NULL

#bin the cell counts for variation estimation
df_long <- df %>% tidyr::pivot_longer(cols = dplyr::all_of(traits), names_to = "Trait")

shapetraits <- c("stom.nsubcells.mean","stom.subsarea.mean","stom.gclength.mean",
                 "stom.AR.mean","stom.butterfly.mean","stom.symmetry.mean",
                 "stom.gclength.sd","stom.AR.sd")

arrtraits <- c("stomatal.density","stomatal.index","dist.between.stom.rows",
               "row.consistency","row.wiggliness","stom.angle.sd",
               "stom.distNN.mean","stom.dist2NN.mean","stom.spacingNN.mean",
               "stom.distNN.sd","stom.dist2NN.sd","stom.spacingNN.sd")

pavezone_traits <- traits[which(grepl("pavezone",traits)==TRUE)]
stomzone_traits <- traits[which(grepl("stomzone",traits)==TRUE)]
polar_traits <- traits[which(grepl("polar",traits)==TRUE)]
stomatal_arr_traits <- traits[which(traits %in% arrtraits)]
stomatal_shape_traits <- traits[which(traits %in% shapetraits)]
pavement_traits <- c(pavezone_traits, stomzone_traits, polar_traits)

p <- list()
datalist <- list()

################################################################################
  #all traits
  d_all <- df_long %>%
     dplyr::mutate(n.cells.binned = Hmisc::cut2(n.cells, g=number.of.bins)) %>%
     dplyr::group_by(Trait, n.cells.binned) %>%
     dplyr::summarise(sd=sd(value), n=n())

  if (max(c(length(stomatal_arr_traits)),
      length(stomatal_shape_traits),
      length(pavezone_traits),
      length(stomzone_traits),
      length(polar_traits))>4){

    plotcols <- c(rep(nyellow(1),times=length(stomatal_arr_traits)),
                  rep(nred(1),times=length(stomatal_shape_traits)),
                  rep(nblue(1),times=length(pavezone_traits)),
                  rep( ngreen(1),times=length(stomzone_traits)),
                  rep( npurple(1),times=length(polar_traits)))

  } else {
    plotcols <- c(nyellow(length(stomatal_arr_traits)),
              nred(length(stomatal_shape_traits)),
              nblue(length(pavezone_traits)),
              ngreen(length(stomzone_traits)),
              npurple(length(polar_traits)))
    }


  d_all$Trait <- factor(d_all$Trait, levels=c(stomatal_arr_traits,stomatal_shape_traits,
                                            pavezone_traits, stomzone_traits,
                                            polar_traits))

  df <- datalist[[length(datalist)+1]] <- d_all
  title <- "All epidermal traits"
  xlabtext <- "Number of cells (binned)"


  n.traits=length(unique(df$Trait))


  (p[[length(p)+1]] <- ggplot(df, aes(x=n.cells.binned, y=sd, colour=Trait, shape=Trait))+
      scale_x_discrete() +
      annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
      geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
      geom_point(size=2,stroke = 1.5)+
      geom_path(aes(group=Trait), size=1)+
      scale_colour_manual(values=plotcols)+
      scale_shape_manual(values=1:n.traits) +
      guides(colour = guide_legend(override.aes = list(shape = 1:n.traits)))+
      expand_limits(y=0)+
      ylab("Variation in subsamples (sd)")+
      xlab(xlabtext)+
      ggtitle(title)+
      ggpubr::theme_pubr()+
      theme(legend.position = "right",
            text = element_text(size = 10),
            plot.title = element_text(size = 12, face="bold"))+
      ggpubr::rotate_x_text(angle=30)

  )


#####################################################

if(length(pavement_traits)!=0){
  #the pavezone traits
  d_pavement <- df_long %>%
     dplyr::mutate(n.pavement.binned = Hmisc::cut2(n.pavement, g=number.of.bins)) %>%
     dplyr::filter(Trait %in% pavement_traits) %>%
     dplyr::group_by(Trait, n.pavement.binned) %>%
     dplyr::summarise(sd=sd(value), n=n())

  if (max(c(length(pavezone_traits),
          length(stomzone_traits),
          length(polar_traits)))>4){

    plotcols <- c(rep(nblue(1),times=length(pavezone_traits)),
                  rep( ngreen(1),times=length(stomzone_traits)),
                  rep( npurple(1),times=length(polar_traits)))

  } else {
    plotcols <- c(nblue(length(pavezone_traits)),
                  ngreen(length(stomzone_traits)),
                  npurple(length(polar_traits)))
  }


  d_pavement$Trait <- factor(d_pavement$Trait, levels=c(pavezone_traits, stomzone_traits,
                                              polar_traits))

  df <- datalist[[length(datalist)+1]] <- d_pavement
  title <- "All pavement cell traits"
  xlabtext <- "Number of pavement cells (binned)"


  n.traits=length(unique(df$Trait))


  (p[[length(p)+1]] <- suppressWarnings(ggplot(df, aes(x=n.pavement.binned, y=sd, colour=Trait, shape=Trait))+
      scale_x_discrete() +
      annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
      geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
      geom_point(size=2,stroke = 1.5)+
      geom_path(aes(group=Trait), size=1)+
      scale_colour_manual(values=plotcols)+
      scale_shape_manual(values=1:n.traits) +
      guides(colour = guide_legend(override.aes = list(shape = 1:n.traits)))+
      scale_y_continuous(limits=c(0, max(df$sd)), expand=expansion(mult = c(0, 0.5),
                                                                      add = c(0, 0)))+
      ylab("Variation in subsamples (sd)")+
      xlab(xlabtext)+
      ggtitle(title)+
        ggpubr::theme_pubr()+
        theme(legend.position = "right",
              text = element_text(size = 10),
              plot.title = element_text(size = 12, face="bold"))+
        ggpubr::rotate_x_text(angle=30)
  ))


}
########################################################
  if(length(stomatal_arr_traits)!=0 & length(stomatal_shape_traits)!=0){
    #all stomatal traits
    d_stomata <- df_long %>%
      dplyr::mutate(n.stomata.binned = Hmisc::cut2(n.stomata, g=number.of.bins)) %>%
      dplyr::filter(Trait %in% c(stomatal_arr_traits,stomatal_shape_traits)) %>%
      dplyr::group_by(Trait, n.stomata.binned) %>%
      dplyr::summarise(sd=sd(value), n=n())

    if (max(c(length(stomatal_arr_traits)),
            length(stomatal_shape_traits))>4){

      plotcols <- c(rep(nyellow(1),times=length(stomatal_arr_traits)),
                    rep(nred(1),times=length(stomatal_shape_traits)))

    } else {
      plotcols <- c(nyellow(length(stomatal_arr_traits)),
                    nred(length(stomatal_shape_traits)))
    }


    d_stomata$Trait <- factor(d_stomata$Trait, levels=c(stomatal_arr_traits,stomatal_shape_traits))

    df <- datalist[[length(datalist)+1]] <- d_stomata
    title <- "All stomatal traits"
    xlabtext <- "Number of stomata (binned)"


    n.traits=length(unique(df$Trait))


    p[[length(p)+1]] <- ggplot(df, aes(x=n.stomata.binned, y=sd, colour=Trait, shape=Trait))+
                          scale_x_discrete() +
                          annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
                          geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
                          geom_point(size=2,stroke = 1.5)+
                          geom_path(aes(group=Trait), size=1)+
                          scale_colour_manual(values=plotcols)+
                          scale_shape_manual(values=1:n.traits) +
                          guides(colour = guide_legend(override.aes = list(shape = 1:n.traits)))+
                          scale_y_continuous(limits=c(0, max(df$sd)), expand=expansion(mult = c(0, 0.5),
                                                                                      add = c(0, 0)))+
                          ylab("Variation in subsamples (sd)")+
                          xlab(xlabtext)+
                          ggtitle(title)+
                          ggpubr::theme_pubr()+
                          theme(legend.position = "right",
                                 text = element_text(size = 10),
                                 plot.title = element_text(size = 12, face="bold"))+
                          ggpubr::rotate_x_text(angle=30)



  }

########################################################

if(length(stomatal_arr_traits)!=0){
  #the stomatal traits
  d_stomata_arr <- df_long %>%
     dplyr::mutate(n.stomata.binned = Hmisc::cut2(n.stomata, g=number.of.bins)) %>%
     dplyr::filter(Trait %in% stomatal_arr_traits) %>%
     dplyr::group_by(Trait, n.stomata.binned) %>%
     dplyr::summarise(sd=sd(value), n=n())

  df <- datalist[[length(datalist)+1]] <- d_stomata_arr
  title <- "Stomatal arrangement traits"
  xlabtext <- "Number of stomata (binned)"


     n.traits=length(unique(df$Trait))
     if(n.traits>4) {
       if(n.traits %in% c(5,6,9)) ncolval=3
       if(n.traits %in% c(7,8,10:16)) ncolval=4
       if(n.traits %in% c(5:8)) {nlineval=2; warning("When plotting >4 traits it may be difficult to distinguish colours/shapes")}
       if(n.traits %in% c(9:12)) {nlineval=3; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}
       if(n.traits %in% c(12:16)) {nlineval=4; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}

       # df$colvals <- rep(1:ncolval, times=nlineval)[1:n]
       # df$linevals <- rep(1:nlineval, each=ncolval)[1:n]
     } else {
       ncolval=n.traits
       nlineval=1
       # test$colvals <- 1:n.traits
       # test$linevals <- 1
     }


     p[[length(p)+1]] <- ggplot(df, aes(x=n.stomata.binned, y=sd, colour=Trait, shape=Trait))+
         scale_x_discrete() +
         annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
         geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
         geom_point(size=2,stroke = 1.5)+
         geom_path(aes(group=Trait, linetype=Trait), size=1)+
         scale_colour_manual(values=rep(nyellow(ncolval), times=nlineval))+
         scale_shape_manual(values=1:n.traits) +
         scale_linetype_manual(values=rep(nlineval, each=ncolval))+
         guides(colour = guide_legend(override.aes = list(shape = 1:n.traits, linetype = 1:nlineval)))+
         scale_y_continuous(limits=c(0, max(df$sd)), expand=expansion(mult = c(0, 0.5),
                                                                      add = c(0, 0)))+
         ylab("Variation in subsamples (sd)")+
         xlab(xlabtext)+
         ggtitle(title)+
       ggpubr::theme_pubr()+
       theme(legend.position = "right",
             text = element_text(size = 10),
             plot.title = element_text(size = 12, face="bold"))+
       ggpubr::rotate_x_text(angle=30)


}

#####################################################################################
if(length(stomatal_shape_traits)!=0){
  #the stomatal traits
  d_stomata_sha <- df_long %>%
    dplyr::mutate(n.stomata.binned = Hmisc::cut2(n.stomata, g=number.of.bins)) %>%
    dplyr::filter(Trait %in% stomatal_shape_traits) %>%
    dplyr::group_by(Trait, n.stomata.binned) %>%
    dplyr::summarise(sd=sd(value), n=n())

  df <- datalist[[length(datalist)+1]] <- d_stomata_sha
  title <- "Stomatal size and shape traits"
  xlabtext <- "Number of stomata (binned)"


  n.traits=length(unique(df$Trait))
  if(n.traits>4) {
    if(n.traits %in% c(5,6,9)) ncolval=3
    if(n.traits %in% c(7,8,10:16)) ncolval=4
    if(n.traits %in% c(5:8)) {nlineval=2; warning("When plotting >4 traits it may be difficult to distinguish colours/shapes")}
    if(n.traits %in% c(9:12)) {nlineval=3; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}
    if(n.traits %in% c(12:16)) {nlineval=4; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}


  } else {
    ncolval=n.traits
    nlineval=1

  }


  p[[length(p)+1]] <- ggplot(df, aes(x=n.stomata.binned, y=sd, colour=Trait, shape=Trait))+
      scale_x_discrete() +
      annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
      geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
      geom_point(size=2,stroke = 1.5)+
      geom_path(aes(group=Trait, linetype=Trait), size=1)+
      scale_colour_manual(values=rep(nred(ncolval), times=nlineval))+
      scale_shape_manual(values=1:n.traits) +
      scale_linetype_manual(values=rep(nlineval, each=ncolval))+
      guides(colour = guide_legend(override.aes = list(shape = 1:n.traits, linetype = 1:nlineval)))+
      scale_y_continuous(limits=c(0, max(df$sd)), expand=expansion(mult = c(0, 0.5),
                                                                   add = c(0, 0)))+
      ylab("Variation in subsamples (sd)")+
      xlab(xlabtext)+
      ggtitle(title)+
    ggpubr::theme_pubr()+
    theme(legend.position = "right",
          text = element_text(size = 10),
          plot.title = element_text(size = 12, face="bold"))+
    ggpubr::rotate_x_text(angle=30)

}

#####################################################################################



if(length(pavezone_traits)!=0){
  #the pavezone traits
 d_pavezone <- df_long %>%
    dplyr::mutate(n.pavezone.binned = Hmisc::cut2(n.pavezone, g=number.of.bins)) %>%
    dplyr::filter(Trait %in% pavezone_traits) %>%
    dplyr::group_by(Trait,n.pavezone.binned) %>%
    dplyr::summarise(sd=sd(value), n=n())

 df <- datalist[[length(datalist)+1]] <- d_pavezone
  title <- "Pavement cells (pavezone) traits"
  xlabtext <- "Number of pavezone cells (binned)"


  n.traits=length(unique(df$Trait))
  if(n.traits>4) {
    if(n.traits %in% c(5,6,9)) ncolval=3
    if(n.traits %in% c(7,8,10:16)) ncolval=4
    if(n.traits %in% c(5:8)) {nlineval=2; warning("When plotting >4 traits it may be difficult to distinguish colours/shapes")}
    if(n.traits %in% c(9:12)) {nlineval=3; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}
    if(n.traits %in% c(12:16)) {nlineval=4; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}


  } else {
    ncolval=n.traits
    nlineval=1

  }


  p[[length(p)+1]] <- ggplot(df, aes(x=n.pavezone.binned, y=sd, colour=Trait, shape=Trait))+
      scale_x_discrete() +
      annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
      geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
      geom_point(size=2,stroke = 1.5)+
      geom_path(aes(group=Trait, linetype=Trait), size=1)+
      scale_colour_manual(values=rep(nblue(ncolval), times=nlineval))+
      scale_shape_manual(values=1:n.traits) +
      scale_linetype_manual(values=rep(nlineval, each=ncolval))+
      guides(colour = guide_legend(override.aes = list(shape = 1:n.traits, linetype = 1:nlineval)))+
      scale_y_continuous(limits=c(0, max(df$sd)), expand=expansion(mult = c(0, 0.5),
                                                                   add = c(0, 0)))+
      ylab("Variation in subsamples (sd)")+
      xlab(xlabtext)+
      ggtitle(title)+
    ggpubr::theme_pubr()+
    theme(legend.position = "right",
          text = element_text(size = 10),
          plot.title = element_text(size = 12, face="bold"))+
    ggpubr::rotate_x_text(angle=30)


}
######################################################################

if(length(stomzone_traits)!=0){
  #the pavezone traits
  d_stomzone <- df_long %>%
     dplyr::mutate(n.stomzone.binned = Hmisc::cut2(n.stomzone, g=number.of.bins)) %>%
     dplyr::filter(Trait %in% stomzone_traits) %>%
     dplyr::group_by(Trait, n.stomzone.binned) %>%
     dplyr::summarise(sd=sd(value), n=n())

   df <- datalist[[length(datalist)+1]] <- d_stomzone
   title <- "Pavement cell (stomzone) traits"
   xlabtext <- "Number of stomzone cells (binned)"


   n.traits=length(unique(df$Trait))
   if(n.traits>4) {
     if(n.traits %in% c(5,6,9)) ncolval=3
     if(n.traits %in% c(7,8,10:16)) ncolval=4
     if(n.traits %in% c(5:8)) {nlineval=2; warning("When plotting >4 traits it may be difficult to distinguish colours/shapes")}
     if(n.traits %in% c(9:12)) {nlineval=3; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}
     if(n.traits %in% c(12:16)) {nlineval=4; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}


   } else {
     ncolval=n.traits
     nlineval=1

   }


   p[[length(p)+1]] <- ggplot(df, aes(x=n.stomzone.binned, y=sd, colour=Trait, shape=Trait))+
       scale_x_discrete() +
       annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
       geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
       geom_point(size=2,stroke = 1.5)+
       geom_path(aes(group=Trait, linetype=Trait), size=1)+
       scale_colour_manual(values=rep(ngreen(ncolval), times=nlineval))+
       scale_shape_manual(values=1:n.traits) +
       scale_linetype_manual(values=rep(nlineval, each=ncolval))+
       guides(colour = guide_legend(override.aes = list(shape = 1:n.traits, linetype = 1:nlineval)))+
       scale_y_continuous(limits=c(0, max(df$sd)), expand=expansion(mult = c(0, 0.5),
                                                                    add = c(0, 0)))+
       ylab("Variation in subsamples (sd)")+
       xlab(xlabtext)+
       ggtitle(title)+
     ggpubr::theme_pubr()+
     theme(legend.position = "right",
           text = element_text(size = 10),
           plot.title = element_text(size = 12, face="bold"))+
     ggpubr::rotate_x_text(angle=30)


   }

########################################################

if(length(polar_traits)!=0){
  #the pavezone traits
  d_polar <- df_long %>%
     dplyr::mutate(n.polar.binned = Hmisc::cut2(n.polar, g=number.of.bins)) %>%
     dplyr::filter(Trait %in% polar_traits) %>%
     dplyr::group_by(Trait, n.polar.binned) %>%
     dplyr::summarise(sd=sd(value), n=n())

   df <- datalist[[length(datalist)+1]] <- d_polar
   title <- "Pavement cell (polar) traits"
   xlabtext <- "Number of polar cells (binned)"


   n.traits=length(unique(df$Trait))
   if(n.traits>4) {
     if(n.traits %in% c(5,6,9)) ncolval=3
     if(n.traits %in% c(7,8,10:16)) ncolval=4
     if(n.traits %in% c(5:8)) {nlineval=2; warning("When plotting >4 traits it may be difficult to distinguish colours/shapes")}
     if(n.traits %in% c(9:12)) {nlineval=3; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}
     if(n.traits %in% c(12:16)) {nlineval=4; warning("When plotting >8 traits it may be difficult to distinguish colours/shapes")}


   } else {
     ncolval=n.traits
     nlineval=1

   }


   p[[length(p)+1]] <- ggplot(df, aes(x=n.polar.binned, y=sd, colour=Trait, shape=Trait))+
       scale_x_discrete() +
       annotate("rect", xmin=0, xmax=Inf, ymin=thresholdvalue, ymax=Inf,fill="gray90")+
       geom_hline(yintercept=thresholdvalue, linetype=2, col="gray70")+
       geom_point(size=2,stroke = 1.5)+
       geom_path(aes(group=Trait, linetype=Trait), size=1)+
       scale_colour_manual(values=rep(npurple(ncolval), times=nlineval))+
       scale_shape_manual(values=1:n.traits) +
       scale_linetype_manual(values=rep(nlineval, each=ncolval))+
       guides(colour = guide_legend(override.aes = list(shape = 1:n.traits, linetype = 1:nlineval)))+
       scale_y_continuous(limits=c(0, max(df$sd)), expand=expansion(mult = c(0, 0.5),
                                                                    add = c(0, 0)))+
       ylab("Variation in subsamples (sd)")+
       xlab(xlabtext)+
       ggtitle(title)+
     ggpubr::theme_pubr()+
     theme(legend.position = "right",
           text = element_text(size = 10),
           plot.title = element_text(size = 12, face="bold"))+
     ggpubr::rotate_x_text(angle=30)


   }

###########################################


return(list(plots=p, data=datalist))

}
