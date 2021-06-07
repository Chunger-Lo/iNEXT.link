# DataInfo.link ----------------------
#' Exhibit basic data information
#'
#' \code{DataInfo.link}: exhibits basic data information
#'
#' @param data a vector/matrix/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence"},
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).#' @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
#' @import tibble
#' @examples
#' \dontrun{
#' data(puerto.rico)
#' DataInfo.link(puerto.rico$data, diversity = 'TD', datatype="abundance")
#' DataInfo.link(puerto.rico$data, diversity = 'PD', datatype="abundance",
#' row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' }

#' @export

DataInfo.link <- function(data, diversity = 'TD', datatype = "abundance", row.tree = NULL,col.tree = NULL){

  if(diversity == 'PD'){
    table <- lapply(data, function(y){datainfphy(data = y, datatype = datatype,
                                                 row.tree = row.tree,col.tree = col.tree)})%>%
      do.call(rbind,.)
    rownames(table) <- names(data)
    table = tibble::rownames_to_column(table, var = "Assemblages")
  }else if(diversity == 'TD'){
    table <- lapply(data, function(y){datainf(data = y, datatype = datatype)})%>%do.call(rbind,.)
    rownames(table) <- names(data)
    table = tibble::rownames_to_column(table, var = "Assemblages")
  }
  return(table)

}

# SC.link ----
#' Sample Completeness main function
#'
#' \code{SC.link} Estimation of Sample Completeness with order q
#'
#' @param data a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr
#' @param q a integer vector for the order of Hill number\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or
#' species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).#'
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50
#' @param conf  positive number < 1 specifying the level of confidence interval, default is 0.95.\cr\cr
#' @return a matrix of estimated sample completeness with order q: \cr\cr
#'
#' @examples
#' data(Norfolk)
#' output = SC.link(Norfolk)
#' ggSC.link(output)
#'
#' @references
#' Chao,A.,Y.Kubota,D.Zelený,C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages.
#' @export
SC.link <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30, conf = 0.95){
  data_long <- lapply(data, function(tab){
    as.matrix(tab)%>%c()}
  )
  res = iNEXT.4steps::SC(data = data_long, q = q, datatype = datatype, nboot = nboot, conf = conf)
  return(res)
}

# ggSC.link -------------------------------------------------------------------
#' ggplot for Sample Completeness
#'
#' \code{ggSC.link} The figure for estimation of Sample Completeness with order q
#'
#' @param outcome a table generated from SC function
#' @return a figure of estimated sample completeness with order q
#'
#' @examples
#' data(Norfolk)
#' output = SC.link(Norfolk)
#' ggSC.link(output)

#' @export
ggSC.link <- function(outcome){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(outcome, aes(x = Order.q, y = Estimate.SC, colour = Assemblage)) +
    geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
    geom_ribbon(aes(ymin = SC.LCL, ymax = SC.UCL, fill = Assemblage),
                alpha = 0.2, linetype = 0) + scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Sample completeness") + theme(text = element_text(size = 12)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"), legend.title = element_blank())
}

# iNEXT.link -------------------------------------------------------------------
#' Interpolation (rarefaction) and extrapolation of Chao et al.’s (2021) network diversity and mean network diversity
#' Function \code{iNEXT.link} computes network diversity estimates for rarefied samples and extrapolated samples
#' along with confidence intervals and related coverage estimates based on Chao et al.’s (2021) network
#' diversity (ND)
#' @param data a matrix/data.frame of species
#' abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
#' Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
#' data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
#' Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
#' and thus N matrices, users must first merge the N matrices by species identity to obtain a large
#' merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
#' data. The row (species) names of data must match the species names in the phylogenetic tree and
#' thus cannot be missing.
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param endpoint a positive integer specifying the endpoint for the rarefaction and extrapolation range.
#' If \code{NULL}, then \code{endpoint} = double of the reference sample size in each assemblage. It is ignored if \code{size} is given.
#' @param knots a positive integer specifying the number of equally-spaced knots between 1 and the \code{endpoint}. Default is 40.
#' @param size a sequence of positive integers specifying the sample sizes for which PD or meanPD estimates will be calculated.
#' If \code{NULL}, then estimates will be calculated for those sample sizes determined by the specified/default \code{endpoint}
#' and \code{knots}.
#' @param PDtype Select phylogenetic diversity type: \code{PDtype = "PD"} for Chao et al. (2010) phylogenetic diversity and
#' \code{PDtype = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number).
#' It will be used when \code{diversity = 'PD'}. Default is \code{"PD"}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param col.tree phylogenetic tree of column assemblage in interaction matrix
#' @param row.tree phylogenetic tree of row assemblage in interaction matrix.
#' @import ape ggplot2 dplyr tidytree stats chaoUtility phytools iNEXT.3D
#' @importFrom phyclust get.rooted.tree.height
#' @return

#' \itemize{
#'  \item{\code{$DataInfo}: A dataframe summarizing data information}
#'  \item{\code{$iNextEst}: coverage-based diversity estimates along with confidence intervals}
#'  (if \code{nboot > 0}) for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#'  \item{\code{$AsyEst}: for
#' showing asymptotic diversity estimates along with related statistics.}
#' }
#'
#' @examples
#' \dontrun{
#' data(puerto.rico)
#' iNEXT.link(puerto.rico$data, diversity = 'TD', datatype="abundance")
#' iNEXT.link(puerto.rico$data, diversity = 'PD', datatype="abundance",
#'            row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' }

#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
#' \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export

iNEXT.link <- function(data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", size = NULL, nT = NULL,
                         endpoint = NULL, knots = 40, conf = 0.95, nboot = 30,
                         row.tree = NULL, col.tree = NULL, PDtype = 'meanPD'
                         ){
  # User interface
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(data)[1]

  if(datatype == "incidence"){
    stop('datatype="incidence" was no longer supported after v2.0.8,
         please try datatype="incidence_freq".')
  }
  if(datatype=="incidence_freq") datatype <- "incidence"

  if(datatype=="incidence_raw"){
    if(class_x=="list"){
      data <- lapply(data, as.incfreq)
    }else{
      data <- as.incfreq(data)
    }
    datatype <- "incidence"
  }
  if ( sum(!(diversity %in% c('TD', 'PD', 'FD', 'AUC')))>0 ){stop("Please select one of below diversity: 'TD', 'PD', 'FD'",
                                                                  call. = FALSE)}

  res = list()
  if(diversity == 'TD'){
    ## 1. datainfo
    datainfo = DataInfo.link(data = data, diversity = diversity, datatype = datatype)
    ## 2. iNterpolation/ Extrapolation
    data_long <- lapply(data, function(tab){
      as.matrix(tab)%>%c()}
    )
    INEXT_est <- iNEXT.3D::iNEXT3D(data_long, diversity = 'TD', q = q,conf = conf,
                                  nboot = nboot, knots = knots, endpoint = endpoint, size = size)

    res[[1]] = datainfo
    res[[2]] = INEXT_est$iNextEst
    res[[3]] = INEXT_est$AsyEst
    names(res) = c("DataInfo", "iNextEst", "AsyEst")

  }else if(diversity == 'PD'){
    ## 1. datainfo
    datainfo = DataInfo.link(data = data, diversity = diversity, datatype = datatype,
                             row.tree = row.tree, col.tree = col.tree)
    ## 2. iNterpolation/ Extrapolation
    data_long <- lapply(data, function(tab){
      as.matrix(tab)%>%c()}
      ## 拉長
    )
    NetiNE <- get.netphydiv_iNE(data = data, q = q,B = nboot,row.tree = row.tree,
                                col.tree = col.tree,conf = conf, knots = knots, PDtype = 'PD')
    ## 3. empirical and asymptotic diversity
    NetDiv <- get.netphydiv(data = data,q = q,B = nboot,row.tree = row.tree,col.tree = col.tree,conf = conf, PDtype = 'PD')

    res[[1]] = datainfo
    res[[2]] = NetiNE
    res[[3]] = NetDiv
    names(res) = c("DataInfo", "iNextEst", "AsyEst")

  }

  return(res)
}

# iNEXTbeta.link ---------------------------
#' Interpolation (rarefaction) and extrapolation of Chao et al.’s (2021) network diversity and mean network diversity

#' Function \code{iNEXTbeta.link} Interpolation and extrapolation of Beta diversity with order q
#'
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
#' Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
#' data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
#' Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
#' and thus N matrices, users must first merge the N matrices by species identity to obtain a large
#' merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
#' data. The row (species) names of data must match the species names in the phylogenetic tree and
#' thus cannot be missing.
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic',
#' and 'FD' = 'Functional' under certain threshold.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import ape
#' @import phytools
#' @import phyclust
#' @import tidytree
#' @import colorRamps
#' @import iNEXT.3D
#' @import future.apply
#' @import ade4
#' @import tidyr
#' @return A list of seven lists with three-diversity and four-dissimilarity.


#' @examples
#' \dontrun{
#' # example
#' data(puerto.rico)
#' beta1 = iNEXTbeta.link(data = puerto.rico$data, level = seq(0.5, 0.9, 0.4), datatype='abundance',q = c(0, 1, 2),
#'                        diversity = 'TD', nboot = 10, conf = 0.95)
#' beta2 = iNEXTbeta.link(networks = puerto.rico$data, level = seq(0.5, 0.9, 0.4), datatype='abundance',q = c(0, 1, 2),
#'                        data = 'PD', nboot = 10, conf = 0.95,
#'                        row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' }
#' @references
#' Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J.(2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters 8, 148-159. (pdf file) Spanish translation in pp. 85-96 of Halffter, G. Soberon, J., Koleff, P. and Melic, A. (eds) 2005 Sobre Diversidad Biologica: el Sognificado de las Diversidades Alfa, Beta y Gamma. m3m-Monografias 3ercer Milenio, vol. 4, SEA, CONABIO, Grupo DIVERSITAS & CONACYT, Zaragoza. IV +242 pp.
#' Chiu, C.-H., Jost, L. and Chao*, A. (2014). Phylogenetic beta diversity, similarity, and differentiation measures based on Hill numbers. Ecological Monographs 84, 21-44.
#' @export

iNEXTbeta.link = function(data, diversity = 'TD', level = seq(0.5, 1, 0.5), datatype=c('abundance', 'incidence_raw'),
                          q = c(0, 1, 2),nboot = 20, conf = 0.95,
                          row.tree = NULL,col.tree = NULL){

  if(class(data[[1]]) == 'data.frame' ){data = list(data); }
  combined_list = lapply(data, function(y){
    long = ready4beta(y)%>%filter_all(any_vars(. != 0))
    rownames(long) = rownames(long)%>%gsub('\\.','_',.)
    return(long)
  })
  if(diversity == 'TD'){
    dissimilarity <- iNEXTbeta(data = combined_list, diversity = 'TD',level = level, datatype = datatype,
                               q =q ,nboot = nboot, conf = conf)
  }
  else if(diversity == 'PD'){
    if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
    if(!is.null(row.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}

    dissimilarity = iNEXTbeta.PDlink(data = combined_list, level = level, datatype = datatype,
                                     q =q ,row.tree = row.tree,col.tree = col.tree, nboot = nboot)
  }

  return(dissimilarity)
}


# ggiNEXT.link -------------------------------------------------------------------
#' ggplot2 extension for outcome from \code{iNEXT.link}
#'
#' \code{ggiNEXT.link}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT.link}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param outcome a list object computed by \code{\link{iNEXT.link}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1});
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param facet.var create a separate plot for each value of a specified variable:
#'  no separation \cr (\code{facet.var="None"});
#'  a separate plot for each diversity order (\code{facet.var="Order.q"});
#'  a separate plot for each assemblage (\code{facet.var="Assemblage"});
#'  a separate plot for each combination of order x assemblage (\code{facet.var="Both"}).
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="None"});
#'  use different colors for diversity orders (\code{color.var="Order.q"});
#'  use different colors for sites (\code{color.var="Assemblage"});
#'  use different colors for combinations of order x assemblage (\code{color.var="Both"}).
#' @param grey a logical variable to display grey and white ggplot2 theme.
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' data(Norfolk)
#' out1 <- iNEXT.link(Norfolk, diversity = 'TD',datatype = "abundance", nboot = 0)
#' ggiNEXT.link(outcome = out1, type = 1)
#' ggiNEXT.link(outcome = out1, type = 2)
#' ggiNEXT.link(outcome = out1, type = 3)
#'
#' #' data(puerto.rico)
#' out2 <- iNEXT.link(puerto.rico$data, diversity = 'PD', datatype="abundance",
#' row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' ggiNEXT.link(outcome = out2, type = 1)
#' ggiNEXT.link(outcome = out2, type = 2)
#' ggiNEXT.link(outcome = out2, type = 3)
#' }
#' @export
ggiNEXT.link <- function(outcome, diversity = 'TD', type = 1,se = TRUE,facet.var = "Assemblage",
                         color.var = "Order.q", text.size = 12, stript.size = 12){
  if(diversity == 'TD'){
    iNEXT.3D::ggiNEXT3D(outcome, type = type)

  }else if(diversity == 'PD'){
    # output = outcome
    # output$iNextEst$size_based = output$iNextEst$size_based%>%
    #   rename('qD'="PD", 'qD.UCL'="PD.UCL",'qD.LCL'="PD.LCL")
    # iNEXT.3D::ggiNEXT3D(output, type = 1)
    iNE <- outcome$iNextEst
    iNE.sub <- iNE[iNE$method == "observed",]
    iNE[iNE$method == "observed",]$method <-  "interpolated"
    ex <- iNE.sub
    ex$method <- "extrapolated"
    iNE <- rbind(iNE,ex)
    iNE$method <- factor(iNE$method,levels = c("interpolated","extrapolated"))
    iNE$Order.q = paste0("q = ", iNE$Order.q)
    iNE.sub$Order.q = paste0("q = ", iNE.sub$Order.q)

    if(type == 1){
      # size-based
      plot <- ggplot(iNE, aes(x = m,y = PD)) + geom_line(aes(color = Region,linetype = method),size = 1.2) +
        facet_wrap(~Order.q, scales = "free") +
        # facet_grid()
        geom_ribbon(aes(x = m,ymax = PD.UCL ,ymin = PD.LCL,fill = Region),alpha = 0.25) + theme_bw()+
        geom_point(aes(x = m,y = PD ,color = Region,shape = Region),size = 5,data = iNE.sub) +
        xlab("Number of individuals")
    }else if(type == 2){
      # output <- outcome$size_based
      # if (length(unique(output$Order.q)) > 1) output <- subset(output, Order.q == unique(output$Order.q)[1])
      # output$y.lwr <- output$SC.LCL
      # output$y.upr <- output$SC.UCL
      # id <- match(c(x_name, "Method", "SC", "SC.LCL", "SC.UCL", "Assemblage", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(output), nomatch = 0)
      # output[,1:10] <- output[, id]
      #
      # xlab_name <- paste0("Number of ", xlab_name)
      # ylab_name <- "Sample Coverage"

    }else if(type == 3){
      # coverage-based
      plot <- ggplot(iNE) + geom_line(aes(x = SC,y = PD,color = Region,linetype = method),size = 1.2) +
        facet_wrap(~Order.q, scales = "free") +
        geom_ribbon(aes(x = SC,ymax = PD.UCL ,ymin = PD.LCL,fill = Region),alpha = 0.25) +
        geom_point(aes(x = SC,y = PD ,color = Region,shape = Region),size = 5,data = iNE.sub) + theme_bw() +
        xlab("Sample coverage")
    }

    plot +
      theme(legend.position = "bottom",
            legend.title=element_blank(), strip.text = element_text(size = stript.size),
            text=element_text(size=text.size),
            legend.key.width = unit(0.8,"cm"))  +
      labs(y = "Phylogenetic network diversity", lty = "Method")
    # if(grey){
    #   g <- g +
    #     scale_fill_grey(start = 0, end = .4) +
    #     scale_colour_grey(start = .2, end = .2)
    # }
  }
}


# ggiNEXT_beta.link -------------------------------------------------------------------
#' ggplot2 extension for outcome from \code{iNEXT.beta.link}
#'
#' \code{ggiNEXTbeta.link}: ggplot for Interpolation and extrapolation of Beta diversity with order q
#'
#' @param outcome the outcome from \code{"iNEXTbeta.link"}
#' @param type selection of plot type : \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ;
#' \code{type = 'D'} for plotting 4 turnover dissimilarities.
#' @param measurement character indicating the label of y-axis.
#' @param scale Are scales shared across all facets (\code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})? Default is \code{"free"}.
#' @param main The title of the plot.
#' @param transp a value between 0 and 1 controlling transparency. \code{transp = 0} is completely transparent, default is 0.4.
#'
#' @return a figure for Beta diversity or dissimilarity diversity.
#'
#'
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' data(Norfolk)
#' beta1 = iNEXTbeta.link(networks = puerto.rico$data, level = seq(0.5, 1, 0.5), datatype='abundance',q = c(0, 1, 2),
#'                        diversity = 'TD', nboot = 10, conf = 0.95)
#' beta2 = iNEXTbeta.link(networks = puerto.rico$data, level = seq(0.5, 1, 0.5), datatype='abundance',q = c(0, 1, 2),
#'                        diversity = 'PD', nboot = 10, conf = 0.95,
#'                        row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#'
#' ggiNEXTbeta.link(beta1,diversity = 'TD', type = 'B')
#' ggiNEXTbeta.link(beta1,diversity = 'TD', type = 'D')
#' ggiNEXTbeta.link(beta2,diversity = 'PD', type = 'B')
#' ggiNEXTbeta.link(beta2,diversity = 'PD', type = 'D')
#' }
#' @export

ggiNEXTbeta.link <- function(outcome, type = c('B', 'D'),
                             diversity = 'TD',
                             scale='free', main=NULL, transp=0.4, stript.size = 11, text.size = 13){
  if (type == 'B'){

    gamma = lapply(outcome, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(outcome, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
    beta =  lapply(outcome, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()
    beta = beta %>% filter(Method != 'Observed')
    beta[beta == 'Observed_alpha'] = 'Observed'

    df = rbind(gamma, alpha, beta)
    for (i in unique(gamma$Order)) df$Order[df$Order==i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))

    id_obs = which(df$Method == 'Observed')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$level = new$level - 0.0001
      new$Method = 'Interpolated'

      newe = df[id_obs[i],]
      newe$level = newe$level + 0.0001
      newe$Method = 'Extrapolated'

      df = rbind(df, new, newe)

    }

    if (diversity=='TD') { ylab = "Taxonomic diversity" }
    if (diversity=='PD') { ylab = "Phylogenetic diversity" }

  }
  if (type=='D'){

    C = lapply(outcome, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
    U = lapply(outcome, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
    V = lapply(outcome, function(y) y[["V"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
    S = lapply(outcome, function(y) y[["S"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()
    C = C %>% filter(Method != 'Observed')
    U = U %>% filter(Method != 'Observed')
    V = V %>% filter(Method != 'Observed')
    S = S %>% filter(Method != 'Observed')
    C[C == 'Observed_alpha'] = U[U == 'Observed_alpha'] = V[V == 'Observed_alpha'] = S[S == 'Observed_alpha'] = 'Observed'

    df = rbind(C, U, V, S)
    for (i in unique(C$Order)) df$Order[df$Order==i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("1-CqN","1-UqN","1-VqN","1-SqN"))

    id_obs = which(df$Method == 'Observed')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$level = new$level - 0.0001
      new$Method = 'Interpolated'

      newe = df[id_obs[i],]
      newe$level = newe$level + 0.0001
      newe$Method = 'Extrapolated'

      df = rbind(df, new, newe)

    }
    if (diversity=='TD') { ylab = "Taxonomic dissimilarity" }
    if (diversity=='PD') { ylab = "Phylogenetic dissimilarity" }

  }
  lty = c(Interpolated = "solid", Extrapolated = "dashed")
  # lty = c(Interpolated = "solid", Extrapolated = "twodash")
  df$Method = factor(df$Method, levels = c('Interpolated', 'Extrapolated', 'Observed'))

  double_size = unique(df[df$Method=="Observed",]$Size)*2
  double_extrapolation = df %>% filter(Method=="Extrapolated" & round(Size) %in% double_size)

  point_size = 2
  ggplot(data = df, aes(x = level, y = Estimate, col = Region)) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=transp) +
    geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) +
    scale_linetype_manual(values = lty) +
    # geom_line(lty=2) +
    geom_point(data = subset(df, Method=='Observed' & div_type=="Gamma"),shape=19, size=point_size) +
    geom_point(data = subset(df, Method=='Observed' & div_type!="Gamma"),shape=1, size=point_size,stroke=1.5)+
    geom_point(data = subset(double_extrapolation, div_type == "Gamma"),shape=17, size=point_size) +
    geom_point(data = subset(double_extrapolation, div_type!="Gamma"),shape=2, size=point_size,stroke=1.5) +
    facet_grid(div_type~Order, scales = scale) +
    # facet_wrap(div_type~Order, scales = scale, switch="both") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(x='Sample coverage', y=ylab, title=main)
}

# ggiNEXT_beta.link <- function(output, type = c('B', 'D'),
#                               diversity = 'TD',
#                               scale='free', main=NULL, transp=0.4, stript.size = 11, text.size = 13){
#   if (type == 'B'){
#
#       gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
#       alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
#       beta =  lapply(output, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()
#       beta = beta %>% filter(Method != 'Observed')
#       beta[beta == 'Observed_alpha'] = 'Observed'
#
#       df = rbind(gamma, alpha, beta)
#       for (i in unique(gamma$Order)) df$Order[df$Order==i] = paste0('q = ', i)
#       df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))
#
#       id_obs = which(df$Method == 'Observed')
#
#       for (i in 1:length(id_obs)) {
#
#         new = df[id_obs[i],]
#         new$level = new$level - 0.0001
#         new$Method = 'Interpolated'
#
#         newe = df[id_obs[i],]
#         newe$level = newe$level + 0.0001
#         newe$Method = 'Extrapolated'
#
#         df = rbind(df, new, newe)
#
#       }
#
#       if (diversity=='TD') { ylab = "Taxonomic diversity" }
#       if (diversity=='PD') { ylab = "Phylogenetic diversity" }
#
#     }
#
#     if (type=='D'){
#
#       C = lapply(output, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
#       U = lapply(output, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
#       V = lapply(output, function(y) y[["V"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
#       S = lapply(output, function(y) y[["S"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()
#       C = C %>% filter(Method != 'Observed')
#       U = U %>% filter(Method != 'Observed')
#       V = V %>% filter(Method != 'Observed')
#       S = S %>% filter(Method != 'Observed')
#       C[C == 'Observed_alpha'] = U[U == 'Observed_alpha'] = V[V == 'Observed_alpha'] = S[S == 'Observed_alpha'] = 'Observed'
#
#       df = rbind(C, U, V, S)
#       for (i in unique(C$Order)) df$Order[df$Order==i] = paste0('q = ', i)
#       df$div_type <- factor(df$div_type, levels = c("1-CqN","1-UqN","1-VqN","1-SqN"))
#
#       id_obs = which(df$Method == 'Observed')
#
#       for (i in 1:length(id_obs)) {
#
#         new = df[id_obs[i],]
#         new$level = new$level - 0.0001
#         new$Method = 'Interpolated'
#
#         newe = df[id_obs[i],]
#         newe$level = newe$level + 0.0001
#         newe$Method = 'Extrapolated'
#
#         df = rbind(df, new, newe)
#
#       }
#
#       if (diversity=='TD') { ylab = "Taxonomic dissimilarity" }
#       if (diversity=='PD') { ylab = "Phylogenetic dissimilarity" }
#
#     }
#
#     lty = c(Interpolated = "solid", Extrapolated = "dotted")
#     df$Method = factor(df$Method, levels = c('Interpolated', 'Extrapolated', 'Observed'))
#
#     double_size = unique(df[df$Method=="Observed",]$Size)*2
#     double_extrepolation = df %>% filter(Method=="Extrapolated" & round(Size) %in% double_size)
#
#     # ggplot(data = df, aes(x = level, y = Estimate)) +
#     plot = ggplot(data = df, aes(x = level, y = Estimate, col = Region)) +
#       # geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=transp) +
#       geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region), alpha=transp) +
#       geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) + scale_linetype_manual(values = lty) +
#       # geom_line(lty=2) +
#       geom_point(data = subset(df, Method=='Observed' & div_type=="Gamma"),shape=19, size=3) +
#       geom_point(data = subset(df, Method=='Observed' & div_type!="Gamma"),shape=1, size=3,stroke=1.5)+
#       geom_point(data = subset(double_extrepolation, div_type == "Gamma"),shape=17, size=3) +
#       geom_point(data = subset(double_extrepolation, div_type!="Gamma"),shape=2, size=3,stroke=1.5) +
#       facet_grid(div_type~Order, scales = scale) +
#       theme_bw() +
#       theme(legend.position = "bottom", legend.title = element_blank(),
#             strip.text = element_text(size = stript.size), text = element_text(size = text.size),
#             legend.key.width = unit(1,"cm")) +
#       labs(x='Sample coverage', y=ylab, title=main)
#     if(length(output) == 1){
#       plot = plot + guides(col=FALSE, fill = FALSE)
#     }
#     return(plot)
# }

# Asy.link -------------------------------------------------------------------
#' Asymptotic diversity q profile
#'
#' \code{Asy.link} The estimated and empirical diversity of order q
#'
#' @param data a matrix/data.frame of species
#' abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
#' Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
#' data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
#' Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
#' and thus N matrices, users must first merge the N matrices by species identity to obtain a large
#' merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
#' data. The row (species) names of data must match the species names in the phylogenetic tree and
#' thus cannot be missing.
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param col.tree phylogenetic tree of column assemblage in interaction matrix
#' @param row.tree phylogenetic tree of row assemblage in interaction matrix.
#' @return a table of Asymptotic network diversity q profile
#'
#' @examples
#' \dontrun{
#' ## Ex.1
#' data(Norfolk)
#' out1 <- Asy.link(Norfolk, diversity = 'TD', datatype = "abundance", nboot = 10)
#' ggAsy.link(out1)
#' ## Ex.2
#' data(puerto.rico)
#' out2 <- Asy.link(puerto.rico$data, diversity = 'PD', datatype = "abundance", nboot = 10,
#' row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' ggAsy.link(out2)
#' }
#' @export
Asy.link <- function(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30, conf = 0.95,
                     row.tree = NULL, col.tree = NULL){
  if(diversity == 'TD'){
    NetDiv <- lapply(1:length(data), function(i) {
      x = data[[i]]
      assemblage = names(data)[[i]]
      # please note that nboot has to larger than 0
      res = MakeTable_Proposeprofile(data = x, B = nboot, q, conf = conf)%>%
        mutate(Network = assemblage, method = "Estimated")%>%
        filter(Target == "Diversity")%>%select(-Target, -`s.e.`)%>%
        rename("qD"="Estimate", "qD.LCL"="LCL", "qD.UCL"="UCL", 'Method' = 'method')

      return(res)
    })%>%do.call("rbind",.)
    return(NetDiv)
  }else if(diversity == 'PD'){
    NetDiv <- get.netphydiv(data = data,q = q,B = nboot,row.tree = row.tree,col.tree = col.tree,conf = conf)%>%
      filter(method == "Estimate")%>%
      mutate(method = ifelse(method == 'Estimate','Estimated',method ))%>%
      dplyr::select(Order.q, Estimate, LCL, UCL, Region, method )%>%
      set_colnames(c('Order.q', 'qD', 'qD.LCL','qD.UCL', 'Network', 'Method'))

    return(NetDiv)
  }
}


# Obs.link -------------------------------------------------------------------
#' Empirical diversity q profile
#'
#' \code{Obs.link} The estimated and empirical diversity of order q
#'
#' @param data a matrix/data.frame of species
#' abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
#' Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
#' data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
#' Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
#' and thus N matrices, users must first merge the N matrices by species identity to obtain a large
#' merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
#' data. The row (species) names of data must match the species names in the phylogenetic tree and
#' thus cannot be missing.
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param col.tree phylogenetic tree of column assemblage in interaction matrix
#' @param row.tree phylogenetic tree of row assemblage in interaction matrix.
#' @return a table of Empirical network diversity q profile
#'
#' @examples
#' \dontrun{
#' ## Example for abundance-based data
#' ## Ex.1
#' data(Norfolk)
#' out1 <- Obs.link(Norfolk, diversity = 'TD', datatype = "abundance", nboot = 30)
#' ggObs.link(out1)
#' ## Ex.2
#' data(puerto.rico)
#' out2 <- Obs.link(data = puerto.rico$data, diversity = 'PD', datatype = "abundance",
#" nboot = 10, row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' ggObs.link(out2)
#' }
#' @export
Obs.link <- function(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95,
                     col.tree = NULL, row.tree = NULL){
  if(diversity == 'TD'){
    NetDiv <- lapply(1:length(data), function(i){
      x = data[[i]]
      assemblage = names(data)[[i]]
      tmp <- c(as.matrix(x))
      ## nboot has to larger than 0
      res = MakeTable_Empericalprofile(data = x, B = nboot, q, conf = conf)%>%
        mutate(Network = assemblage, method = "Empirical")%>%
        filter(Target == "Diversity")%>%select(-Target, -`s.e.`)%>%
        rename("qD"="Emperical", "qD.LCL"="LCL", "qD.UCL"="UCL", 'Method' = 'method')
      return(res)
    })%>%do.call("rbind",.)
    return(NetDiv)
  }

  else if(diversity == 'PD'){
    NetDiv <- get.netphydiv(data = data,q = q,B = nboot,row.tree = row.tree,
                            col.tree = col.tree,conf = conf)%>%
      filter(method == "Empirical")%>%
      dplyr::select(Order.q, Estimate, LCL, UCL, Region, method )%>%
      set_colnames(c('Order.q', 'qD', 'qD.LCL','qD.UCL', 'Network', 'Method'))

    return(NetDiv)
  }
}


# ggObs.link -------------------------------------------------------------------
#' ggplot for Empirical Network diversity
#'
#' \code{ggObs.link} Plots q-profile based on the outcome of \code{Obs.link} using the ggplot2 package.\cr
#' @param outcome the outcome of the functions \code{Obs.link} .\cr
#' @param text.size control the text size of the output plot.
#' @return a figure of estimated sample completeness with order q\cr\cr
#'
#' @examples
#' \dontrun{
#' ## Ex.1
#' data(Norfolk)
#' out1 <- Obs.link(Norfolk, diversity = 'TD', datatype = "abundance", nboot = 10)
#' ggObs.link(out1)
#' ## Ex.2
#' data(puerto.rico)
#' out2 <- Obs.link(puerto.rico$data, diversity = 'PD', datatype = "abundance",
#' nboot = 10, row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' ggObs.link(out2)
#' }

#' @export
ggObs.link <- function(outcome, text.size = 14){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  # if (sum(unique(outcome$method) %in% c("Estimated", "Empirical")) == 0)
  #   stop("Please use the outcome from specified function 'AsyD'")

  ggplot(outcome, aes(x = Order.q, y = qD, colour = Network, lty = Method)) +
    geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = outcome[outcome$method == "Empirical", ],
                aes(ymin = qD.LCL, ymax = qD.UCL, fill = Network), alpha = 0.2, linetype = 0) +
    scale_fill_manual(values = cbPalette) +
    scale_linetype_manual(values = c(Estimated = 1, Empirical = 2)) +
    labs(x = "Order q", y = "Network diversity") + theme(text = element_text(size = 10)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.1, "cm"), legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-10,-10, -5, -10),
          text = element_text(size = text.size)
    )

  # +labs(x = "Order q", y = "Network phylogenetic diversity", lty = "Method") + scale_linetype_manual(values=c("dashed","solid"))
}
# ggAsy.link -------------------------------------------------------------------
#' ggplot for Asymptotic Network diversity
#'
#' \code{ggAsy.link} Plots q-profile based on the outcome of \code{Asy.link} using the ggplot2 package.\cr
#'
#' @param outcome the outcome of the functions \code{Asy.link} .\cr
#' @param text.size control the text size of the output plot.
#' @return a figure of estimated sample completeness with order q\cr\cr
#'
#' @examples
#' \dontrun{
#' ## Ex.1
#' data(Norfolk)
#' out1 <- Asy.link(Norfolk, diversity = 'TD', datatype = "abundance", nboot = 10)
#' ggAsy.link(out1)
#' ## Ex.2
#' data(puerto.rico)
#' out2 <- Asy.link(puerto.rico$data, diversity = 'PD', datatype = "abundance",
#'                  nboot = 10, row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#' ggAsy.link(out2)
#' }
#' @export
ggAsy.link <- function(outcome, text.size = 14){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  # if (sum(unique(outcome$method) %in% c("Estimate", "Empirical")) == 0)
  #   stop("Please use the outcome from specified function 'AsyD'")

  ggplot(outcome, aes(x = Order.q, y = qD, colour = Network, lty = Method)) +
    geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = outcome[outcome$method == "Estimated", ],
                aes(ymin = qD.LCL, ymax = qD.UCL, fill = Network), alpha = 0.2, linetype = 0) +
    scale_fill_manual(values = cbPalette) +
    scale_linetype_manual(values = c(Estimated = 1, Empirical = 2)) +
    labs(x = "Order q", y = "Network diversity") + theme(text = element_text(size = 10)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.1, "cm"), legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-10,-10, -5, -10),
          text = element_text(size = text.size)
    )

  # +labs(x = "Order q", y = "Network phylogenetic diversity", lty = "Method") + scale_linetype_manual(values=c("dashed","solid"))
}

# estimateD.link  -------------------------------------------------------------------
#' Compute species diversity with a particular of sample size/coverage
#'
#' \code{estimateD.link} computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#'
#' @param data a \code{matrix}, \code{data.frame} (species by assemblages), or \code{list} of species abundance/incidence raw data.\cr
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold. Besides,'AUC' is the fourth choice which
#' integrates several threshold functional diversity to get diversity.
#' @param q a numerical vector of the order of Hill number.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50
#' @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1).
#' If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes.
#' If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param col.tree phylogenetic tree of column assemblage in interaction matrix
#' @param row.tree phylogenetic tree of row assemblage in interaction matrix.
#' @return a data.frame of species diversity table including the sample size, sample coverage, method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#'
#' @examples
#' \dontrun{
#' data(Norfolk)
#' out1 <- estimateD.link(Norfolk, diversity = 'TD',datatype="abundance",
#'                        base="coverage", level=0.7, nboot = 30)
#' out2 <- estimateD.link(Norfolk, diversity = 'TD',datatype="abundance",
#'                        base="size", level=0.7, nboot = 30)
#' }
#' @export
estimateD.link = function(data, diversity = 'TD', q = c(0, 1, 2),datatype = "abundance",base = "size",
                          level = NULL,nboot = 50,conf = 0.95,
                          row.tree = NULL, col.tree = NULL){
  if(diversity == 'TD'){
    div = lapply(1:length(data), function(i){
      x = data[[i]]
      assemblage = names(data)[[i]]
      long = as.matrix(x)%>%c()
      iNEXT.3D::estimate3D(long, q=q,datatype=datatype, base=base,
                          diversity = 'TD', nboot = nboot,conf=conf)%>%
        mutate(Assemblage = assemblage)
    })%>%do.call("rbind",.)

    return(div)
  }else if(diversity == 'PD'){

    if(datatype=='abundance'){

      if(class(data)=="data.frame" | class(data)=="matrix"| class(data)=="integer" ) data = list(Region_1 = data)

      if(class(data)== "list"){
        if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
        Ns = sapply(data, ncol)
        data_list = data
      }

    }
    if(is.null(conf)) conf = 0.95
    tmp = qnorm(1 - (1 - conf)/2)

    for_each_region = function(data_2d, region_name, N){
      if (datatype=='abundance') {

        n = sum(data_2d)
        if(base == 'coverage'){
          size_m = sapply(level, function(i) iNEXT.beta:::coverage_to_size(data_2d, i, datatype='abundance'))
        }else if(base == 'size'){
          size_m = level
        }


        ref= iNEXT.3D:::Coverage(data_2d,m= n, datatype = 'abundance')
        #
        aL_table = create.aili(data_2d, row.tree = row.tree, col.tree = col.tree) %>%
          select(branch.abun, branch.length, tgroup)%>%
          filter(branch.abun>0)

        qPDm <- iNEXTPD2:::PhD.m.est(ai = aL_table$branch.abun,
                                     Lis = aL_table$branch.length%>%as.matrix(),
                                     m = size_m,
                                     q = q,nt = n,cal = 'PD')%>%
          as.vector()

        ## boot
        tbar <- sum(aL_table$branch.length*aL_table$branch.abun)/n

        if(nboot >1 ){
          boot.sam <- sample.boot.phy(data_2d,nboot,row.tree = row.tree,col.tree = col.tree)
          PD.sd <- lapply(boot.sam, function(aL_boot){
            tmp = iNEXTPD2:::PhD.m.est(ai = aL_boot$branch.abun,
                                       Lis = aL_boot$branch.length%>%as.matrix(),
                                       m = size_m,
                                       q = q,nt = n,cal = 'PD')%>%
              as.vector()%>%as.data.frame()
            return(tmp)
          })%>%
            abind(along=3) %>% apply(1:2, sd)%>%as.vector()
        }else{
          PD.sd = c(0,0,0)
        }






        ##
        len = length(q)
        res = data.frame(Assemblage = rep(region_name,len),
                   goalSC = rep(level, len),
                   SC = rep(iNEXT.3D:::Coverage(data_2d, datatype='abundance', size_m), len),
                   m = rep(size_m,len),
                   Method = ifelse(level > ref, 'Extrapolated', 'Interpolated'),
                   Order.q = q,
                   qPD = qPDm,
                   qPD.UCL = qPDm+1.96*PD.sd,
                   qPD.LCL = qPDm-1.96*PD.sd
                   )
        return(res)
      }
    }
    print(length(data))
    output = lapply(1:length(data), function(i) for_each_region(data = data_list[[i]],
                                                                region_name = region_names[i], N = Ns[i]))%>%
      do.call('rbind',.)

    return(output)
  }

}
### old version
# estimateD.link = function(data, diversity = 'TD', q = c(0, 1, 2),datatype = "abundance",base = "size",
#                           level = NULL,nboot = 50,conf = 0.95,
#                           row.tree = NULL, col.tree = NULL){
#   if(diversity == 'TD'){
#     div = lapply(1:length(data), function(i){
#       x = data[[i]]
#       assemblage = names(data)[[i]]
#       long = as.matrix(x)%>%c()
#       iNEXT.3D::estimate3D(long, q=q,datatype=datatype, base=base,
#                           diversity = 'TD', nboot = nboot,conf=conf)%>%
#         mutate(Assemblage = assemblage)
#     })%>%do.call("rbind",.)
#
#     return(div)
#   }else if(diversity == 'PD'){
#
#     q <- unique(ceiling(q))
#     ci <- qnorm(conf/2+0.5)
#
#     if (is.null(level)) {
#       if (datatype == "abundance") {
#         level <- sapply(data, function(x) {
#           ni <- sum(x)
#           Coverage(data = x, datatype = datatype, m = 2 * ni, nt = ni)
#         })
#       }
#       else if (datatype == "incidence_raw") {
#         # level <- sapply(data, function(x) {
#         #   ni <- ncol(x)
#         #   Coverage(data = x, datatype = datatype, m = 2 *
#         #              ni, nt = ni)
#         # })
#       }
#       level <- min(level)
#     }
#     # level =0.7
#     res = lapply(1:length(data), function(i){
#       x = data[[i]]
#       assemblage = names(data)[[i]]
#       long = as.matrix(x)%>%c()
#       m_target = coverage_to_size(x, C = level, datatype = 'abundance')
#       # if(base == "size"){m_target = level}
#       inex <- function(data,m,q,B,row.tree = NULL,col.tree = NULL) {
#         data <- as.matrix(data)
#         n <- sum(data)
#         phydata <- create.aili(data,row.tree = row.tree,col.tree = col.tree)
#         ## why tbar can be calculated by this?
#         ## (ai * Li) / n
#         tbar <- sum(phydata$branch.length*phydata$branch.abun)/n
#         boot.sam <- sample.boot.phy(data,B,row.tree = row.tree,col.tree = col.tree)
#         sc <- PhD:::Coverage(data,datatype = "abundance",m,nt =n)
#         # sc <- coverage(data,m)
#         sc.sd <- lapply(boot.sam,function(x){
#           x <- x[x$tgroup == "Tip",]$branch.abun
#           PhD:::Coverage(x,datatype = "abundance",m,nt =n)
#         })
#         sc.sd <- do.call(cbind,sc.sd)
#         sc.sd <- sapply(1:length(m), function(x){
#           sd(sc.sd[x,])
#         })
#         sc.table <- data.frame(m=m,SC = sc, SC.UCL = sc+ci * sc.sd,SC.LCL = sc - ci * sc.sd)
#         out <- lapply(q, function(x){
#           PD <- lapply(m,function(y){
#             my_PhD.m.est(ai = phydata$branch.abun, Lis = phydata$branch.length, m = y, q = x, nt = n, cal = 'PD')/tbar
#             # PhD:::PhD.m.est(phydata,y,x,datatype = "abundance",nt = n)/tbar
#           })%>%unlist()
#           PD.sd <- lapply(boot.sam, function(z){
#             tmp <- lapply(m,function(y){
#               # PhD:::PhD.m.est(z,y,x,datatype = "abundance",nt = n)/tbar
#               my_PhD.m.est(ai = z$branch.abun, Lis = z$branch.length, m = y, q = x, nt = n, cal = 'PD')/tbar
#             })
#             unlist(tmp)
#           })
#           PD.sd <- do.call(cbind,PD.sd)
#           PD.sd <- sapply(1:length(m), function(x){
#             sd(PD.sd[x,])
#           })
#           PD.table <- data.frame(m=m,method = ifelse(m<n,"interpolated",ifelse(n == m,"observed","extrapolated")),
#                                  Order.q = x,PD = PD, PD.UCL = PD+ci * PD.sd,PD.LCL = PD - ci * PD.sd)
#           out <- left_join(PD.table,sc.table)
#           out
#         })
#         do.call(rbind,out)
#       }
#
#       inex(data = x,m = m_target, q,B = nboot,row.tree,col.tree)%>%
#         mutate(Assemblage = assemblage)
#
#     })%>%do.call("rbind",.)
#     return(res)
#   }
# }

# Spec.link  -------------------------------------------------------------------
#' Specialization Estimation of Evenness with order q
#'
#' \code{Spec.link} computes Evenness Estimation of Evenness with order q.
#'
#' @param outcome the outcome of the functions \code{ObsND} .\cr
#' @return A list of estimated(empirical) evenness with order q.
#' Different lists represents different classes of Evenness.
#' Each list is combined with order.q and sites.
#' If "method" is estimated, then fist list will be named "C" which means the maximum standardized coverage between all double reference sample size.
#
# $summary individual summary of 4 steps of data.
#'
#' @examples
#' \dontrun{
#' data(Norfolk)
#' Est <- Spec.link(x = Norfolk, datatype = "abundance", q = c(0,1,2),
#' nboot = 30, method = "Estimated")
#' Emp <- Spec.link(x = Norfolk, datatype = "abundance", q = c(0,1,2),
#' nboot = 30, method = "Empirical")
#' Est
#' Emp
#' ggSpec(Est)
#' ggSpec(Emp)
#' }
#' @export

Spec.link <- function(data, q = seq(0, 2, 0.2),
                      diversity = 'TD',
                      datatype = "abundance",
                      method = "Estimated",
                      nboot = 30,
                      conf = 0.95,
                      E.class = c(1:5),
                      C = NULL){
  if (diversity == 'TD'){
    long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})

    Spec <- lapply(E.class, function(e){
      each_class = lapply(seq_along(long), function(i){
        res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                    method = method, nboot=nboot, E.class = e, C = C)
        res['Coverage'] = NULL
        res = lapply(res, function(each_class){
          each_class%>%
            mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
            rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
            mutate(Assemblage = names(long)[[i]])
        })
        # if(method == "Empirical") index = 1
        # if(method == "Estimated") index = 2
        # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
        return(res[[1]])
      })%>%do.call("rbind",.)

      each_class%>%mutate(class = paste0("E",e))
    })
    names(Spec) = paste0("E",E.class)
    return(Spec)

  }else if (diversity == 'PD'){


    long = lapply(data, function(da){
      da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%column_to_rownames('col_sp')
    })
    names(long) = names(data)

    Spec <- lapply(E.class, function(e){
      each_class = lapply(seq_along(data), function(i){
        res = Spec.PD(data[[i]], q = q,datatype = datatype,
                      method = method, nboot=nboot, E.class = e, C = C)

        res['Coverage'] = NULL
        res = lapply(res, function(each_class){
          each_class%>%
            mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
            rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
            mutate(Assemblage = names(long)[[i]])
        })
        # if(method == "Empirical") index = 1
        # if(method == "Estimated") index = 2
        # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
        return(res[[1]])
      })%>%do.call("rbind",.)

      each_class%>%mutate(class = paste0("E",e))
    })



    spec_PD <- function(){
      if (datatype == "abundance") {
        qD <- Evenness.profile(data, q, "abundance", method,
                               E.class, C)
        qD <- map(qD, as.vector)
        if (nboot > 1) {
          Prob.hat <- lapply(1:length(data), function(i) iNEXT.3D:::EstiBootComm.Ind(data[[i]]))
          Abun.Mat <- lapply(1:length(data), function(i) rmultinom(nboot,
                                                                   sum(data[[i]]), Prob.hat[[i]]))
          error = apply(matrix(sapply(1:nboot, function(b) {
            dat = lapply(1:length(Abun.Mat), function(j) Abun.Mat[[j]][,
                                                                       b])
            names(dat) = paste("Site", 1:length(dat), sep = "")
            dat.qD = Evenness.profile(dat, q, "abundance",
                                      method, E.class, C)
            unlist(dat.qD)
          }), nrow = length(q) * length(E.class) * length(Abun.Mat)),
          1, sd, na.rm = TRUE)
          error = matrix(error, ncol = length(E.class))
          se = split(error, col(error))
        }
        else {
          se = lapply(1:length(E.class), function(x) NA)
        }
        out <- lapply(1:length(E.class), function(k) {
          tmp = data.frame(Order.q = rep(q, length(data)),
                           Evenness = as.vector(qD[[k]]), s.e. = as.vector(se[[k]]),
                           Even.LCL = as.vector(qD[[k]] - qnorm(1 - (1 -
                                                                       conf)/2) * se[[k]]), Even.UCL = as.vector(qD[[k]] +
                                                                                                                   qnorm(1 - (1 - conf)/2) * se[[k]]), Assemblage = rep(names(data),
                                                                                                                                                                        each = length(q)), Method = rep(method, length(q) *
                                                                                                                                                                                                          length(data)))
          tmp$Even.LCL[tmp$Even.LCL < 0] <- 0
          tmp
        })
        if (is.null(C) == TRUE)
          C = unique(estimate3D(data, diversity = "TD", q = 0,
                                datatype = "abundance", base = "coverage", nboot = 0)$goalSC)
        if (method == "Estimated") {
          out <- append(C, out)
        }
      }
    }

    Spec <- lapply(E.class, function(e){
      each_class = lapply(seq_along(long), function(i){
        res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                     method = method, nboot=nboot, E.class = e, C = C)
        res['Coverage'] = NULL
        res = lapply(res, function(each_class){
          each_class%>%
            mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
            rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
            mutate(Assemblage = names(long)[[i]])
        })
        # if(method == "Empirical") index = 1
        # if(method == "Estimated") index = 2
        # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
        return(res[[1]])
      })%>%do.call("rbind",.)

      each_class%>%mutate(class = paste0("E",e))
    })
    # ### not finished yet
    # long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})
    #
    # Spec <- lapply(E.class, function(e){
    #   each_class = lapply(seq_along(long), function(i){
    #     res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
    #                                 method = method, nboot=nboot, E.class = e, C = C)
    #     res['Coverage'] = NULL
    #     res = lapply(res, function(each_class){
    #       each_class%>%
    #         mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
    #         rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
    #         mutate(Assemblage = names(long)[[i]])
    #     })
    #     # if(method == "Empirical") index = 1
    #     # if(method == "Estimated") index = 2
    #     # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
    #     return(res[[1]])
    #   })%>%do.call("rbind",.)
    #
    #   each_class%>%mutate(class = paste0("E",e))
    # })
    # names(Spec) = paste0("E",E.class)
    # return(Spec)
  }
}

# ggSpec.link -------------------------------------------------------------------
#' ggplot for Evenness
#
#' \code{ggSpec.link} The figure for estimation of Evenness with order q\cr
#'
#' @param outcome a table generated from \code{Spec.link} function\cr
#' @return a figure of estimated sample completeness with order q\cr
#'
#' @examples
#' \dontrun{
#' data(Norfolk)
#' Est <- Spec.link(data = Norfolk, diversity = 'TD', datatype = "abundance", q = c(0,1,2),
#'                  nboot = 30, method = "Estimated")
#' Emp <- Spec.link(data = Norfolk, diversity = 'TD', datatype = "abundance", q = c(0,1,2),
#'                  nboot = 30, method = "Empirical")
#' ggSpec(output = Est)
#' ggSpec(output = Emp)
#' }
#' @references
#' Chao,A.and Ricotta,C.(2019).Quantifying evenness and linking it to diversity, beta diversity, and similarity.
#' @export

ggSpec.link <- function(outcome){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  classdata = cbind(do.call(rbind, outcome))

  fig = classdata%>%
    ggplot(aes(x=Order.q, y=Specialization, colour=Assemblage, lty = Method)) +
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = classdata %>% filter(Method=="Estimated"),
                aes(ymin=Spec.LCL, ymax=Spec.UCL, fill=Assemblage),
                alpha=0.2, linetype=0) +
    geom_ribbon(data = classdata %>% filter(Method=="Empirical"),
                aes(ymin=Spec.LCL, ymax=Spec.UCL, fill=Assemblage),
                alpha=0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    labs(x="Order q", y="Specialization") +
    # theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2,"cm"),
          legend.title = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-5,-10),
          text = element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt")
    )

  if (length(outcome) != 1) fig = fig +
    facet_wrap(~class) +
    theme(strip.text.x = element_text(size=12, colour = "purple", face="bold"))

  return(fig)
}


even.class = function (q, qD, S, E.class)
{
  tmp = c()
  if (E.class == 1)
    tmp = ifelse(q != 1, (1 - qD^(1 - q))/(1 - S^(1 - q)),
                 log(qD)/log(S))
  if (E.class == 2)
    tmp = ifelse(q != 1, (1 - qD^(q - 1))/(1 - S^(q - 1)),
                 log(qD)/log(S))
  if (E.class == 3)
    tmp = (qD - 1)/(S - 1)
  if (E.class == 4)
    tmp = (1 - 1/qD)/(1 - 1/S)
  if (E.class == 5)
    tmp = log(qD)/log(S)
  return(tmp)
}

# iNEXTbeta.PDlink ----
iNEXTbeta.PDlink <- function(data, level, datatype='abundance', q = c(0, 1, 2),
                             nboot = 20, conf = 0.95,
                             row.tree = NULL,col.tree = NULL){
  max_alpha_coverage = F

  if(datatype=='abundance'){

    if(class(data)=="data.frame" | class(data)=="matrix" ) data = list(Region_1 = data)

    if(class(data)== "list"){
      if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
      Ns = sapply(data, ncol)
      data_list = data
    }

  }

  if(datatype=='incidence_raw'){

    if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
    Ns = sapply(data, length)
    data_list = data

  }

  if(is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)

  for_each_region = function(data, region_name, N){

    #data
    if (datatype=='abundance') {

      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector

      # ref_gamma = iNEXT.3D:::Chat.Ind(data_gamma, n)
      # ref_alpha = iNEXT.3D:::Chat.Ind(data_alpha, n)
      #
      # ref_alpha_max = iNEXT.3D:::Chat.Ind(data_alpha, n * 2)
      # ref_gamma_max = iNEXT.3D:::Chat.Ind(data_gamma, n * 2)

      ref_gamma = iNEXT.3D:::Coverage(data_gamma, n, datatype = 'abundance')
      ref_alpha = iNEXT.3D:::Coverage(data_alpha, n, datatype = 'abundance')

      ref_alpha_max = iNEXT.3D:::Coverage(data_gamma, n*2, datatype = 'abundance')
      ref_gamma_max = iNEXT.3D:::Coverage(data_alpha, n*2, datatype = 'abundance')

      level = level[level<1]
      level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique

      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
    }

    if (datatype=='abundance') {
      ### 1. aL_table_gamma
      data_gamma_2d = long_to_wide(data_gamma)

      aL_table_gamma = create.aili(data_gamma_2d, row.tree = row.tree, col.tree = col.tree) %>%
        select(branch.abun, branch.length, tgroup)%>%
        filter(branch.abun>0)

      ### 2. aL_table_alpha
      plan(multisession)
      aL_table_alpha =  future_lapply(1:N, function(i){
        x = data[data[,i]>0,i]
        names(x) = rownames(data)[data[,i]>0]

        aL_table = create.aili(data = x%>%long_to_wide(),row.tree = row.tree, col.tree=col.tree )%>%
          select(branch.abun, branch.length, tgroup)%>%
          filter(branch.abun>0)
        return(aL_table)
      }, future.seed = TRUE)%>%do.call("rbind",.)


      get_phylogenetic_alpha_gamma <- function(aL_table_gamma, aL_table_alpha, n, m_gamma, m_alpha){
        ## 1. gamma
        qPDm <- iNEXTPD2:::PhD.m.est(ai = aL_table_gamma$branch.abun,
                                     Lis = aL_table_gamma$branch.length%>%as.matrix(),m = m_gamma,
                                     q = q,nt = n,cal = 'PD')
        gamma = qPDm %>% t %>%as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level=rep(level, 3), Coverage_real=rep(iNEXT.3D:::Coverage(data_gamma, m_gamma, datatype = 'abundance'), 3),
                 Size=rep(m_gamma, 3))%>%
          mutate(Method = ifelse(level>=ref_gamma, ifelse(level==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))

        ## 2. alpha
        qPDm = iNEXTPD2:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = aL_table_alpha$branch.length%>%as.matrix(),
                               m = m_alpha,q = q,nt = n,cal = 'PD')

        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level=rep(level, 3), Coverage_real=rep(iNEXT.3D:::Coverage(data_alpha, m_alpha, datatype = 'abundance'), 3), Size=rep(m_alpha, 3))%>%
          mutate(Method = ifelse(level>=ref_gamma, ifelse(level==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))
        res = list()
        res[['gamma']] = gamma
        res[['alpha']] = alpha
        return(res)
      }

      PD_results = get_phylogenetic_alpha_gamma(aL_table_gamma, aL_table_alpha,
                                                n = sum(data_gamma), m_gamma, m_alpha)
      gamma = PD_results$gamma
      alpha = PD_results$alpha
    }

    gamma = (gamma %>%
               mutate(Method = ifelse(level>=ref_gamma,
                                      ifelse(level==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
      set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

    if (max_alpha_coverage==T) {
      under_max_alpha = !((gamma$Order==0) & (gamma$level>ref_alpha_max))
    }else{
      under_max_alpha = gamma$level>0
    }
    gamma = gamma[under_max_alpha,]
    gamma$Order = as.numeric(gamma$Order)


    alpha = (alpha %>%
               mutate(Method = ifelse(level>=ref_alpha, ifelse(level==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
      set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

    alpha = alpha[under_max_alpha,]
    alpha$Order = as.numeric(alpha$Order)

    beta = alpha
    beta$Estimate = gamma$Estimate/alpha$Estimate
    beta[beta == "Observed"] = "Observed_alpha"
    beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate),
                                   Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))

    C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
    U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
    V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
    S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

    if(nboot>1){
      ### step0. get information from data (across samples)
      get_f0_hat = function(dat){
        dat <- as.matrix(dat)
        n <- sum(dat)
        f1 <- sum(dat == 1)
        f2 <- sum(dat == 2)
        f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
        return(f0)
      }
      data_to_pi <- function(data) {
        f1 <- sum(data == 1)
        f2 <- sum(data == 2)
        n <- sum(sapply(unique(data), function(x){x*sum(data == x)}))
        C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
        lambda <- ifelse(C.hat != 1,(1-C.hat)/sum(data/n*(1-data/n)^n),0)
        f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
        p.seen <- data/n*(1-lambda*(1-data/n)^n)
        # p.seen <- p.seen[p.seen>0]
        p.unseen <- (1-C.hat)/f0
        list(seen = p.seen,unseen = p.unseen)
      }

      f0_pool = get_f0_hat(data_gamma)
      pi <- lapply(data, function(x){data_to_pi(x)})

      g0_hat = sapply(1:ncol(data), function(k){
        k_net = data[,k]
        names(k_net) = rownames(data)
        n = sum(k_net)
        f1 = sum(k_net==1)
        f2 = sum(k_net==2)

        # aL = phyBranchAL_Abu(phylo = phy_tree, data = x, rootExtend = T, refT = reft)
        # aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL = create.aili(k_net%>%long_to_wide(), row.tree, col.tree)%>%
          select(branch.abun,branch.length)
        g1 = aL$branch.length[aL$branch.abun==1] %>% sum
        g2 = aL$branch.length[aL$branch.abun==2] %>% sum
        g0_hat = ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) ,
                         ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))))
        return(g0_hat)

      })

      plan(multisession)
      se = future_lapply(1:nboot, function(b){
        ## for each bootstrap samples
        ### 1. get population (random assign position from candidates)
        piLi_k = lapply(1:ncol(data), function(k){
          kth_net = data[,k]
          names(kth_net) = rownames(data)

          f1 <- sum(kth_net == 1)
          f2 <- sum(kth_net == 2)
          f0_k <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
          C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
          ## phy
          phy <- create.aili(kth_net%>%long_to_wide(),row.tree,col.tree)
          g1 <- sum(phy$branch.length[phy$branch.abun==1])
          g2 <- sum(phy$branch.length[phy$branch.abun==2])
          g0_k <- ifelse(g2>g1*f2/2/f1,(n-1)/n*g1^2/2/g2,(n-1)/n*g1*(f1-1)/2/(f2+1))/f0_k
          ##
          pi = data_to_pi(kth_net)
          ## pi_matrix: S1*S2 (42*98)
          pi_matrix = pi$seen%>%long_to_wide()
          ## S1*S2 -> B1*B2 (8736)
          seen_interaction_aili = create.aili(pi_matrix,row.tree,col.tree)
          # unseen_interaction_aili = data.frame(branch.abun = rep(pi$unseen,f0_k), branch.length = g0 / f0_k,
          #                                      tgroup = "Tip",interaction = "unseen")
          unseen_interaction_aili = data.frame(branch.abun = rep(0,f0_pool), branch.length = rep(0,f0_pool),
                                               tgroup = "Tip",interaction = "unseen")%>%
            mutate(interaction = paste0(interaction, row_number()))
          seen_unseen = rbind(seen_interaction_aili, unseen_interaction_aili)
          ## 8268 candidates out of 8977
          candidates = which(seen_unseen$branch.abun == 0)
          assigned_position = sample(candidates, size = min(f0_k, length(candidates)) , replace = F)
          seen_unseen$branch.abun[assigned_position] = (1-C.hat) / f0_k
          # seen_unseen <- seen_unseen%>%filter(branch.abun != 0)

          return(seen_unseen)
        })

        length_bt = lapply(piLi_k, function(pili){
          pili[, c("branch.length", "interaction")]%>%column_to_rownames("interaction")
        })%>%do.call("cbind",.)
        colnames(length_bt) = paste0("net_", 1:ncol(length_bt))

        p_bt = lapply(piLi_k, function(pili){
          pili[, c("branch.abun", "interaction")]%>%column_to_rownames("interaction")
        })%>%do.call("cbind",.)
        colnames(p_bt) = paste0("net_", 1:ncol(p_bt))

        p_bt = p_bt%>%cbind(tgroup = piLi_k[[1]]$tgroup)%>%
          filter_all(any_vars(. != 0))
        ### 2. get samples(x_bt)
        x_bt = sapply(1:(ncol(p_bt)-1), function(k){
          n = sum(data[,k])
          sapply(p_bt[,k], function(p){rbinom(1,size=n, prob = p)})

        })
        rownames(x_bt) = rownames(p_bt)

        L0_hat = sapply(1:ncol(data), function(k){

          compare = data.frame(interaction = names(data_gamma), raw = data_gamma)%>%
            inner_join(x_bt[,k]%>%as.data.frame()%>%rownames_to_column('interaction'),
                       by = 'interaction')%>%rename('sample'=".")
          know_unseen_interaction = which(and((compare$raw==0) ,(compare$sample!=0)))

          if(length(know_unseen_interaction) == 0){ used_length = 0}else{
            used_length = compare%>%filter(raw == 0, sample != 0 )%>%
              left_join(cbind(length_bt[,c('net_1')], interaction = rownames(length_bt))%>%
                          as.data.frame()%>%rename("length"="V1"),
                        by = 'interaction')%>%pull(length)%>%sum()
          }
          L0_hat_k = g0_hat[k] - used_length
          return(L0_hat_k)
        })

        for(k in 1:ncol(data)){
          length_bt[,k] = c(length_bt[,k]%>%.[1:(length(.)-f0_pool)],
                            rep(L0_hat[k]/ f0_pool, f0_pool))

        }
        length_bt_average = length_bt%>%rowMeans()%>%as.data.frame()%>%rownames_to_column("interaction")
        ###############
        ### 1. aL_table_gamma
        bootstrap_data_gamma = rowSums(x_bt)
        bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]

        bt_table_data_gamma = bootstrap_data_gamma%>%
          as.data.frame()%>%rownames_to_column("interaction")

        bt_aL_table_gamma = bt_table_data_gamma%>%
          left_join(length_bt_average , by = 'interaction')%>%
          set_colnames(c("interaction", 'branch.abun', 'branch.length'))

        ### 2. aL_table_alpha
        plan(multisession)

        bt_aL_table_alpha =  future_lapply(1:N, function(k){
          bootstrap_data_alpha = cbind(interaction = rownames(x_bt)[x_bt[,k]>0],
                                       branch.abun =  x_bt[x_bt[,k]>0,k])%>%as.data.frame()%>%
            mutate(branch.abun = as.integer(branch.abun))

          aL_table = bootstrap_data_alpha%>%
            left_join(length_bt_average , by = 'interaction')%>%
            set_colnames(c("interaction", 'branch.abun', 'branch.length'))

          return(aL_table)
        }, future.seed = TRUE)%>%do.call("rbind",.)


        bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
        bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]


        m_gamma_bt = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma, i, datatype='abundance'))
        m_alpha_bt = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha, i, datatype='abundance'))


        get_phylogenetic_alpha_gamma_bootstrap <- function(aL_table_gamma, aL_table_alpha,
                                                           n, m_gamma, m_alpha){
          ## 1. gamma
          gamma <- iNEXTPD2:::PhD.m.est(ai = aL_table_gamma$branch.abun,
                                        Lis = aL_table_gamma$branch.length%>%as.matrix(),m = m_gamma,
                                        q = q,nt = n,cal = 'PD')%>%
            t%>% as.vector()
          ## 2. alpha
          alpha =(PhD:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = aL_table_alpha$branch.length%>%as.matrix(),
                                  m = m_alpha,q = q,nt = n,cal = 'PD')/N)%>%
            t%>% as.vector()
          beta_obs = (iNEXT.3D:::PD.Tprofile(ai=aL_table_gamma$branch.abun,
                                             Lis=as.matrix(aL_table_gamma$branch.length),
                                             q=q, nt=n, cal="PD") /
                        (iNEXT.3D:::PD.Tprofile(ai=aL_table_alpha$branch.abun,
                                                Lis=as.matrix(aL_table_alpha$branch.length),
                                                q=q, nt=n, cal="PD") / N)) %>% unlist()%>%as.vector()

          res = list()
          res[['gamma']] = gamma
          res[['alpha']] = alpha
          res[['beta_obs']] = beta_obs
          return(res)
        }

        PD_bootstrap_results = get_phylogenetic_alpha_gamma_bootstrap(aL_table_gamma = bt_aL_table_gamma, aL_table_alpha = bt_aL_table_alpha,
                                                                      n = sum(bootstrap_data_gamma), m_gamma_bt, m_alpha_bt)
        gamma = PD_bootstrap_results$gamma
        alpha = PD_bootstrap_results$alpha
        beta_obs = PD_bootstrap_results$beta_obs

        ##
        gamma = gamma[under_max_alpha]
        alpha = alpha[under_max_alpha]
        beta = gamma/alpha

        gamma = c(gamma, rep(0, length(q)))
        alpha = c(alpha, rep(0, length(q)))
        beta = c(beta, beta_obs)

        order = rep(q, each=length(level) + 1)[under_max_alpha]

        beta = data.frame(Estimate=beta, order)

        C = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(1-order) - 1)/(N^(1-order)-1))))$Estimate
        U = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(order-1) - 1)/(N^(order-1)-1))))$Estimate
        V = (beta %>% mutate(Estimate = (Estimate-1)/(N-1)))$Estimate
        S = (beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1)))$Estimate

        beta = beta$Estimate

        diversity = cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
        return(diversity)

      }, future.seed = T)%>%
        abind(along=3) %>% apply(1:2, sd)


    } else {

      se = matrix(0, ncol = 7, nrow = nrow(gamma))
      colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
      se = as.data.frame(se)

    }
    se = as.data.frame(se)

    gamma = gamma %>% mutate(LCL = Estimate - tmp * se$gamma[1:(nrow(se)-length(q))],
                             UCL = Estimate + tmp * se$gamma[1:(nrow(se)-length(q))],
                             Region = region_name)

    alpha = alpha %>% mutate(LCL = Estimate - tmp * se$alpha[1:(nrow(se)-length(q))],
                             UCL = Estimate + tmp * se$alpha[1:(nrow(se)-length(q))],
                             Region = region_name)

    beta = beta %>% mutate(  LCL = Estimate - tmp * se$beta,
                             UCL = Estimate + tmp * se$beta,
                             Region = region_name)

    C = C %>% mutate(        LCL = Estimate - tmp * se$C,
                             UCL = Estimate + tmp * se$C,
                             Region = region_name)


    U = U %>% mutate(        LCL = Estimate - tmp * se$U,
                             UCL = Estimate + tmp * se$U,
                             Region = region_name)

    V = V %>% mutate(        LCL = Estimate - tmp * se$V,
                             UCL = Estimate + tmp * se$V,
                             Region = region_name)

    S = S %>% mutate(        LCL = Estimate - tmp * se$S,
                             UCL = Estimate + tmp * se$S,
                             Region = region_name)

    list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)
  }

  output = lapply(1:length(data), function(i) for_each_region(data = data_list[[i]],
                                                              region_name = region_names[i], N = Ns[i]))
  names(output) = region_names
  return(output)
}
