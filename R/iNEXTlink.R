#' Sample Completeness main function
#'
#' \code{NetSC} Estimation of Sample Completeness with order q
#'
#' @param x a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr
#' @param q a integer vector for the order of Hill number\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).\cr
#' @param nboot an integer specifying the number of bootstrap replications, default is 30.\cr
#' @param conf  positive number < 1 specifying the level of confidence interval, default is 0.95.\cr\cr
#' @return a matrix of estimated sample completeness with order q: \cr\cr
#'
#' @examples
#' data(Norfolk)
#' output = NetSC(Norfolk)
#' ggNetSC(output)
#'
#' @references
#' Chao,A.,Y.Kubota,D.Zelený,C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages.
#' @export
NetSC <- function(x, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30, conf = 0.95){
  data_long <- lapply(x, function(tab){
    as.matrix(tab)%>%c()}
  )
  res = iNEXT4steps::SC(x = data_long, q = q, datatype = datatype, nboot = nboot, conf = conf)
}

# ggSC -------------------------------------------------------------------
#' ggplot for Sample Completeness
#'
#' \code{ggNetSC} The figure for estimation of Sample Completeness with order q
#'
#' @param output a table generated from SC function
#' @return a figure of estimated sample completeness with order q
#'
#' @examples
#' data(Norfolk)
#' output = NetSC(Norfolk)
#' ggNetSC(output)
#' @references
#' Chao,A.,Y.Kubota,D.Zelený,C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages.
#' @export
ggNetSC <- function(output){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(output, aes(x = Order.q, y = Estimate.SC, colour = Assemblage)) +
    geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
    geom_ribbon(aes(ymin = SC.LCL, ymax = SC.UCL, fill = Assemblage),
                alpha = 0.2, linetype = 0) + scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Sample completeness") + theme(text = element_text(size = 18)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"), legend.title = element_blank())
}

#' Interpolation (rarefaction) and extrapolation of Chao et al.’s (2021) network diversity and mean network diversity
#' Function \code{iNEXTPD} computes network diversity estimates for rarefied samples and extrapolated samples
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
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param endpoint a positive integer specifying the endpoint for the rarefaction and extrapolation range.
#' If \code{NULL}, then \code{endpoint} = double of the reference sample size in each assemblage. It is ignored if \code{size} is given.
#' @param knots a positive integer specifying the number of equally-spaced knots between 1 and the \code{endpoint}. Default is 40.
#' @param size a sequence of positive integers specifying the sample sizes for which PD or meanPD estimates will be calculated.
#' If \code{NULL}, then estimates will be calculated for those sample sizes determined by the specified/default \code{endpoint}
#' and \code{knots}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import ape ggplot2 dplyr tidytree stats
#' @importFrom phyclust get.rooted.tree.height
#' @return
#'
#'
#' \itemize{
#'  \item{\code{$DataInfo}: A dataframe summarizing data information}
#'  \item{\code{$iNextEst}: coverage-based diversity estimates along with confidence intervals}
#'  (if \code{nboot > 0}) for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#'  \item{\code{$AsyEst}: for
#' showing asymptotic diversity estimates along with related statistics.}
#' }
#' @examples
#' \dontrun{
#' # Datatype: abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- iNEXT_ND(data = data, q = c(0, 1, 2), nboot = 30)
#' out
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
#' \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export

iNEXT_ND <- function(x, q = c(0,1,2), datatype = "abundance", size = NULL,
                     endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = 30){
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
  ## calculate
  res = list()
  ## 1. datainfo
  datainfo = NDinfo(x = x, datatype = datatype)

  ## 2. iNterpolation/ Extrapolation
  data_long <- lapply(x, function(tab){
    as.matrix(tab)%>%c()}
    ## 拉長
  )
  INEXT_est <- iNEXT::iNEXT(data_long,q = q,conf = conf,nboot = nboot, knots = knots, se = se, endpoint = endpoint, size = size)

  ## 3. empirical and asymptotic diversity
  # abund <- lapply(dat, function(x) {
  #   tmp <- c(as.matrix(x))
  #   tmp = x
  #   AsyND(data = tmp, B = 10, q, conf = 0.95)
  # })
  # mle = lapply(dat, function(x){
  #   tmp <- c(as.matrix(x))
  #   ObsND(data = tmp, B = 10, q, conf = 0.95)
  # })
  #
  # plot.m <- c()
  # for ( i in 1:length(dat)){
  #   temp <- abund[[i]][abund[[i]]$Target == "Diversity",-c(2,4)]
  #   mle1 <- mle[[i]][mle[[i]]$Target == "Diversity",-c(2,4)]
  #   names(mle1)[2] <- "Estimate"
  #
  #   out_est = data.frame(temp, method = "Estimate", Region = names(dat)[i])
  #   out_emp = data.frame(mle1, method = "Empirical", Region = names(dat)[i])
  #   plot.m = rbind(plot.m,out_est,out_emp)
  # }
  # plot.m$Region <- factor(plot.m$Region,levels = levels(factor(names(dat))))
  res[[1]] = datainfo
  res[[2]] = INEXT_est$iNextEst
  res[[3]] = INEXT_est$AsyEst
  names(res) = c("DataInfo", "iNextEst", "AsyEst")
  class(res) <- c("iNEXT")
  return(res)
}

#' Exhibit basic data information
#'
#' \code{DataInfo}: exhibits basic data information
#'
#' @param x a vector/matrix/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence"},
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
#' @examples
#' \dontrun{
#' data(Norfolk)
#' DataInfo(Norfolk, datatype="abundance")
#' }

#' @export
NDinfo <- function(x, datatype = "abundance"){
  table <- lapply(x, function(y){datainf(data = y, datatype = datatype)})%>%do.call(rbind,.)
  rownames(table) <- names(x)
  table = rownames_to_column(table, var = "Assemblages")
  return(table)
}
#' Exhibit basic data information
#'
#' \code{DataInfo}: exhibits basic data information
#'
#' @param x a vector/matrix/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence"},
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
#' @examples
#' \dontrun{
#' data(Norfolk)
#' PNDinfo(Norfolk, datatype="abundance", row.tree = rowtree, col.tree = coltree)
#' }

#' @export
PNDinfo <- function(data, datatype = "abundance", row.tree = NULL,col.tree = NULL){
  table <- lapply(data, function(y){datainfphy(data = y, datatype = datatype,
                                            row.tree = row.tree,col.tree = col.tree)})%>%
    do.call(rbind,.)
  rownames(table) <- names(data)
  table = rownames_to_column(table, var = "Assemblages")
  return(table)
}

#' Interpolation (rarefaction) and extrapolation of Chao et al.’s (2021) network diversity and mean network diversity

#' Function \code{iNEXTPD} computes network diversity estimates for rarefied samples and extrapolated samples
#' along with confidence intervals and related coverage estimates based on Chao et al.’s (2021) network
#' diversity (ND)
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
#' Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
#' data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
#' Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
#' and thus N matrices, users must first merge the N matrices by species identity to obtain a large
#' merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
#' data. The row (species) names of data must match the species names in the phylogenetic tree and
#' thus cannot be missing.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param endpoint a positive integer specifying the endpoint for the rarefaction and extrapolation range.
#' If \code{NULL}, then \code{endpoint} = double of the reference sample size in each assemblage. It is ignored if \code{size} is given.
#' @param knots a positive integer specifying the number of equally-spaced knots between 1 and the \code{endpoint}. Default is 40.
#' @param size a sequence of positive integers specifying the sample sizes for which PD or meanPD estimates will be calculated.
#' If \code{NULL}, then estimates will be calculated for those sample sizes determined by the specified/default \code{endpoint}
#' and \code{knots}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import tidytree
#' @importFrom stats rbinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom phyclust get.rooted.tree.height
#' @return
#' \itemize{
#'  \item{\code{$DataInfo}: A dataframe summarizing data information}
#'  \item{\code{$iNextEst}: coverage-based diversity estimates along with confidence intervals}
#'  (if \code{nboot > 0}) for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#'  \item{\code{$AsyEst}: for
#' showing asymptotic diversity estimates along with related statistics.}
#' }


#' @examples
#' \dontrun{
#' # Datatype: abundance data with phylogenetic tree
#' data(puerto.rico)
#' data <- puerto.rico$data
#' rowtree <- puerto.rico$rowtree
#' coltree <- puerto.rico$coltree
#' out <- iNEXT_PND(data = data, row.tree = rowtree, col.tree = coltree, q = c(0, 1, 2))
#' out
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
#' \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export
iNEXT_PND <- function(x, row.tree = NULL,col.tree = NULL,q = c(0,1,2), datatype = "abundance",
                      size = NULL, endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = 30){
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
  ## calculate
  res = list()
  ## 1. datainfo
  datainfo = PNDinfo(data = x, datatype = datatype, row.tree = row.tree,col.tree = col.tree)

  ## 2. iNterpolation/ Extrapolation
  data_long <- lapply(x, function(tab){
    as.matrix(tab)%>%c()}
    ## 拉長
  )

  NetiNE <- get.netphydiv_iNE(data = x, q = q,B = nboot,row.tree = row.tree,col.tree = col.tree,conf = conf, knots = knots)
  ## 3. empirical and asymptotic diversity
  NetDiv <- get.netphydiv(data = x,q = q,B = nboot,row.tree = row.tree,col.tree = col.tree,conf = conf)

  res[[1]] = datainfo
  res[[2]] = NetiNE
  res[[3]] = NetDiv
  names(res) = c("DataInfo", "iNextEst", "AsyEst")
  class(res) <- c("iNEXT")
  return(res)
}
#' Interpolation (rarefaction) and extrapolation of Chao et al.’s (2021) network diversity and mean network diversity

#' Function \code{iNEXTPD} computes network diversity estimates for rarefied samples and extrapolated samples
#' along with confidence intervals and related coverage estimates based on Chao et al.’s (2021) beta diversity
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
#' Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
#' data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
#' Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
#' and thus N matrices, users must first merge the N matrices by species identity to obtain a large
#' merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
#' data. The row (species) names of data must match the species names in the phylogenetic tree and
#' thus cannot be missing.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.

#' @param knots a positive integer specifying the number of equally-spaced knots between 1 and the \code{endpoint}. Default is 40.
#' @param size a sequence of positive integers specifying the sample sizes for which PD or meanPD estimates will be calculated.
#' If \code{NULL}, then estimates will be calculated for those sample sizes determined by the specified/default \code{endpoint}
#' and \code{knots}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import tidytree
#' @importFrom stats rbinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom phyclust get.rooted.tree.height
#' @return
#' \itemize{
#'  \item{\code{$DataInfo}: A dataframe summarizing data information}
#'  \item{\code{$iNextEst}: coverage-based diversity estimates along with confidence intervals}
#'  (if \code{nboot > 0}) for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#'  \item{\code{$AsyEst}: for
#' showing asymptotic diversity estimates along with related statistics.}
#' }


#' @examples
#' \dontrun{
#' # example
#' data(Norfolk)
#' iNEXT_beta_link(Norfolk, coverage_expected, data_type=c('abundance', 'incidence_raw'), q = c(0, 1, 2),
#' level=c('taxonomic', 'phylogenetic', 'functional'), nboot = 20, conf = 0.95, max_alpha_coverage=F,
#'  by=c('coverage', 'size'),phy_tree=NULL, reftime = NULL)
#'
#'  dissimilarity1 = iNEXT_beta_link(Norfolk, coverage_expected = seq(0.5,1,0.05),
#'  data_type='abundance', q = c(0, 1, 2),level='taxonomic',
#'  nboot = 20, conf = 0.95, max_alpha_coverage=F, by='coverage')

#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
#' \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export

# row.tree = rowtree
# col.tree = coltree
# x = puerto.rico
# coverage_expected = seq(0.5,1,0.05)
# level = 'phylogenetic'; data_type = "abundance";conf=0.95;by = 'coverage';nboot=20
# level = 'taxonomic'

iNEXT_beta_link = function(x, coverage_expected = seq(0.5, 1, 0.5), data_type=c('abundance', 'incidence_raw'), q = c(0, 1, 2), level=c('taxonomic', 'phylogenetic', 'functional'),
                           nboot = 20, conf = 0.95, max_alpha_coverage=F, by=c('coverage', 'size'),
                           row.tree = NULL,col.tree = NULL){

  combined = ready4beta(x)
  if(level == 'taxonomic'){
    # dissimilarity <- iNEXT_beta(x = combined, coverage_expected = coverage_expected, data_type = data_type, level = 'taxonomic',
    dissimilarity <- iNEXT_beta(x = combined, coverage_expected = coverage_expected, data_type = data_type, level = 'taxonomic',
                                nboot = nboot, conf = conf, max_alpha_coverage = max_alpha_coverage, by = by)
  }
  else if(level == 'phylogenetic'){

    dissimilarity = iNEXT_link_phybeta(x = combined, coverage_expected =coverage_expected, "abundance", level = 'phylogenetic',
                                       row.tree = rowtree,col.tree = coltree,
                                       nboot = 0, by = 'coverage')
  }

  return(dissimilarity)

}



# iNEXT_beta_link(x = puerto.rico, coverage_expected = seq(0.5,1,0.05), "abundance", level = 'phylogenetic',
#                 row.tree = rowtree,col.tree = coltree,
#                    nboot = 0, by = 'coverage')
#
#
# str(dissimilarity)

iNEXT_link_phybeta <- function(x, coverage_expected, data_type=c('abundance', 'incidence_raw'), q = c(0, 1, 2), level=c('taxonomic', 'phylogenetic', 'functional'),
                               nboot = 20, conf = 0.95, max_alpha_coverage=F, by=c('coverage', 'size'),
                               row.tree = NULL,col.tree = NULL){
  if(data_type=='abundance'){

    if( class(x)=="data.frame" | class(x)=="matrix" ) x = list(Region_1 = x)

    if(class(x)== "list"){
      if(is.null(names(x))) region_names = paste0("Region_", 1:length(x)) else region_names = names(x)
      Ns = sapply(x, ncol)
      data_list = x
    }

  }

  if(data_type=='incidence_raw'){

    if(is.null(names(x))) region_names = paste0("Region_", 1:length(x)) else region_names = names(x)
    Ns = sapply(x, length)
    data_list = x

  }

  if(is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)

  if(level=='phylogenetic'){

    if (data_type=='abundance') {
      pool.data = do.call(cbind, data_list) %>% rowSums
    }

    else if (data_type=='incidence_raw') {
      pool.data = do.call(cbind, data_list[[1]]) %>% rowSums

    }
    pool.name = names(pool.data[pool.data>0])
    # tip = phy_tree$tip.label[-match(pool.name, phy_tree$tip.label)]
    # mytree = drop.tip(phy_tree, tip)
    # H_max = get.rooted.tree.height(mytree)
    # if(is.null(reftime)) { reft = H_max
    # } else if (reftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
    # } else { reft = reftime }
  }

  for_each_region = function(data, region_name, N){

    #data
    if (data_type=='abundance') {

      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector

      ref_gamma = iNEXT:::Chat.Ind(data_gamma, n)
      if (by=='size') ref_alpha = ref_gamma
      if (by=='coverage') ref_alpha = iNEXT:::Chat.Ind(data_alpha, n)
      ref_alpha_max = iNEXT:::Chat.Ind(data_alpha, n*2)

      coverage_expected = c(coverage_expected, ref_gamma, ref_alpha, ref_alpha_max) %>% sort %>% unique
      # coverage_expected = coverage_expected[coverage_expected<1]

      m_gamma = sapply(coverage_expected, function(i) coverage_to_size(data_gamma, i, data_type='abundance'))
      if (by=='size') m_alpha = m_gamma
      if (by=='coverage') m_alpha = sapply(coverage_expected, function(i) coverage_to_size(data_alpha, i, data_type='abundance'))

    }
    if (data_type=='incidence_raw') {

      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)

      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))

      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)

      # data_gamma_freq = data_gamma_freq[data_gamma_freq>0]
      # data_alpha_freq = data_alpha_freq[data_alpha_freq>0]

      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

      ref_gamma = iNEXT:::Chat.Sam(data_gamma_freq, n)
      if (by=='size') ref_alpha = ref_gamma
      if (by=='coverage') ref_alpha = iNEXT:::Chat.Sam(data_alpha_freq, n)
      ref_alpha_max = iNEXT:::Chat.Sam(data_alpha_freq, n*2)

      coverage_expected = c(coverage_expected, ref_gamma, ref_alpha, ref_alpha_max) %>% sort %>% unique
      # coverage_expected = coverage_expected[coverage_expected<1]

      m_gamma = sapply(coverage_expected, function(i) coverage_to_size(data_gamma_freq, i, data_type='incidence_freq'))
      if (by=='size') m_alpha = m_gamma
      if (by=='coverage') m_alpha = sapply(coverage_expected, function(i) coverage_to_size(data_alpha_freq, i, data_type='incidence_raw'))

    }

    if (level=='phylogenetic') {

      if (data_type=='abundance') {
        data_gamma_2d = long_to_wide(data_gamma)
        aL_table_gamma = create.aili(data_gamma_2d, row.tree = row.tree, col.tree = col.tree) %>%
          select(branch.abun, branch.length, tgroup)

        n <- sum(data_gamma)
        cal = 'PD'

        qPDm <- iNEXTPD2:::PhD.m.est(ai = aL_table_gamma$branch.abun,Lis = aL_table_gamma$branch.length%>%as.matrix(),m = m_gamma,
                                     q = c(0,1,2),nt = n,cal = cal)
        # %>% as.numeric()

        # qPDm <- iNEXTPD2:::PhD.m.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,m = m_gamma,
        #                              q = c(0,1,2),nt = n,cal = cal) %>% as.numeric()

        # iNEXTPD2:::PhD.m.est(ai = aL_table_gamma, lis, m = m_gamma , q = c(0,1,2), nt = n, cal)
        gamma = qPDm %>% t %>%
          as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3),
                 Size=rep(m_gamma, 3))%>%
          mutate(Method = ifelse(Coverage_expected>=ref_gamma, ifelse(Coverage_expected==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))
        # gamma = iNEXTPD2:::PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
        #                              m=m[[i]],q=q,nt = n,cal = cal) %>% t %>%
        #   as.data.frame %>%
        #   set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
        #   mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3),
        #          Size=rep(m_gamma, 3))
        aL_table_alpha = c()

        for (i in 1:N){

          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]

          # aL = phyBranchAL_Abu(phylo = phy_tree, data = x, rootExtend = T, refT = reft)
          # aL$treeNabu$branch.length = aL$BLbyT[,1]
          # aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          #
          # aL_table_alpha = rbind(aL_table_alpha, aL_table)
          aL_table = create.aili(data = x%>%long_to_wide(),row.tree = row.tree, col.tree=col.tree )%>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }

        # qPDm = PhD:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='abundance', nt=n)
        qPDm = PhD:::PhD.m.est(ai = aL_table_alpha$branch.abun,Lis = aL_table_alpha$branch.length%>%as.matrix(),m = m_alpha,
        q = c(0,1,2),nt = n,cal = cal)
        # qPDm = iNEXTPD2:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='abundance', nt=n)
        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))%>%
          mutate(Method = ifelse(Coverage_expected>=ref_gamma, ifelse(Coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))

      }

      if (data_type=='incidence_raw') {

        aL = phyBranchAL_Inc(phylo=phy_tree, data=as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

        gamma = PhD:::PhD.m.est(aL = aL_table_gamma, m = m_gamma, Q = c(0,1,2), datatype = 'incidence_raw', nt = n) %>% t %>% as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Sam(data_gamma_freq, m_gamma), 3), Size=rep(m_gamma, 3))

        aL_table_alpha = c()

        for (i in 1:N){

          x = data[[i]]

          aL = phyBranchAL_Inc(phylo = phy_tree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }

        alpha = (PhD:::PhD.m.est(aL = aL_table_alpha, m = m_alpha, Q = c(0,1,2), datatype = 'incidence_raw', nt = n)/N) %>% t %>% as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))


      }

      gamma = (gamma %>%
                 mutate(Method = ifelse(Coverage_expected>=ref_gamma,
                                        ifelse(Coverage_expected==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      if (max_alpha_coverage==T) under_max_alpha = !((gamma$Order==0) & (gamma$Coverage_expected>ref_alpha_max)) else under_max_alpha = gamma$Coverage_expected>0
      gamma = gamma[under_max_alpha,]
      gamma$Order = as.numeric(gamma$Order)


      alpha = (alpha %>%
                 mutate(Method = ifelse(Coverage_expected>=ref_alpha, ifelse(Coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      alpha = alpha[under_max_alpha,]
      alpha$Order = as.numeric(alpha$Order)

      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate
      beta[beta == "Observed"] = "Observed_alpha"
      beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate),
                                     Order = q, Method = "Observed", Coverage_expected = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))

      C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

      if(nboot>1){
        se = future_lapply(1:nboot, function(i){

          if (data_type=='abundance') {
            # tree_bt = phy_tree
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol=ncol(data))

            if ( nrow(p_bt) > nrow(data) & sum(unseen_p)>0 ){

              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)

              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample

              rownames(x_bt) = rownames(p_bt)

              if ( sum(x_bt[-(1:nrow(data)),])>0 ){
                # g0_hat = apply(data, 2, function(x){
                #
                #   n = sum(x)
                #   f1 = sum(x==1)
                #   f2 = sum(x==2)
                #
                #   # aL = phyBranchAL_Abu(phylo = phy_tree, data = x, rootExtend = T, refT = reft)
                #   # aL$treeNabu$branch.length = aL$BLbyT[,1]
                #   aL = create.aili(x%>%long_to_wide(), row.tree, col.tree)
                #
                #   aL = aL$treeNabu %>% select(branch.abun,branch.length)
                #   g1 = aL$branch.length[aL$branch.abun==1] %>% sum
                #   g2 = aL$branch.length[aL$branch.abun==2] %>% sum
                #   g0_hat = ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                #   g0_hat
                #
                # })
                g0_hat = sapply(1:ncol(data), function(j){
                  x = data[,j]
                  names(x) = rownames(data)
                  n = sum(x)
                  f1 = sum(x==1)
                  f2 = sum(x==2)

                  # aL = phyBranchAL_Abu(phylo = phy_tree, data = x, rootExtend = T, refT = reft)
                  # aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = create.aili(x%>%long_to_wide(), row.tree, col.tree)%>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun==1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun==2] %>% sum
                  g0_hat = ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  g0_hat

                })

                te = (x_bt[1:nrow(data),]*(data==0))>0
                used_length = sapply(1:ncol(data), function(i) {

                  if (sum(te[,i])==0) return(0) else {

                    aili = create.aili(x_bt[1:nrow(data),i]%>%long_to_wide(), row.tree, col.tree)
                    aili%>%
                      subset(interaction %in% (names(which(te[,i]==TRUE))%>%str_replace("_", "-")) ) %>%
                      select(branch.length) %>% sum

                    # phyBranchAL_Abu(phylo = phy_tree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                    #   subset(label %in% names(which(te[,i]==TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                g0_hat = g0_hat-used_length
                g0_hat[g0_hat<0] = 0

                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol=ncol(x_bt), byrow=T)

                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i]>0)>0) (g0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow=T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample)==0)] = 0

                # for (i in 1:length(L0_hat)){
                #
                #   tip = list(edge=matrix(c(2,1),1,2),
                #              tip.label=unseen_name[i],
                #              edge.length=L0_hat[i],
                #              Nnode=1)
                #   class(tip) = "phylo"
                #
                #   tree_bt = tree_bt + tip
                #
                # }

              } else {

                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]

              }

            } else {

              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)

            }

            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]

            m_gamma = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_gamma, i, data_type='abundance'))
            m_alpha = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_alpha, i, data_type='abundance'))

            # aL = phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            # aL$treeNabu$branch.length = aL$BLbyT[,1]
            # aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            ### 不確定如何處理unseen species in bootstraped sample, (long to wide)
            aL_table_gamma = create.aili(data = bootstrap_data_gamma%>%long_to_wide(), row.tree = row.tree, col.tree)%>%
              select(branch.abun, branch.length, tgroup)

            gamma = iNEXTPD2:::PhD.m.est(ai = aL_table_gamma$branch.abun,Lis = aL_table_gamma$branch.length%>%as.matrix(),m = m_gamma,
                                 q = c(0,1,2),nt = n,cal = cal) %>% t %>%as.data.frame%>%
              set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
              mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3),
                     Size=rep(m_gamma, 3))

            # gamma = as.vector(PhD:::PhD.m.est(aL=aL_table_gamma, m=m_gamma, Q=c(0,1,2), datatype='abundance', nt=n) %>% t)


            aL_table_alpha = c()

            # for (i in 1:N){
            #
            #   # x = x_bt[x_bt[,i]>0,i]
            #   # names(x) = rownames(p_bt)[x_bt[,i]>0]
            #
            #   x = x_bt[,i]
            #   names(x) = rownames(p_bt)
            #   x = x[x_bt[,i]>0]
            #
            #   aL = phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
            #   aL$treeNabu$branch.length = aL$BLbyT[,1]
            #   aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            #
            #   aL_table_alpha = rbind(aL_table_alpha, aL_table)
            #
            # }
            for (i in 1:N){

              x = x_bt[x_bt[,i]>0,i]
              names(x) = rownames(x_bt)[x_bt[,i]>0]
              aL_table = create.aili(data = x%>%long_to_wide(),row.tree = row.tree, col.tree=col.tree)%>% select(branch.abun, branch.length, tgroup)
              aL_table_alpha = rbind(aL_table_alpha, aL_table)

            }
            qPDm = PhD:::PhD.m.est(ai = aL_table_alpha$branch.abun,Lis = aL_table_alpha$branch.length%>%as.matrix(),m = m_alpha,
                                   q = c(0,1,2),nt = n,cal = cal)
            qPDm = qPDm/N
            alpha = qPDm %>% t %>% as.data.frame %>%
              set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
              mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))
            # alpha = as.vector((PhD:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='abundance', nt=n)/N) %>% t)

          }

          if (data_type=='incidence_raw') {

            tree_bt = phy_tree

            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)

            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){

              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)

              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){

                R0_hat = sapply(data, function(x){

                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)

                  aL = phyBranchAL_Inc(phylo = phy_tree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)

                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun==1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun==2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  R0_hat

                })

                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums)==0))>0
                used_length = sapply(1:N, function(i) {

                  if (sum(te[,i])==0) return(0) else {

                    phyBranchAL_Inc(phylo = phy_tree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i]==TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                R0_hat = R0_hat-used_length
                R0_hat[R0_hat<0] = 0

                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol=N, byrow=T)

                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i]>0)>0) (R0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow=T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample)==0)] = 0

                for (i in 1:length(L0_hat)){

                  tip = list(edge=matrix(c(2,1),1,2),
                             tip.label=unseen_name[i],
                             edge.length=L0_hat[i],
                             Nnode=1)
                  class(tip) = "phylo"

                  tree_bt = tree_bt + tip

                }

              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])

            } else {

              p_bt = p_bt[1:nrow(data),]
              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n=nT, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

            }
            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))

            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)

            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq>0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq>0]

            m_gamma = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_gamma_freq, i, data_type='incidence'))
            m_alpha = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_alpha_freq, i, data_type='incidence'))

            aL = phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            gamma = as.vector(PhD:::PhD.m.est(aL=aL_table_gamma, m=m_gamma, Q=c(0,1,2), datatype='incidence_raw', nt=n) %>% t)

            aL_table_alpha = c()

            for (i in 1:N){
              x = raw[[i]]
              aL = phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
            }
            alpha = as.vector((PhD:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='abundance', nt=n)/N) %>% t)
            alpha = as.vector((PhD:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='incidence_raw', nt=n)/N) %>% t)

          }

          alpha$Order = alpha$Order%>%as.numeric()
          beta = alpha
          beta$Estimate = gamma$Estimate/alpha$Estimate

          C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
          U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
          V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
          S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))


          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along=3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

    }

    if (level=='functional') {

      FD_by_tau = function(data, distance_matrix, tau, coverage_expected, data_type, by) {

        if (data_type=='abundance') {

          zik = data
          zik = zik[rowSums(data)>0,]

          dij = distance_matrix
          dij = dij[rowSums(data)>0, rowSums(data)>0]

          dij[which(dij>tau, arr.ind = T)] = tau
          aik = (1 - dij/tau) %*% as.matrix(zik)
          positive_id = rowSums(aik)>0



          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          # gamma_a = ifelse(gamma_a<1, 1, round(gamma_a))

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, c(0,1,2), n) %>% as.vector


          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          # alpha_a = ifelse(alpha_a<1, 1, round(alpha_a))

          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))

          if (by=='size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, c(0,1,2), n)/N) %>% as.vector
          if (by=='coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, c(0,1,2), n)/N) %>% as.vector

        }

        if (data_type=='incidence_raw') {

          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D

          gamma_Y = data_gamma_freq[-1]

          dij = distance_matrix
          dij = dij[gamma_Y>0, gamma_Y>0]
          gamma_Y = gamma_Y[gamma_Y>0]

          dij[which(dij>tau, arr.ind = T)] = tau
          gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          gamma_a[gamma_a>n] = n

          gamma_v = gamma_Y/gamma_a

          # gamma_a = ifelse(gamma_a<1, 1, round(gamma_a))
          # gamma_a[gamma_a>n] = n

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, c(0,1,2), n) %>% as.vector


          alpha_Y = data_2D[-1,]

          dij = distance_matrix
          dij = dij[rowSums(data_2D[-1,])>0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]

          dij[which(dij>tau, arr.ind = T)] = tau
          alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)

          # alpha_a = ifelse(alpha_a<1, 1, round(alpha_a))
          alpha_a[alpha_a>n] = n
          alpha_a = as.vector(alpha_a)

          # alpha_v = rep(gamma_v, N)
          # alpha_v = alpha_v[alpha_a>0]
          alpha_v = as.vector(as.matrix(alpha_Y))/alpha_a
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))

          if (by=='size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, c(0,1,2), n)/N) %>% as.vector
          if (by=='coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, c(0,1,2), n)/N) %>% as.vector

        }

        return(data.frame(gamma,alpha))

      }

      if (tau_type=='single'){

        if (data_type=='abundance') {

          if (by=='size') {

            output = FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='size')
            gamma = output$gamma
            alpha = output$alpha

          }

          if (by=='coverage') {

            output = FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='coverage')
            gamma = output$gamma
            alpha = output$alpha

          }

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected>=ref_gamma, ifelse(coverage_expected==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected>=ref_alpha, ifelse(coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

        }

        if (data_type=='incidence_raw') {

          if (by=='size') {

            output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='size')
            gamma = output$gamma
            alpha = output$alpha

          }

          if (by=='coverage') {

            output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='coverage')
            gamma = output$gamma
            alpha = output$alpha

          }

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected>=ref_gamma, ifelse(coverage_expected==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_gamma_freq, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected>=ref_alpha, ifelse(coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

        }

      }

      if (tau_type=='AUC'){

        cut = seq(0.00000001, 1, length.out = cut_number)
        width = diff(cut)

        if (data_type=='abundance') {

          if (by=='size') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='size')

            })

          }

          if (by=='coverage') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='coverage')

            })

          }
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected>=ref_gamma, ifelse(coverage_expected==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected>=ref_alpha, ifelse(coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = data.frame(coverage_expected, beta) %>%
            mutate(Method = ifelse(coverage_expected>=ref_alpha, ifelse(coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))


        }

        if (data_type=='incidence_raw') {

          if (by=='size') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='size')

            })

          }

          if (by=='coverage') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='coverage')

            })

          }
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected>=ref_gamma, ifelse(coverage_expected==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_gamma_freq, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected>=ref_alpha, ifelse(coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = data.frame(coverage_expected, beta) %>%
            mutate(Method = ifelse(coverage_expected>=ref_alpha, ifelse(coverage_expected==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))

        }

      }

      gamma = gamma[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      if (max_alpha_coverage==T) under_max_alpha = !((gamma$Order==0) & (gamma$Coverage_expected>ref_alpha_max)) else under_max_alpha = gamma$Coverage_expected>0
      gamma = gamma[under_max_alpha,]



      alpha = alpha[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      alpha = alpha[under_max_alpha,]

      beta = beta[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      beta = beta[under_max_alpha,]

      C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

      if(nboot>1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"coverage_expected","N",'under_max_alpha',
        #                     'data_type', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess, workers=7)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (data_type=='abundance') {

            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)

            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), distance_matrix, f0_hat, 'abundance')

            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))

            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector

            if (tau_type=='single'){

              if (by=='size') {

                output = FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='size')
                gamma = output$gamma
                alpha = output$alpha

              }

              if (by=='coverage') {

                output = FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='coverage')
                gamma = output$gamma
                alpha = output$alpha

              }

              beta=gamma/alpha

            }

            if (tau_type=='AUC'){


              if (by=='size') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='size')

                })

              }

              if (by=='coverage') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='coverage')

                })

              }
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

            }

          }

          if (data_type=='incidence_raw') {

            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])

            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), distance_matrix, f0_hat, 'incidence_freq')

            # p_bt = p_bt[rowSums(p_bt)>0,]

            raw = lapply(1:ncol(p_bt), function(j){

              lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

            })

            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))

            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)

            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt>0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt>0]

            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

            if (tau_type=='single'){

              if (by=='size') {

                output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='size')
                gamma = output$gamma
                alpha = output$alpha

              }

              if (by=='coverage') {

                output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='coverage')
                gamma = output$gamma
                alpha = output$alpha

              }

              beta = gamma/alpha
            }

            if (tau_type=='AUC'){


              if (by=='size') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='size')

                })

              }

              if (by=='coverage') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='coverage')

                })

              }
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

            }

          }

          gamma = gamma[under_max_alpha]
          alpha = alpha[under_max_alpha]
          beta = beta[under_max_alpha]

          order = rep(c(0,1,2), each=length(coverage_expected))[under_max_alpha]

          beta = data.frame(Estimate=beta, order)

          C = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(1-order) - 1)/(N^(1-order)-1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(order-1) - 1)/(N^(order-1)-1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate-1)/(N-1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1)))$Estimate

          beta = beta$Estimate

          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along=3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

    }


    se = as.data.frame(se)
    # gamma <- gamma%>%rename("LCL"="qD.LCL", "UCL"="qD.UCL")
    # alpha <- alpha%>%rename("LCL"="qD.LCL", "UCL"="qD.UCL")
    # beta <- beta%>%rename("LCL"="qD.LCL", "UCL"="qD.UCL")

    gamma = gamma %>% mutate(LCL = Estimate - 1.96 * se$gamma,
                             UCL = Estimate + 1.96 * se$gamma,
                             Region = region_name)

    alpha = alpha %>% mutate(LCL = Estimate - 1.96 * se$alpha,
                             UCL = Estimate + 1.96 * se$alpha,
                             Region = region_name)

    beta = beta %>% mutate(  LCL = Estimate - 1.96 * se$beta,
                             UCL = Estimate + 1.96 * se$beta,
                             Region = region_name)

    C = C %>% mutate(        LCL = Estimate - 1.96 * se$C,
                             UCL = Estimate + 1.96 * se$C,
                             Region = region_name)


    U = U %>% mutate(        LCL = Estimate - 1.96 * se$U,
                             UCL = Estimate + 1.96 * se$U,
                             Region = region_name)

    V = V %>% mutate(        LCL = Estimate - 1.96 * se$V,
                             UCL = Estimate + 1.96 * se$V,
                             Region = region_name)

    S = S %>% mutate(        LCL = Estimate - 1.96 * se$S,
                             UCL = Estimate + 1.96 * se$S,
                             Region = region_name)

    list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)

  }

  output = lapply(1:length(x), function(i) for_each_region(data = data_list[[i]],
                                                           region_name = region_names[i], N = Ns[i]))
  names(output) = region_names
  return(output)
}
# ggiNEXT -------------------------------------------------------------------
#' ggplot2 extension for an iNEXT object
#'
#' \code{ggiNEXT_ND}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
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
#' #' data(Norfolk)
#' out1 <- iNEXT_ND(Norfolk, datatype = "abundance")
#' ggiNEXT_ND(x = out1, type = 1)
#' ggiNEXT_ND(x = out1, type = 2)
#' ggiNEXT_ND(x = out1, type = 3)
#' }

#'
#' @export
ggiNEXT_ND <- function(outcome,type = 1,se = TRUE,facet.var = "None",color.var = "Assemblage",grey = FALSE, text.size = 18){

  iNEXT::ggiNEXT(outcome,type = type,facet.var = facet.var, color.var = color.var, se = se, grey = grey) + theme_bw() +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          legend.box.spacing = unit(0.4, "cm"),
          text=element_text(size= text.size),
          legend.key.width = unit(1,"cm")) + ylab("Network diversity")
}
# ggiNEXT -------------------------------------------------------------------
#' ggplot2 extension for an iNEXT object
#'
#' \code{ggiNEXT_PND}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXTlink}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
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
#' out1 <- iNEXT_ND(Norfolk, datatype = "abundance")
#' ggiNEXT_ND(x = out1, type = 1)
#' ggiNEXT_ND(x = out1, type = 2)
#' ggiNEXT_ND(x = out1, type = 3)
#' }

#' @export
ggiNEXT_PND <- function(outcome,type = 1, stript.size = 14, text.size = 14){
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
    ggplot(iNE, aes(x = m,y = PD)) + geom_line(aes(color = Region,linetype = method),size = 1.2) + facet_wrap(~Order.q) +
      geom_ribbon(aes(x = m,ymax = PD.UCL ,ymin = PD.LCL,fill = Region),alpha = 0.25) +
      geom_point(aes(x = m,y = PD ,color = Region,shape = Region),size = 5,data = iNE.sub) +
      theme(legend.position = "bottom",
            legend.title=element_blank(), strip.text = element_text(size = stript.size),
            text=element_text(size=text.size),
            legend.key.width = unit(0.8,"cm"))  +
      labs(x = "Number of individuals", y = "Network phylogenetic diversity", lty = "Method")+ theme_bw()
  }else if(type == 3){
    # coverage-based
    ggplot(iNE) + geom_line(aes(x = SC,y = PD,color = Region,linetype = method),size = 1.2) + facet_wrap(~Order.q) +
      geom_ribbon(aes(x = SC,ymax = PD.UCL ,ymin = PD.LCL,fill = Region),alpha = 0.25) +
      geom_point(aes(x = SC,y = PD ,color = Region,shape = Region),size = 5,data = iNE.sub) + theme_bw() +
      theme(legend.position = "bottom",
            legend.title=element_blank(), strip.text = element_text(size = stript.size),
            text=element_text(size=text.size),
            legend.key.width = unit(0.8,"cm"))  +
      labs(x = "Sample coverage", y = "Phylogenetic network diversity", lty = "Method")
  }
  # if(grey){
  #   g <- g +
  #     scale_fill_grey(start = 0, end = .4) +
  #     scale_colour_grey(start = .2, end = .2)
  # }

}
# ggiNEXT -------------------------------------------------------------------
#' ggplot2 extension for an iNEXT object
#'
#' \code{ggiNEXT_beta_link}: the \code{\link[ggplot2]{ggplot}} extension for
#' \code{\link{iNEXTlink}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
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
#' out = iNEXT_beta_link(Norfolk, coverage_expected, data_type=c('abundance', 'incidence_raw'), q = c(0, 1, 2),
#' level=c('taxonomic', 'phylogenetic', 'functional'), nboot = 20, conf = 0.95, max_alpha_coverage=F,
#'  by=c('coverage', 'size'),phy_tree=NULL, reftime = NULL)
#'
#' ggiNEXT_beta_link(out, type = 1,se = TRUE,facet = "None",color = "Assemblage",grey = FALSE, text_size = 10)
#' }


#' @export
ggiNEXT_beta_link <- function(output, type = c('B', 'D'), measurement = c('T', 'P', 'F_tau', 'F_AUC'),
                              scale='free', main=NULL, transp=0.4, stript.size = 11, text.size = 13){
  # if(length(outcome) == 1){ outcome = outcome}
  if (type == 'B'){

      gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
      alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
      beta =  lapply(output, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()
      beta = beta %>% filter(Method != 'Observed')
      beta[beta == 'Observed_alpha'] = 'Observed'

      df = rbind(gamma, alpha, beta)
      for (i in unique(gamma$Order)) df$Order[df$Order==i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))

      id_obs = which(df$Method == 'Observed')

      for (i in 1:length(id_obs)) {

        new = df[id_obs[i],]
        new$Coverage_expected = new$Coverage_expected - 0.0001
        new$Method = 'Interpolated'

        newe = df[id_obs[i],]
        newe$Coverage_expected = newe$Coverage_expected + 0.0001
        newe$Method = 'Extrapolated'

        df = rbind(df, new, newe)

      }

      if (measurement=='T') { ylab = "Taxonomic diversity" }
      if (measurement=='P') { ylab = "Phylogenetic Hill number" }
      if (measurement=='F_tau') { ylab = "Functional diversity (given tau)" }
      if (measurement=='F_AUC') { ylab = "Functional diversity (AUC)" }

    }

    if (type=='D'){

      C = lapply(output, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
      U = lapply(output, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
      V = lapply(output, function(y) y[["V"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
      S = lapply(output, function(y) y[["S"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()
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
        new$Coverage_expected = new$Coverage_expected - 0.0001
        new$Method = 'Interpolated'

        newe = df[id_obs[i],]
        newe$Coverage_expected = newe$Coverage_expected + 0.0001
        newe$Method = 'Extrapolated'

        df = rbind(df, new, newe)

      }

      if (measurement=='T') { ylab = "Taxonomic dissimilarity" }
      if (measurement=='P') { ylab = "Phylogenetic dissimilarity" }
      if (measurement=='F_tau') { ylab = "Functional dissimilarity (given tau)" }
      if (measurement=='F_AUC') { ylab = "Functional dissimilarity (AUC)" }

    }

    lty = c(Interpolated = "solid", Extrapolated = "dotted")
    df$Method = factor(df$Method, levels = c('Interpolated', 'Extrapolated', 'Observed'))

    double_size = unique(df[df$Method=="Observed",]$Size)*2
    double_extrepolation = df %>% filter(Method=="Extrapolated" & round(Size) %in% double_size)

    # ggplot(data = df, aes(x = Coverage_expected, y = Estimate)) +
    plot = ggplot(data = df, aes(x = Coverage_expected, y = Estimate, col = Region)) +
      # geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=transp) +
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region), alpha=transp) +
      geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) + scale_linetype_manual(values = lty) +
      # geom_line(lty=2) +
      geom_point(data = subset(df, Method=='Observed' & div_type=="Gamma"),shape=19, size=3) +
      geom_point(data = subset(df, Method=='Observed' & div_type!="Gamma"),shape=1, size=3,stroke=1.5)+
      geom_point(data = subset(double_extrepolation, div_type == "Gamma"),shape=17, size=3) +
      geom_point(data = subset(double_extrepolation, div_type!="Gamma"),shape=2, size=3,stroke=1.5) +
      facet_grid(div_type~Order, scales = scale) +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank(),
            strip.text = element_text(size = stript.size), text = element_text(size = text.size),
            legend.key.width = unit(1,"cm")) +
      labs(x='Sample coverage', y=ylab, title=main)
    if(length(output) == 1){
      plot = plot + guides(col=FALSE, fill = FALSE)
    }
    return(plot)
}
# ggAsyD -------------------------------------------------------------------
#' ggplot for Asymptotic Network diversity
#'
#' \code{ggAsyND} Plots q-profile based on the outcome of \code{AsyND} using the ggplot2 package.\cr
#' It will only show the confidence interval of 'Estimated'.
#'
#' @param outcome the outcome of the functions \code{AsyND} .\cr
#' @return a figure of estimated sample completeness with order q\cr\cr
#'
#' @examples
#' \dontrun{
#' ## example for abundance-based data
#' ## Ex.1
#' data(Norfolk)
#' out1 <- AsyND(Norfolk, datatype = "abundance")
#' ggAsyND(out1)
#' }

#' @export
ggAsyND <- function(outcome, text.size = 14){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  if (sum(unique(outcome$method) %in% c("Estimated", "Empirical")) == 0)
    stop("Please use the outcome from specified function 'AsyD'")

  ggplot(outcome, aes(x = Order.q, y = qD, colour = Assemblage, lty = method)) +
    geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = outcome[outcome$method == "Estimated", ],
                aes(ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage), alpha = 0.2, linetype = 0) +
    scale_fill_manual(values = cbPalette) +
    scale_linetype_manual(values = c(Estimated = 1, Empirical = 2)) +
    labs(x = "Order q", y = "Network diversity") + theme(text = element_text(size = 10)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.1, "cm"), legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-10,-10, -5, -10),
          text = element_text(size = text.size)
          )
}

# ggAsyD -------------------------------------------------------------------
#' ggplot for Asymptotic diversity
#'
#' \code{ggAsyPND} Plots q-profile based on the outcome of \code{AsyD} using the ggplot2 package.\cr
#' It will only show the confidence interval of 'Estimated'.
#'
#' @param outcome the outcome of the functions \code{AsyD} .\cr
#' @return a figure of estimated sample completeness with order q\cr\cr
#'
#' @examples
#' \dontrun{
#' data(puerto.rico)
#' out1 <- AsyPND(puerto.rico$data, datatype = "abundance")
#' ggAsyPND(out1)
#' }
#' @export
ggAsyPND <- function(outcome, text_size = 14){
  table = outcome
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  if (sum(unique(table$method) %in% c("Estimated", "Empirical")) == 0)
    stop("Please use the outcome from specified function 'AsyD'")

  ggplot(table) +
    geom_line(aes(x = Order.q, y = Estimate,lty = method, color = Region),lwd = 1.4) +
    geom_ribbon(aes(x = Order.q, ymin = LCL, ymax = UCL,fill = Region,lty = method), alpha = 0.25) + theme_bw() +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          text=element_text(size=text_size),
          legend.key.width = unit(1,"cm"))  +
    labs(x = "Order q", y = "Network phylogenetic diversity", lty = "Method") + scale_linetype_manual(values=c("dashed","solid"))
}

# AsyND -------------------------------------------------------------------
#' Asymptotic diversity q profile
#'
#' \code{AsyND} The estimated and empirical diversity of order q
#'
#' @param outcome the outcome of the functions \code{AsyD} .\cr
#' @return a table of Asymptoti network diversity q profile
#'
#' @examples
#' \dontrun{
#' data(Norfolk)
#' out1 <- AsyND(Norfolk,datatype = "abundance")
#' ggAsyPND(out1)
#' }
#'
#' @export
AsyND <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95){
  lapply(1:length(data), function(i) {
    x = data[[i]]
    assemblage = names(data)[[i]]
    # tmp <- c(as.matrix(x))
    # tmp = x
    res = MakeTable_Proposeprofile(data = x, B = nboot, q, conf = conf)%>%
      rename("qD"="Estimate", "qD.LCL"="LCL", "qD.UCL"="UCL")%>%
      mutate(Assemblage = assemblage, method = "Estimated")%>%filter(Target == "Diversity")%>%select(-Target)
    return(res)
    })%>%do.call("rbind",.)
}
# AsyPND -------------------------------------------------------------------
#' Asymptotic Phylogenetic diversity q profile
#'
#' \code{AsyPND} The estimated and empirical diversity of order q
#'
#' @param outcome the outcome of the functions \code{AsyD} .\cr
#' @return a table of Asymptotic network diversity q profile
#'
#' @examples
#' \dontrun{
#' ## Type (1) example for abundance-based data
#' ## Ex.1
#' data(puerto.rico$data)
#' out1 <- AsyND(puerto.rico$data, datatype = "abundance", row.tree = row.tree,col.tree = col.tree)
#' ggAsyPND(out1)
#' }
#' @export
AsyPND <- function(data = puerto.rico$data, q = seq(0, 2, 0.2), datatype = "abundance",
                   row.tree = NULL, col.tree = NULL,
                   nboot = 50, conf = 0.95){

  NetDiv <- get.netphydiv(data = data,q = q,B = nboot,row.tree = row.tree,col.tree = col.tree,conf = conf)%>%
    filter(method == "Estimate")

  return(NetDiv)
}



# ObsND -------------------------------------------------------------------
#' Empirical diversity q profile
#'
#' \code{ObsND} The estimated and empirical diversity of order q
#'
#' @param outcome the outcome of the functions \code{ObsND} .\cr
#' @return a table of Asymptotic network diversity q profile
#'
#' @examples
#' \dontrun{
#' ## Example for abundance-based data
#' data(Norfolk)
#' out1 <- ObsND(Norfolk, datatype = "abundance", nboot = 30)
#' ggObsND(out1)
#' }
#' @export
ObsND <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95){

  lapply(1:length(data), function(i){
    x = data[[i]]
    assemblage = names(data)[[i]]
    tmp <- c(as.matrix(x))
    ## nboot has to larger than 0
    res = MakeTable_Empericalprofile(data = x, B = nboot, q, conf = conf)%>%
      rename("qD"="Emperical", "qD.LCL"="LCL", "qD.UCL"="UCL")%>%
      mutate(Assemblage = assemblage, method = "Empirical")%>%filter(Target == "Diversity")%>%select(-Target)
    return(res)
  })%>%do.call("rbind",.)
}
# ObsPND -------------------------------------------------------------------
#' Empirical diversity q profile
#'
#' \code{ObsND} The estimated and empirical diversity of order q
#'
#' @param outcome the outcome of the functions \code{ObsND} .\cr
#' @return a table of Asymptotic network diversity q profile
#'
#' @examples
#' \dontrun{
#' ## Example for abundance-based data
#' data(puerto.rico$data)
#' out1 <- ObsPND(puerto.rico$data, datatype = "abundance", row.tree = row.tree,col.tree = col.tree)
#' ggObsPND(out1)
#' }
#' @export
ObsPND <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, row.tree = NULL, col.tree = NULL){
  NetDiv <- get.netphydiv(data = data,q = q,B = nboot,row.tree = row.tree,col.tree = col.tree,conf = conf)%>%
    filter(method == "Empirical")

  return(NetDiv)
}

# estimateD  -------------------------------------------------------------------
#' Compute species diversity with a particular of sample size/coverage
#'
#' \code{estimateND} computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#'
#' @param outcome the outcome of the functions \code{ObsND} .\cr
#' @return a data.frame of species diversity table including the sample size, sample coverage, method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#'
#' @examples
#' \dontrun{
#' data(Norfolk)
#' out1 <- estimateND(Norfolk, datatype="abundance", base="coverage",level=0.7, nboot = 30,conf=0.95)
#' out2 <- estimateND(Norfolk, datatype="abundance", base="size",level=0.7, nboot = 30,conf=0.95)
#' }

#'
#' @export
estimateND = function(dat,q = c(0, 1, 2),datatype = "abundance",base = "size",
                      level = NULL,nboot = 50,conf = 0.95){
  lapply(1:length(dat), function(i){
    x = dat[[i]]
    assemblage = names(dat)[[i]]
    long = as.matrix(x)%>%c()
    estimateD(long, q=q,datatype=datatype, base=base,level=level, nboot = nboot,conf=conf)%>%
      mutate(Assemblage = assemblage)
  })%>%do.call("rbind",.)
}
# estimatedPND  -------------------------------------------------------------------
#' Compute species diversity with a particular of sample size/coverage
#'
#' \code{estimatedPND} computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#'
#' @param outcome the outcome of the functions \code{ObsND} .\cr
#' @return a data.frame of species diversity table including the sample size, sample coverage, method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#'
#' @examples
#' \dontrun{
#' data(puerto.rico$data)
#' out <- estimatedPND(puerto.rico$data, q = c(0,1,2), datatype = "abundance", row.tree = rowtree, col.tree = coltree)
#' out
#' }

#'
#' @export
estimatePND <- function(data,q = c(0, 1, 2),datatype = "abundance",
                         row.tree = NULL, col.tree =NULL,
                         level = NULL,nboot = 30,conf = 0.95){
  q <- unique(ceiling(q))
  ci <- qnorm(conf/2+0.5)

  if (is.null(level)) {
    if (datatype == "abundance") {
      level <- sapply(data, function(x) {
        ni <- sum(x)
        iNEXTPD2:::Coverage(data = x, datatype = datatype, m = 2 * ni, nt = ni)
      })
    }
    else if (datatype == "incidence_raw") {
      level <- sapply(data, function(x) {
        ni <- ncol(x)
        Coverage(data = x, datatype = datatype, m = 2 *
                   ni, nt = ni)
      })
    }
    level <- min(level)
  }
  # level =0.7
  res = lapply(1:length(data), function(i){
    x = data[[i]]
    assemblage = names(data)[[i]]
    long = as.matrix(x)%>%c()
    m_target = Coverage_to_size(x, C = level)
    # if(base == "size"){m_target = level}
    inex <- function(data,m,q,B,row.tree = NULL,col.tree = NULL) {
      data <- as.matrix(data)
      n <- sum(data)
      phydata <- create.aili(data,row.tree = row.tree,col.tree = col.tree)
      tbar <- sum(phydata$branch.length*phydata$branch.abun)/n
      boot.sam <- sample.boot.phy(data,B,row.tree = row.tree,col.tree = col.tree)
      sc <- PhD:::Coverage(data,datatype = "abundance",m,nt =n)
      # sc <- coverage(data,m)
      sc.sd <- lapply(boot.sam,function(x){
        x <- x[x$tgroup == "Tip",]$branch.abun
        PhD:::Coverage(x,datatype = "abundance",m,nt =n)
      })
      sc.sd <- do.call(cbind,sc.sd)
      sc.sd <- sapply(1:length(m), function(x){
        sd(sc.sd[x,])
      })
      sc.table <- data.frame(m=m,SC = sc, SC.UCL = sc+ci * sc.sd,SC.LCL = sc - ci * sc.sd)
      out <- lapply(q, function(x){
        PD <- lapply(m,function(y){
          my_PhD.m.est(ai = phydata$branch.abun, Lis = phydata$branch.length, m = y, q = x, nt = n, cal = 'PD')/tbar
          # PhD:::PhD.m.est(phydata,y,x,datatype = "abundance",nt = n)/tbar
        })%>%unlist()
        PD.sd <- lapply(boot.sam, function(z){
          tmp <- lapply(m,function(y){
            # PhD:::PhD.m.est(z,y,x,datatype = "abundance",nt = n)/tbar
            my_PhD.m.est(ai = z$branch.abun, Lis = z$branch.length, m = y, q = x, nt = n, cal = 'PD')/tbar
          })
          unlist(tmp)
        })
        PD.sd <- do.call(cbind,PD.sd)
        PD.sd <- sapply(1:length(m), function(x){
          sd(PD.sd[x,])
        })
        PD.table <- data.frame(m=m,method = ifelse(m<n,"interpolated",ifelse(n == m,"observed","extrapolated")),
                               Order.q = x,PD = PD, PD.UCL = PD+ci * PD.sd,PD.LCL = PD - ci * PD.sd)
        out <- left_join(PD.table,sc.table)
        out
      })
      do.call(rbind,out)
    }

    inex(data = x,m = m_target, q,B = nboot,row.tree,col.tree)%>%
      mutate(Assemblage = assemblage)

  })%>%do.call("rbind",.)
  return(res)
}


# Evenness  -------------------------------------------------------------------
#' Evenness Estimation of Evenness with order q
#'
#' \code{Evenness} computes Evenness Estimation of Evenness with order q.
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
#' Est <- NetEvenness(x = Norfolk, datatype = "abundance", q = c(0,1,2), nboot = 30, method = "Estimated")
#' Emp <- NetEvenness(x = Norfolk, datatype = "abundance", q = c(0,1,2), nboot = 30, method = "Empirical")
#' Est
#' Emp
#' ggNetEven(Est)
#' ggNetEven(Est)
#' @export
# x = Norfolk
# tmp1 = NetEvenness(Norfolk, q=  seq(0,2,0.2), E = 1:5, method = "Estimated", C = 0.9)
# tmp1 = NetEvenness(Norfolk, q=  seq(0,2,0.2), E = 1:5, method = "Empirical")
NetEvenness <- function(x,q = seq(0, 2, 0.2),
                     datatype = "abundance",
                     method = "Estimated",
                     nboot = 30,
                     conf = 0.95,
                     E.class = c(1:5),
                     C = NULL){

  long = lapply(x, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})
  # EVEN <- lapply(seq_along(long),
  #                           function(i){
  #                             res = iNEXT4steps::Evenness(long[[i]], q = q,datatype = datatype,
  #                                                         method = method, nboot=nboot, E.class = E, C = C)
  #                             res[["Coverage"]] = NULL
  #                             tab = lapply(E.class, function(j){
  #                               tmp = res[[j]]%>%mutate(class = paste0("E", j))
  #                               return(tmp)
  #                             })%>%do.call("rbind",.)%>%
  #                               mutate(Assemblage = names(long)[[i]])
  #                             return(tab)
  #                           })%>%do.call("rbind",.)
  EVEN <- lapply(E.class, function(e){
    each_class = lapply(seq_along(long), function(i){
      res = iNEXT4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                  method = method, nboot=nboot, E.class = e, C = C)
      if(method == "Empirical") index = 1
      if(method == "Estimated") index = 2
      return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
    })%>%do.call("rbind",.)

    each_class%>%mutate(class = paste0("E",e))
  })
  names(EVEN) = paste0("E",E.class)
  return(EVEN)
}
# ggNetEven -------------------------------------------------------------------
#' ggplot for Evenness
#
#' \code{ggNetEven} The figure for estimation of Evenness with order q\cr
#'
#' @param output a table generated from Evenness function\cr
#' @return a figure of estimated sample completeness with order q\cr
#'
#' @examples
#' \dontrun{
#' data(Norfolk)
#' Est <- NetEvenness(x = Norfolk, datatype = "abundance", q = c(0,1,2), nboot = 30, method = "Estimated")
#' Emp <- NetEvenness(x = Norfolk, datatype = "abundance", q = c(0,1,2), nboot = 30, method = "Empirical")
#' Est
#' Emp
#' ggNetEven(Est)
#' ggNetEven(Est)
#' @references
#' Chao,A.and Ricotta,C.(2019).Quantifying evenness and linking it to diversity, beta diversity, and similarity.
#' @export
ggNetEven <- function(output){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  classdata = cbind(do.call(rbind, output))%>%
    rename("Method"="method")

  fig = classdata%>%
    ggplot(aes(x=Order.q, y=Evenness, colour=Assemblage, lty = Method)) +
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = classdata %>% filter(Method=="Estimated"),
                aes(ymin=Even.LCL, ymax=Even.UCL, fill=Assemblage),
                alpha=0.2, linetype=0) +
    geom_ribbon(data = classdata %>% filter(Method=="Empirical"),
                aes(ymin=Even.LCL, ymax=Even.UCL, fill=Assemblage),
                alpha=0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    labs(x="Order q", y="Evenness") +
    # theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2,"cm"),
          # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
          legend.title = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-5,-10),
          text = element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt")
    )

  if (length(output) != 1) fig = fig +
    facet_wrap(~class) +
    theme(strip.text.x = element_text(size=12, colour = "purple", face="bold"))

  return(fig)
}

#' Evenness main function
#'
#' \code{phyNetEvenness} Estimation (Empirical) of Evenness with order q
#'
#' R scipts "Evenness" for Chao and Ricotta (2019) Ecology paper.
#' This R code is for computing Figures 2, 3 and 4 of Chao and Ricotta (2019) paper.
#' installed and loaded before running the scripts.

#' @param x a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr
#' @param q a integer vector of the order of Hill number\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).\cr
#' @param method a binary calculation method with 'Estimated' or 'Empirical'\cr
#' @param nboot an integer specifying the number of bootstrap replications, default is 30.\cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.\cr
#' @param E.class a integer vector between 1 to 6
#' @param C a standardized coverage for calculating evenness index
#' @return A list of estimated(empirical) evenness with order q.\cr
#'         Different lists represents different classes of Evenness.\cr
#'         Each list is combined with order.q and sites.\cr
#'         If "method" is estimated, then fist list will be named "C" which means the
#'         maximum standardized coverage between all double reference sample size.\cr\cr
#' \code{$summary} individual summary of 4 steps of data. \cr\cr
#'
#' @examples
#' \dontrun{
#' data(
#' Est <- phyNetEvenness(x = puerto.rico$data, q = c(0,1,2), nboot = 30, method = "Estimated", col.tree = puerto.rico$col.tree, row.tree = puerto.rico$row.tree)
#' Emp <- phyNetEvenness(x = puerto.rico$data, q = c(0,1,2), nboot = 30, method = "Empirical", col.tree = puerto.rico$col.tree, row.tree = puerto.rico$row.tree)
#' Est
#' Emp
#' ggNetEven(Est)
#' ggNetEven(Est)
#' @references
#' Chao,A.and Ricotta,C.(2019).Quantifying evenness and linking it to diversity, beta diversity, and similarity.
#' @export
phyNetEvenness <- function(x,row.tree = NULL,col.tree = NULL,
                           q = seq(0, 2, 0.2),
                           datatype = "abundance",
                           method = "Estimated",
                           nboot = 30,
                           conf = 0.95,
                           E.class = c(1:5),
                           C = NULL){
  long = lapply(x, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})



  EVEN <- lapply(seq_along(long),
                            function(i){
                              res = iNEXT4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                                          method = method, nboot=nboot, E.class = 1:5, C = C)
                              res[["Coverage"]] = NULL
                              tab = lapply(E.class, function(j){
                                tmp = res[[j]]%>%mutate(class = paste0("E", j))
                                return(tmp)
                              })%>%do.call("rbind",.)%>%
                                mutate(Assemblage = names(long)[[i]])
                              return(tab)
                            })%>%do.call("rbind",.)
  return(EVEN)
}


# ggEven <- function(output = tmp1) {
#   if (names(output[1]) == "Coverage")  output = output[-1]
#   cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
#                      "#330066", "#CC79A7", "#0072B2", "#D55E00"))
#   classdata = cbind(do.call(rbind, output),
#                     class = rep(names(output), each=nrow(output[[1]])))%>%
#     rename("Method"="method")
#
#   fig = classdata%>%
#     ggplot(aes(x=Order.q, y=Evenness, colour=Assemblage, lty = Method)) +
#     geom_line(size=1.2) +
#     scale_colour_manual(values = cbPalette) +
#     geom_ribbon(data = classdata %>% filter(Method=="Estimated"),
#                 aes(ymin=Even.LCL, ymax=Even.UCL, fill=Assemblage),
#                 alpha=0.2, linetype=0) +
#     geom_ribbon(data = classdata %>% filter(Method=="Empirical"),
#                 aes(ymin=Even.LCL, ymax=Even.UCL, fill=Assemblage),
#                 alpha=0.2, linetype=0) +
#     scale_fill_manual(values = cbPalette) +
#     labs(x="Order q", y="Evenness") +
#     # theme_bw(base_size = 18) +
#     theme(text=element_text(size=18)) +
#     theme(legend.position = "bottom", legend.box = "vertical",
#           legend.key.width = unit(1.2,"cm"),
#           # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
#           legend.title = element_blank(),
#           legend.margin = margin(0,0,0,0),
#           legend.box.margin = margin(-10,-10,-5,-10),
#           text = element_text(size=12),
#           plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt")
#     )
#
#   if (length(output) != 1) fig = fig +
#     facet_wrap(~class) +
#     theme(strip.text.x = element_text(size=12, colour = "purple", face="bold"))
#
#   return(fig)
# }

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
# MakeTable_Proposeprofile = function(data, B, q, conf){
#   Diversity = bootstrap_forq(data, B, q, conf, Diversity_profile)
#   Entropy = bootstrap_forq(data, B, q, conf, Diversity_Tsallis)
#   # tmp <- Diversity_Tsallis(Diversity[,1],q)
#   # Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
#   output = rbind(data.frame("Order.q" = q,"Target"="Diversity","Estimate"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4]),
#                  data.frame("Order.q" = q,"Target"="Entropy","Estimate"=Entropy[,1],"s.e."=Entropy[,2],"LCL"=Entropy[,3],"UCL"=Entropy[,4]))
#   output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)
#
#   return(output)
# }
# MakeTable_Empericalprofile = function(data, B, q, conf){
#   Diversity = bootstrap_forq( data, B,q,conf,Diversity_profile_MLE)
#   Entropy = bootstrap_forq( data, B,q,conf,Diversity_Tsallis_MLE)
#   # tmp <- Diversity_Tsallis(Diversity[,1],q)
#   # Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
#   output = rbind(data.frame("Order.q" = q,"Target"="Diversity","Emperical"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4]),
#                  data.frame("Order.q" = q,"Target"="Entropy","Emperical"=Entropy[,1],"s.e."=Entropy[,2],"LCL"=Entropy[,3],"UCL"=Entropy[,4]))
#   output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)
#
#   return(output)
# }
# bootstrap_forq = function(data,B,q,conf,FUNNAME){
#   data <- data[data!=0]
#   n <- sum(data)
#   f1 = sum(data==1); f2 = sum(data==2)
#   f0 = ceiling(ifelse( f2>0, (n-1)*f1^2/n/2/f2, (n-1)*f1*(f1-1)/2/n ))
#   C_hat = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
#   lamda_hat = (1-C_hat)/sum((data/n)*(1-data/n)^n)
#   pi_hat = (data/n)*(1-lamda_hat*((1-data/n)^n))
#   p_hat = c( pi_hat, rep( (1-C_hat)/f0, f0 ))
#   random = rmultinom( B, n, p_hat )
#   #Bt_estimate <- sapply(c(1:B),function(i) FUNNAME(random[,i],q))
#   Bt_estimate <- apply(random,MARGIN = 2,function(i) FUNNAME(i,q))
#   estimate <- FUNNAME(data,q)
#   #Interval_mean = apply( Bt_estimate, 1, mean)
#   Interval_mean = rowMeans(Bt_estimate)
#   Interval_sd = apply(Bt_estimate, 1, sd)
#   Interval_quantileL = apply( Bt_estimate, 1, quantile, p=(1-conf)/2)
#   Interval_quantileU = apply( Bt_estimate, 1, quantile, p=1-(1-conf)/2)
#   Upper_bound = estimate+Interval_quantileU-Interval_mean
#   Lower_bound = estimate+Interval_quantileL-Interval_mean
#   result <- cbind("estimate"=estimate,"sd"=Interval_sd,"LCL"=Lower_bound,"UCL"=Upper_bound)
#   result
# }
