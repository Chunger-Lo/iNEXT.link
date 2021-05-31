#' #' iNterpolation and EXTrapolation Beta diversity
#' #'
#' #' \code{iNEXTbeta}: Interpolation and extrapolation of Beta diversity with order q
#' #'
#' #' @param data a \code{matrix}, \code{data.frame} (species by assemblages), or \code{list} of species abundance/incidence raw data.\cr
#' #' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' #' @param q a numerical vector specifying the diversity orders.
#' #' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' #' @param level a sequence specifying the particular sample coverages (between 0 and 1). Calculating asymptotic diversity if \code{C}=1. \cr
#' #' @param nboot an integer specifying the number of replications for bootstrapping.
#' #' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' #' @param PDtree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. It is necessary when \code{diversity = 'PD'}.
#' #' @param PDreftime Select a reference time point for \code{diversity = 'PD'}. Default is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in the pooled assemblage. Only required when \code{diversity = 'PD'}.
#' #' @param PDtype Select phylogenetic diversity type: \code{PDtype = "PD"} for Chao et al. (2010) phylogenetic diversity and \code{PDtype = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Only required when \code{diversity = 'PD'}. Default is \code{"meanPD"}.
#' #' @param FDdistM a pair wise distance matrix for all pairs of observed species in the pooled assemblage. It is necessary when diversity = 'FD'.
#' #' @param FDtype a binary selection for FD. \code{FDtype = "single"} computes diversity under a certain threshold value. \code{FDtype = "AUC"} computes an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' #' @param FDtau a value between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}. Only required when \code{diversity = 'FD'} and \code{FDtype = "single"}.
#' #' @param FDcut_number a numeric number to split zero to one into several equal-spaced length. Only required when \code{diversity = "FD"} and \code{FDtype = "AUC"}. Default is 30.
#' #'
#' #' @import tidyverse
#' #' @import magrittr
#' #' @import ggplot2
#' #' @import abind
#' #' @import ape
#' #' @import phytools
#' #' @import phyclust
#' #' @import tidytree
#' #' @import colorRamps
#' #' @import iNEXT.3D
#' #' @import future.apply
#' #' @import ade4
#' #' @import tidyr
#' #'
#' #' @return A list of seven lists with three-diversity and four-dissimilarity.
#' #'
#' #' @examples
#' #' ## Taxonomic diversity for abundance data
#' #' data(ant)
#' #' output1 = iNEXTbeta(data = ant, diversity = 'TD', datatype = 'abundance', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 20, conf = 0.95)
#' #' output1
#' #'
#' #'
#' #' ## Taxonomic diversity for incidence data
#' #' data(fish)
#' #' output2 = iNEXTbeta(data = fish, diversity = 'TD', datatype = 'incidence_raw', level = seq(0.7, 0.95, 0.05),
#' #'                     nboot = 20, conf = 0.95)
#' #' output2
#' #'
#' #'
#' #' ## Phylogenetic diversity for abundance data
#' #' data(ant)
#' #' data(ant_tree)
#' #' output3 = iNEXTbeta(data = ant, diversity = 'PD', datatype = 'abundance',  level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 20, conf = 0.95, PDtree = ant_tree, PDreftime = NULL, PDtype = 'PD')
#' #' output3
#' #'
#' #'
#' #' ## Phylogenetic diversity for incidence data
#' #' data(fish)
#' #' data(fish_tree)
#' #' output4 = iNEXTbeta(data = fish, diversity = 'PD', datatype = 'incidence_raw', level = seq(0.7, 0.95, 0.05),
#' #'                     nboot = 0, conf = 0.95, PDtree = fish_tree, PDreftime = NULL, PDtype = 'PD')
#' #' output4
#' #'
#' #'
#' #' ## Functional diversity for abundance data under single threshold
#' #' data(ant)
#' #' data(ant_distM)
#' #' output5 = iNEXTbeta(data = ant, diversity = 'FD', datatype = 'abundance', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 20, conf = 0.95, FDdistM = ant_distM, FDtype = 'single', FDtau = NULL)
#' #' output5
#' #'
#' #'
#' #' ## Functional diversity for incidence data under single threshold
#' #' data(fish)
#' #' data(fish_distM)
#' #' output6 = iNEXTbeta(data = fish, diversity = 'FD', datatype = 'incidence_raw', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 20, conf = 0.95, FDdistM = fish_distM, FDtype = 'single', FDtau = NULL)
#' #' output6
#' #'
#' #'
#' #' ## Functional diversity for abundance data with thresholds integrating from 0 to 1
#' #' data(ant)
#' #' data(ant_distM)
#' #' output7 = iNEXTbeta(data = ant, diversity = 'FD', datatype = 'abundance', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 10, conf = 0.95, FDdistM = ant_distM, FDtype = 'AUC', FDcut_number = 30)
#' #' output7
#' #'
#' #'
#' #' ## Functional diversity for incidence data with thresholds integrating from 0 to 1
#' #' data(fish)
#' #' data(fish_distM)
#' #' output8 = iNEXTbeta(data = fish, diversity = 'FD', datatype = 'incidence_raw', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 10, conf = 0.95, FDdistM = fish_distM, FDtype = 'AUC', FDcut_number = 30)
#' #' output8
#' #'
#' #' @export
#' iNEXTbeta = function(data, diversity = c('TD', 'PD', 'FD'), q = c(0, 1, 2), datatype = c('abundance', 'incidence_raw'),
#'                      level, nboot = 20, conf = 0.95, PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD',
#'                      FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, FDcut_number = 30) {
#'   max_alpha_coverage = F
#'   if (datatype == 'abundance') {
#'
#'     if( class(data) == "data.frame" | class(data) == "matrix" ) data = list(Region_1 = data)
#'
#'     if(class(data) == "list"){
#'       if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
#'       Ns = sapply(data, ncol)
#'       data_list = data
#'     }
#'
#'   }
#'
#'   if (datatype == 'incidence_raw') {
#'
#'     if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
#'     Ns = sapply(data, length)
#'     data_list = data
#'
#'   }
#'
#'
#'   if (is.null(conf)) conf = 0.95
#'   tmp = qnorm(1 - (1 - conf)/2)
#'
#'   if (diversity == 'FD' & FDtype == 'single' & is.null(FDtau) == T) {
#'     if (datatype == 'abundance') {
#'       pdata <- lapply(data_list, rowSums) %>% do.call(cbind,.) %>% apply(., 2, function(x) x/sum(x)) %>% rowMeans()
#'     } else if (datatype == 'incidence_raw') {
#'       pdata <- sapply(data_list, function(x) {tmp = Reduce('+', x); tmp[tmp > 1] = 1; rowSums(tmp) }) %>% apply(., 2, function(x) x/sum(x)) %>% rowMeans()
#'     }
#'     FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
#'   }
#'
#'   if (diversity == 'PD') {
#'
#'     if (datatype == 'abundance') pool.data = do.call(cbind, data_list) %>% rowSums
#'
#'     if (datatype == 'incidence_raw') pool.data = do.call(cbind,lapply(data_list, function(x) do.call(cbind,x)) ) %>% rowSums
#'
#'     pool.name = names(pool.data[pool.data>0])
#'     tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
#'     mytree = drop.tip(PDtree, tip)
#'     H_max = get.rooted.tree.height(mytree)
#'
#'     if(is.null(PDreftime)) { reft = H_max
#'     } else if (PDreftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
#'     } else { reft = PDreftime }
#'
#'   }
#'
#'   for_each_region = function(data, region_name, N) {
#'
#'     #data
#'     if (datatype == 'abundance') {
#'
#'       n = sum(data)
#'       data_gamma = rowSums(data)
#'       data_gamma = data_gamma[data_gamma>0]
#'       data_alpha = as.matrix(data) %>% as.vector
#'
#'       ref_gamma = iNEXT.3D:::Chat.Ind(data_gamma, n)
#'       ref_alpha = iNEXT.3D:::Chat.Ind(data_alpha, n)
#'       ref_alpha_max = iNEXT.3D:::Chat.Ind(data_alpha, n*2)
#'       ref_gamma_max = iNEXT.3D:::Chat.Ind(data_gamma, n*2)
#'
#'       level = c(level, ref_gamma, ref_alpha, ref_alpha_max,ref_gamma_max) %>% sort %>% unique
#'       # level = level[level<1]
#'
#'       m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
#'       m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
#'
#'     }
#'
#'     if (datatype == 'incidence_raw') {
#'
#'       sampling_units = sapply(data, ncol)
#'       if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
#'       if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
#'
#'       gamma = Reduce('+', data)
#'       gamma[gamma>1] = 1
#'       data_gamma_raw = gamma
#'       data_gamma_freq = c(n, rowSums(gamma))
#'
#'       data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
#'
#'       # data_gamma_freq = data_gamma_freq[data_gamma_freq>0]
#'       # data_alpha_freq = data_alpha_freq[data_alpha_freq>0]
#'
#'       data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
#'
#'       ref_gamma = iNEXT.3D:::Chat.Sam(data_gamma_freq, n)
#'       ref_alpha = iNEXT.3D:::Chat.Sam(data_alpha_freq, n)
#'       ref_alpha_max = iNEXT.3D:::Chat.Sam(data_alpha_freq, n*2)
#'       ref_gamma_max = iNEXT.3D:::Chat.Sam(data_gamma_freq, n*2)
#'
#'       level = c(level, ref_gamma, ref_alpha, ref_alpha_max) %>% sort %>% unique
#'       # level = level[level < 1]
#'
#'       m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq, i, datatype='incidence_freq'))
#'       m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq, i, datatype='incidence_raw'))
#'
#'     }
#'
#'
#'
#'     if (diversity == 'TD') {
#'
#'       if (datatype == 'abundance') {
#'
#'         gamma = lapply(1:length(level), function(i){
#'           estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
#'         }) %>% do.call(rbind,.)
#'
#'         alpha = lapply(1:length(level), function(i){
#'           estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
#'         }) %>% do.call(rbind,.)
#'
#'       }
#'
#'       if (datatype == 'incidence_raw') {
#'
#'         gamma = lapply(1:length(level), function(i){
#'           estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
#'         }) %>% do.call(rbind,.)
#'
#'         alpha = lapply(1:length(level), function(i){
#'           estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
#'         }) %>% do.call(rbind,.)
#'
#'
#'
#'       }
#'
#'       gamma = (cbind(level = rep(level, each=length(q)), gamma[,-c(1,2,8,9)]) %>%
#'                  mutate(Method = ifelse(level>=ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))
#'       )[,c(6,5,4,1,2,3)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
#'
#'       # for (i in 0:2) gamma$Order[gamma$Order==paste0('q = ', i)] = i
#'       # gamma$Order = as.numeric(gamma$Order)
#'
#'       if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order == 0) & (gamma$level > ref_alpha_max)) else under_max_alpha = gamma$level > 0
#'       gamma = gamma[under_max_alpha,]
#'
#'
#'
#'       alpha = (cbind(level = rep(level, each = length(q)), alpha[,-c(1,2,8,9)]) %>%
#'                  mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))
#'       )[,c(6,5,4,1,2,3)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
#'
#'       alpha$Estimate = alpha$Estimate / N
#'
#'       # for (i in 0:2) alpha$Order[alpha$Order == paste0('q = ', i)] = i
#'       # alpha$Order = as.numeric(alpha$Order)
#'
#'       alpha = alpha[under_max_alpha,]
#'
#'
#'
#'       beta = alpha
#'       beta$Estimate = gamma$Estimate/alpha$Estimate
#'       beta[beta == "Observed"] = "Observed_alpha"
#'       beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate),
#'                                      Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))
#'
#'       C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
#'       U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
#'       V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
#'       S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))
#'
#'       if(nboot>1){
#'
#'
#'         # start = Sys.time()
#'         se = future_lapply(1:nboot, function(i){
#'
#'           if (datatype == 'abundance') {
#'
#'             bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
#'             bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
#'
#'             bootstrap_data_gamma = rowSums(bootstrap_sample)
#'             bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
#'             bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
#'             bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
#'
#'             m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma, i, datatype='abundance'))
#'
#'             gamma = lapply(1:length(level), function(i){
#'               estimate3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = m_gamma[i], nboot = 0)
#'             }) %>% do.call(rbind,.)
#'
#'             alpha = lapply(1:length(level), function(i){
#'               estimate3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
#'             }) %>% do.call(rbind,.)
#'
#'             beta_obs = (obs3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", nboot = 0) %>% select(qD) /
#'                           (obs3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", nboot = 0) %>% select(qD) / N)) %>% unlist()
#'
#'           }
#'
#'           if (datatype == 'incidence_raw') {
#'
#'             bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
#'
#'             raw = lapply(1:ncol(bootstrap_population), function(j){
#'
#'               lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
#'
#'             })
#'
#'             gamma = Reduce('+', raw)
#'             gamma[gamma > 1] = 1
#'             bootstrap_data_gamma_freq = c(n, rowSums(gamma))
#'
#'             bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
#'
#'             bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
#'             bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
#'
#'             m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma_freq, i, datatype='incidence_freq'))
#'
#'             gamma = lapply(1:length(level), function(i){
#'               estimate3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = m_gamma[i], nboot = 0)
#'             }) %>% do.call(rbind,.)
#'
#'             alpha = lapply(1:length(level), function(i){
#'               estimate3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
#'             }) %>% do.call(rbind,.)
#'
#'             beta_obs = (obs3D(as.numeric(bootstrap_data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0) %>% select(qD) /
#'                           (obs3D(as.numeric(bootstrap_data_alpha_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0) %>% select(qD) / N)) %>% unlist()
#'
#'           }
#'
#'           gamma = gamma[,c(6,3,7)]$qD[under_max_alpha]
#'
#'           alpha = alpha[,c(6,3,7)]$qD[under_max_alpha]
#'           alpha = alpha / N
#'
#'           beta = gamma/alpha
#'
#'           gamma = c(gamma, rep(0, length(q)))
#'           alpha = c(alpha, rep(0, length(q)))
#'           beta = c(beta, beta_obs)
#'
#'           order = rep(q, length(level) + 1)[under_max_alpha]
#'
#'           beta = data.frame(Estimate=beta, order)
#'
#'           C = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
#'           U = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
#'           V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
#'           S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
#'
#'           beta = beta$Estimate
#'
#'           cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
#'
#'           # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
#'         }) %>% abind(along = 3) %>% apply(1:2, sd)
#'         # end = Sys.time()
#'         # end - start
#'
#'         # stopCluster(cl)
#'         # plan(sequential)
#'
#'       } else {
#'
#'         se = matrix(0, ncol = 7, nrow = nrow(beta))
#'         colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
#'         se = as.data.frame(se)
#'
#'       }
#'
#'     }
#'
#'     if (diversity == 'PD') {
#'
#'       if (datatype == 'abundance') {
#'
#'         aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
#'         aL$treeNabu$branch.length = aL$BLbyT[,1]
#'         aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'
#'         gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>%
#'           set_colnames(q) %>% gather(Order, Estimate) %>%
#'           mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Ind(data_gamma, m_gamma), length(q)), Size = rep(m_gamma, length(q)))
#'
#'
#'         aL_table_alpha = c()
#'
#'         for (i in 1:N){
#'
#'           x = data[data[,i]>0,i]
#'           names(x) = rownames(data)[data[,i]>0]
#'
#'           aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
#'           aL$treeNabu$branch.length = aL$BLbyT[,1]
#'           aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'
#'           aL_table_alpha = rbind(aL_table_alpha, aL_table)
#'
#'         }
#'
#'
#'         qPDm = iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")
#'         qPDm = qPDm/N
#'         alpha = qPDm %>% t %>% as.data.frame %>%
#'           set_colnames(q) %>% gather(Order, Estimate) %>%
#'           mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Ind(data_alpha, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'       }
#'
#'       if (datatype == 'incidence_raw') {
#'
#'         aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
#'         aL$treeNabu$branch.length = aL$BLbyT[,1]
#'         aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'         gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>%
#'           set_colnames(q) %>% gather(Order, Estimate) %>%
#'           mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Sam(data_gamma_freq, m_gamma), length(q)), Size = rep(m_gamma, length(q)))
#'
#'         aL_table_alpha = c()
#'
#'         for (i in 1:N){
#'
#'           x = data[[i]]
#'
#'           aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
#'           aL$treeNabu$branch.length = aL$BLbyT[,1]
#'           aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'
#'           aL_table_alpha = rbind(aL_table_alpha, aL_table)
#'
#'         }
#'
#'         alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")/N) %>% t %>% as.data.frame %>%
#'           set_colnames(q) %>% gather(Order, Estimate) %>%
#'           mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Sam(data_alpha_freq, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'
#'       }
#'
#'       gamma = (gamma %>%
#'                  mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
#'         set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
#'
#'       if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order == 0) & (gamma$level > ref_alpha_max)) else under_max_alpha = gamma$level>0
#'       gamma = gamma[under_max_alpha,]
#'       gamma$Order = as.numeric(gamma$Order)
#'
#'
#'       alpha = (alpha %>%
#'                  mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
#'         set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
#'
#'       alpha = alpha[under_max_alpha,]
#'       alpha$Order = as.numeric(alpha$Order)
#'
#'       if (PDtype == 'meanPD') {
#'         gamma$Estimate = gamma$Estimate/reft
#'         alpha$Estimate = alpha$Estimate/reft
#'       }
#'
#'       beta = alpha
#'       beta$Estimate = gamma$Estimate/alpha$Estimate
#'       beta[beta == "Observed"] = "Observed_alpha"
#'       beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate),
#'                                      Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))
#'
#'       C = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(1 - Order) - 1)/(N^(1 - Order) - 1)))
#'       U = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(Order - 1) - 1)/(N^(Order - 1) - 1)))
#'       V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
#'       S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
#'
#'       if(nboot>1){
#'
#'         # cl = makeCluster(cluster_numbers)
#'         # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
#'         #                     'datatype', 'data_2D'))
#'         # clusterEvalQ(cl, library(tidyverse, magrittr))
#'
#'         # plan(sequential)
#'         # plan(multiprocess)
#'
#'         # se = parSapply(cl, 1:nboot, function(i){
#'
#'         # start = Sys.time()
#'         se = future_lapply(1:nboot, function(i){
#'
#'           if (datatype == 'abundance') {
#'
#'             tree_bt = PDtree
#'
#'             bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
#'             p_bt = bootstrap_population
#'             unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))
#'
#'             if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){
#'
#'               unseen = unseen_p[which(rowSums(unseen_p) > 0),]
#'               unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
#'               p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
#'               unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
#'               rownames(p_bt) = c(rownames(data), unseen_name)
#'
#'               bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
#'               x_bt = bootstrap_sample
#'
#'               rownames(x_bt) = rownames(p_bt)
#'
#'               if ( sum(x_bt[-(1:nrow(data)),])>0 ){
#'
#'                 g0_hat = apply(data, 2, function(x){
#'
#'                   n = sum(x)
#'                   f1 = sum(x == 1)
#'                   f2 = sum(x == 2)
#'
#'                   aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
#'
#'                   aL$treeNabu$branch.length = aL$BLbyT[,1]
#'                   aL = aL$treeNabu %>% select(branch.abun,branch.length)
#'                   g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
#'                   g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
#'                   g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
#'                   g0_hat
#'
#'                 })
#'
#'                 te = (x_bt[1:nrow(data),]*(data == 0))>0
#'                 used_length = sapply(1:ncol(data), function(i) {
#'
#'                   if (sum(te[,i]) == 0) return(0) else {
#'
#'                     iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
#'                       subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
#'
#'                   }
#'
#'                 })
#'
#'                 g0_hat = g0_hat - used_length
#'                 g0_hat[g0_hat < 0] = 0
#'
#'                 unseen_sample = x_bt[-(1:nrow(data)),]
#'                 if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt), byrow = T)
#'
#'                 L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )
#'
#'                 L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
#'                 L0_hat[which(rowSums(unseen_sample) == 0)] = 0
#'
#'                 for (i in 1:length(L0_hat)){
#'
#'                   tip = list(edge = matrix(c(2,1),1,2),
#'                              tip.label = unseen_name[i],
#'                              edge.length = L0_hat[i],
#'                              Nnode = 1)
#'                   class(tip) = "phylo"
#'
#'                   tree_bt = tree_bt + tip
#'
#'                 }
#'
#'               } else {
#'
#'                 x_bt = x_bt[1:nrow(data),]
#'                 p_bt = p_bt[1:nrow(data),]
#'
#'               }
#'
#'             } else {
#'
#'               p_bt = p_bt[1:nrow(data),]
#'               x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
#'               rownames(x_bt) = rownames(data)
#'
#'             }
#'
#'             bootstrap_data_gamma = rowSums(x_bt)
#'             bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
#'             bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
#'             bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]
#'
#'             m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma, i, datatype='abundance'))
#'             m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha, i, datatype='abundance'))
#'
#'             aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
#'             aL$treeNabu$branch.length = aL$BLbyT[,1]
#'             aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'
#'             gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t)
#'
#'
#'             aL_table_alpha = c()
#'
#'             for (i in 1:N){
#'
#'               # x = x_bt[x_bt[,i]>0,i]
#'               # names(x) = rownames(p_bt)[x_bt[,i]>0]
#'
#'               x = x_bt[,i]
#'               names(x) = rownames(p_bt)
#'               x = x[x_bt[,i]>0]
#'
#'               aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
#'               aL$treeNabu$branch.length = aL$BLbyT[,1]
#'               aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'
#'               aL_table_alpha = rbind(aL_table_alpha, aL_table)
#'
#'             }
#'
#'             alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")/N) %>% t)
#'
#'             beta_obs = (iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, nt = n, cal = "PD") /
#'                           (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, nt = n, cal = "PD") / N)) %>% unlist()
#'           }
#'
#'           if (datatype == 'incidence_raw') {
#'
#'             tree_bt = PDtree
#'
#'             bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
#'             p_bt = bootstrap_population
#'             unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)
#'
#'             if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){
#'
#'               unseen = unseen_p[which(rowSums(unseen_p)>0),]
#'               unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
#'               p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
#'               unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
#'               rownames(p_bt) = c(rownames(data[[1]]), unseen_name)
#'
#'               raw = lapply(1:ncol(p_bt), function(j){
#'
#'                 lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)
#'
#'               })
#'
#'               for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
#'
#'               if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){
#'
#'                 R0_hat = sapply(data, function(x){
#'
#'                   nT = ncol(x)
#'                   Q1 = sum(rowSums(x)==1)
#'                   Q2 = sum(rowSums(x)==2)
#'
#'                   aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
#'
#'                   aL$treeNabu$branch.length = aL$BLbyT[,1]
#'                   aL = aL$treeNabu %>% select(branch.abun,branch.length)
#'                   R1 = aL$branch.length[aL$branch.abun == 1] %>% sum
#'                   R2 = aL$branch.length[aL$branch.abun == 2] %>% sum
#'                   R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
#'                   R0_hat
#'
#'                 })
#'
#'                 te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums) == 0)) > 0
#'                 used_length = sapply(1:N, function(i) {
#'
#'                   if (sum(te[,i]) == 0) return(0) else {
#'
#'                     iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
#'                       subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
#'
#'                   }
#'
#'                 })
#'
#'                 R0_hat = R0_hat - used_length
#'                 R0_hat[R0_hat < 0] = 0
#'
#'                 unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
#'                 if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = N, byrow = T)
#'
#'                 L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (R0_hat[i] / nrow(unseen)) else 0 )
#'
#'                 L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
#'                 L0_hat[which(rowSums(unseen_sample) == 0)] = 0
#'
#'                 for (i in 1:length(L0_hat)){
#'
#'                   tip = list(edge = matrix(c(2,1), 1, 2),
#'                              tip.label = unseen_name[i],
#'                              edge.length = L0_hat[i],
#'                              Nnode = 1)
#'                   class(tip) = "phylo"
#'
#'                   tree_bt = tree_bt + tip
#'
#'                 }
#'
#'               } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])
#'
#'             } else {
#'
#'               p_bt = p_bt[1:nrow(data[[1]]),]
#'               raw = lapply(1:ncol(p_bt), function(j){
#'
#'                 lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
#'
#'               })
#'
#'               for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
#'
#'             }
#'
#'             gamma = Reduce('+', raw)
#'             gamma[gamma > 1] = 1
#'             bootstrap_data_gamma_raw = gamma
#'             bootstrap_data_gamma_freq = c(n, rowSums(gamma))
#'
#'             bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
#'
#'             bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
#'             bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
#'
#'             m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma_freq, i, datatype = 'incidence'))
#'             m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha_freq, i, datatype = 'incidence'))
#'
#'             aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
#'             aL$treeNabu$branch.length = aL$BLbyT[,1]
#'             aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'
#'             gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal="PD") %>% t)
#'
#'
#'             aL_table_alpha = c()
#'
#'             for (i in 1:N){
#'
#'               x = raw[[i]]
#'
#'               aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
#'               aL$treeNabu$branch.length = aL$BLbyT[,1]
#'               aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#'
#'               aL_table_alpha = rbind(aL_table_alpha, aL_table)
#'
#'             }
#'
#'             alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q=q, nt = n, cal = "PD")/N) %>% t)
#'
#'             beta_obs = (iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, nt = n, cal = "PD") /
#'                           (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, nt = n, cal = "PD") / N)) %>% unlist()
#'           }
#'
#'           gamma = gamma[under_max_alpha]
#'
#'           alpha = alpha[under_max_alpha]
#'
#'           if (PDtype == 'meanPD') {
#'             gamma = gamma/reft
#'             alpha = alpha/reft
#'           }
#'
#'           beta = gamma/alpha
#'
#'           gamma = c(gamma, rep(0, length(q)))
#'           alpha = c(alpha, rep(0, length(q)))
#'           beta = c(beta, beta_obs)
#'
#'           order = rep(q, each = length(level) + 1)[under_max_alpha]
#'
#'           beta = data.frame(Estimate = beta, order)
#'
#'           C = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
#'           U = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
#'           V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
#'           S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
#'
#'           beta = beta$Estimate
#'
#'           cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
#'
#'           # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
#'         }) %>% abind(along = 3) %>% apply(1:2, sd)
#'         # end = Sys.time()
#'         # end - start
#'
#'         # stopCluster(cl)
#'         # plan(sequential)
#'
#'       } else {
#'
#'         se = matrix(0, ncol = 7, nrow = nrow(beta))
#'         colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
#'         se = as.data.frame(se)
#'
#'       }
#'
#'     }
#'
#'     if (diversity == 'FD') {
#'
#'       FD_by_tau = function(data, distM, tau, level, datatype, by, m_gamma, m_alpha) {
#'
#'         if (datatype == 'abundance') {
#'
#'           zik = data
#'           zik = zik[rowSums(data)>0,]
#'
#'           dij = distM
#'           dij = dij[rowSums(data)>0, rowSums(data)>0]
#'
#'           dij[which(dij>tau, arr.ind = T)] = tau
#'           aik = (1 - dij/tau) %*% as.matrix(zik)
#'           positive_id = rowSums(aik)>0
#'
#'
#'           gamma_x = rowSums(zik)[positive_id]
#'           gamma_a = rowSums(aik)[positive_id]
#'           gamma_v = gamma_x/gamma_a
#'           # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
#'
#'           ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
#'
#'           gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
#'           # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector
#'
#'
#'           alpha_x = as.vector(as.matrix(zik))
#'           alpha_a = as.vector(aik)
#'           # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
#'
#'           alpha_v = alpha_x/alpha_a
#'           alpha_v = alpha_v[alpha_a > 0]
#'           alpha_a = alpha_a[alpha_a > 0]
#'
#'           ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
#'
#'           if (by == 'size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
#'           if (by == 'coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
#'           # if (by == 'size') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
#'           # if (by == 'coverage') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
#'
#'         }
#'
#'         if (datatype == 'incidence_raw') {
#'
#'           data_gamma_freq = data$data_gamma_freq
#'           data_2D = data$data_2D
#'
#'           gamma_Y = data_gamma_freq[-1]
#'
#'           dij = distM
#'           dij = dij[gamma_Y > 0, gamma_Y > 0]
#'           gamma_Y = gamma_Y[gamma_Y > 0]
#'
#'           dij[which(dij>tau, arr.ind = T)] = tau
#'           gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
#'           gamma_a[gamma_a > n] = n
#'
#'           gamma_v = gamma_Y/gamma_a
#'
#'           # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
#'           # gamma_a[gamma_a > n] = n
#'
#'           ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
#'
#'           gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
#'           # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector
#'
#'
#'           alpha_Y = data_2D[-1,]
#'
#'           dij = distM
#'           dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
#'           alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]
#'
#'           dij[which(dij>tau, arr.ind = T)] = tau
#'           alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)
#'
#'           # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
#'           alpha_a[alpha_a > n] = n
#'           alpha_a = as.vector(alpha_a)
#'
#'           # alpha_v = rep(gamma_v, N)
#'           # alpha_v = alpha_v[alpha_a > 0]
#'           alpha_v = as.vector(as.matrix(alpha_Y))/alpha_a
#'           alpha_v = alpha_v[alpha_a > 0]
#'           alpha_a = alpha_a[alpha_a > 0]
#'
#'           ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
#'
#'           if (by == 'size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
#'           if (by == 'coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
#'           # if (by == 'size') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
#'           # if (by == 'coverage') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
#'
#'         }
#'
#'         return(data.frame(gamma,alpha))
#'
#'       }
#'
#'       if (FDtype == 'single'){
#'
#'         if (datatype == 'abundance') {
#'
#'           output = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
#'           gamma = output$gamma
#'           alpha = output$alpha
#'
#'           gamma = data.frame(level, gamma) %>%
#'             mutate(Method = ifelse(level>=ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each=length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Ind(data_gamma, m_gamma), length(q)), Size = rep(m_gamma, length(q)))
#'
#'           alpha = data.frame(level, alpha) %>%
#'             mutate(Method = ifelse(level>=ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each=length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Ind(data_alpha, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'           beta = alpha
#'           beta$alpha = gamma$gamma/alpha$alpha
#'
#'           ## Observed Beta ##
#'           output_obs = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))
#'           gamma_obs = output_obs$gamma
#'           alpha_obs = output_obs$alpha
#'           obs_beta = gamma_obs/alpha_obs
#'
#'         }
#'
#'         if (datatype == 'incidence_raw') {
#'
#'           output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype='incidence_raw', by = 'coverage', m_gamma=m_gamma, m_alpha=m_alpha)
#'           gamma = output$gamma
#'           alpha = output$alpha
#'
#'           gamma = data.frame(level, gamma) %>%
#'             mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Sam(data_gamma_freq, m_gamma), length(q)), Size = rep(m_gamma, length(q)))
#'
#'           alpha = data.frame(level, alpha) %>%
#'             mutate(Method = ifelse(level>=ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Sam(data_alpha_freq, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'           beta = alpha
#'           beta$alpha = gamma$gamma/alpha$alpha
#'
#'           ## Observed Beta ##
#'           output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_gamma_freq[1])
#'           gamma_obs = output_obs$gamma
#'           alpha_obs = output_obs$alpha
#'           obs_beta = gamma_obs/alpha_obs
#'
#'         }
#'
#'       }
#'
#'       if (FDtype == 'AUC'){
#'
#'         cut = seq(0.00000001, 1, length.out = FDcut_number)
#'         width = diff(cut)
#'
#'         if (datatype == 'abundance') {
#'
#'           gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'             FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
#'
#'           })
#'
#'           gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
#'
#'           left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
#'           right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
#'
#'           gamma = colSums((left_limit + right_limit)/2)
#'
#'           alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
#'
#'           left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
#'           right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
#'
#'           alpha = colSums((left_limit + right_limit)/2)
#'
#'           beta_over_tau = gamma_over_tau/alpha_over_tau
#'
#'           left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
#'           right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
#'
#'           beta = colSums((left_limit + right_limit)/2)
#'
#'           gamma = data.frame(level, gamma) %>%
#'             mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Ind(data_gamma, m_gamma), length(q)), Size = rep(m_gamma, length(q)))
#'
#'           alpha = data.frame(level, alpha) %>%
#'             mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Ind(data_alpha, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'           beta = data.frame(level, beta) %>%
#'             mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Ind(data_alpha, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'           ## Observed Beta ##
#'           obs_gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'             FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))
#'
#'           })
#'
#'           obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)
#'
#'           obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)
#'
#'         }
#'
#'         if (datatype == 'incidence_raw') {
#'
#'           gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'             FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
#'
#'           })
#'
#'           gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
#'
#'           left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
#'           right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
#'
#'           gamma = colSums((left_limit + right_limit)/2)
#'
#'           alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
#'
#'           left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
#'           right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
#'
#'           alpha = colSums((left_limit + right_limit)/2)
#'
#'           beta_over_tau = gamma_over_tau/alpha_over_tau
#'
#'           left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
#'           right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
#'
#'           beta = colSums((left_limit + right_limit)/2)
#'
#'           gamma = data.frame(level, gamma) %>%
#'             mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Sam(data_gamma_freq, m_gamma), length(q)), Size = rep(m_gamma, length(q)))
#'
#'           alpha = data.frame(level, alpha) %>%
#'             mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Sam(data_alpha_freq, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'           beta = data.frame(level, beta) %>%
#'             mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
#'                    Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Chat.Sam(data_alpha_freq, m_alpha), length(q)), Size = rep(m_alpha, length(q)))
#'
#'           ## Observed Beta ##
#'           obs_gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'             FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_alpha_freq[1])
#'
#'           })
#'
#'           obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)
#'
#'           obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)
#'
#'         }
#'
#'       }
#'
#'       gamma = gamma[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
#'
#'       if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order == 0) & (gamma$level > ref_alpha_max)) else under_max_alpha = gamma$level > 0
#'       gamma = gamma[under_max_alpha,]
#'
#'
#'       alpha = alpha[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
#'
#'       alpha = alpha[under_max_alpha,]
#'
#'       beta = beta[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
#'
#'       beta = beta[under_max_alpha,]
#'
#'       beta[beta == "Observed"] = "Observed_alpha"
#'       beta = beta %>%
#'         rbind(., data.frame(Estimate = obs_beta, Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))
#'
#'       C = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(1 - Order) - 1)/(N^(1 - Order) - 1)))
#'       U = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(Order - 1) - 1)/(N^(Order - 1) - 1)))
#'       V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
#'       S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
#'
#'       if(nboot > 1){
#'
#'         # cl = makeCluster(cluster_numbers)
#'         # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
#'         #                     'datatype', 'data_2D'))
#'         # clusterEvalQ(cl, library(tidyverse, magrittr))
#'
#'         # plan(sequential)
#'         # plan(multiprocess, workers=7)
#'
#'         # se = parSapply(cl, 1:nboot, function(i){
#'
#'         # start = Sys.time()
#'         se = future_lapply(1:nboot, function(i){
#'
#'           if (datatype == 'abundance') {
#'
#'             p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
#'             f0_hat = nrow(p_bt) - nrow(data)
#'
#'             distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')
#'
#'             data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
#'
#'             data_gamma = rowSums(data_bt)
#'             data_gamma = data_gamma[data_gamma>0]
#'             data_alpha = as.matrix(data_bt) %>% as.vector
#'
#'             if (FDtype == 'single'){
#'
#'               output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype='abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
#'               gamma = output$gamma
#'               alpha = output$alpha
#'
#'               beta=gamma/alpha
#'
#'               ## Observed Beta ##
#'               output_obs = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))
#'               gamma_obs = output_obs$gamma
#'               alpha_obs = output_obs$alpha
#'               beta_obs = gamma_obs/alpha_obs
#'
#'             }
#'
#'             if (FDtype == 'AUC'){
#'
#'               gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'                 FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
#'
#'               })
#'
#'               gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
#'
#'               left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
#'
#'               gamma = colSums((left_limit + right_limit)/2)
#'
#'               alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
#'
#'               left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
#'
#'               alpha = colSums((left_limit + right_limit)/2)
#'
#'               beta_over_tau = gamma_over_tau/alpha_over_tau
#'
#'               left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
#'
#'               beta = colSums((left_limit + right_limit)/2)
#'
#'               ## Observed Beta ##
#'               gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'                 FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))
#'
#'               })
#'
#'               gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
#'
#'               alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
#'
#'               beta_over_tau = gamma_over_tau/alpha_over_tau
#'
#'               left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
#'
#'               beta_obs = colSums((left_limit + right_limit)/2)
#'
#'             }
#'
#'           }
#'
#'           if (datatype == 'incidence_raw') {
#'
#'             p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
#'             f0_hat = nrow(p_bt) - nrow(data_2D[-1,])
#'
#'             distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')
#'
#'             # p_bt = p_bt[rowSums(p_bt)>0,]
#'
#'             raw = lapply(1:ncol(p_bt), function(j){
#'
#'               lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
#'
#'             })
#'
#'             gamma = Reduce('+', raw)
#'             gamma[gamma>1] = 1
#'             data_gamma_raw_bt = gamma
#'             data_gamma_freq_bt = c(n, rowSums(gamma))
#'
#'             data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)
#'
#'             # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt > 0]
#'             # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt > 0]
#'
#'             data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
#'
#'             if (FDtype == 'single'){
#'
#'               output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
#'               gamma = output$gamma
#'               alpha = output$alpha
#'
#'               beta = gamma/alpha
#'
#'               ## Observed Beta ##
#'               output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])
#'               gamma_obs = output_obs$gamma
#'               alpha_obs = output_obs$alpha
#'               beta_obs = gamma_obs/alpha_obs
#'
#'             }
#'
#'             if (FDtype == 'AUC'){
#'
#'               gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'                 FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
#'
#'               })
#'
#'               gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
#'
#'               left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
#'
#'               gamma = colSums((left_limit + right_limit)/2)
#'
#'               alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
#'
#'               left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
#'
#'               alpha = colSums((left_limit + right_limit)/2)
#'
#'               beta_over_tau = gamma_over_tau/alpha_over_tau
#'
#'               left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
#'
#'               beta = colSums((left_limit + right_limit)/2)
#'
#'               ## Observed Beta ##
#'               gamma_alpha_over_tau = lapply(cut, function(tau) {
#'
#'                 FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])
#'
#'               })
#'
#'               gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
#'
#'               alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
#'
#'               beta_over_tau = gamma_over_tau/alpha_over_tau
#'
#'               left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
#'               right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
#'
#'               beta_obs = colSums((left_limit + right_limit)/2)
#'
#'             }
#'
#'           }
#'
#'           gamma = gamma[under_max_alpha]
#'           alpha = alpha[under_max_alpha]
#'           beta = beta[under_max_alpha]
#'
#'           beta = gamma/alpha
#'
#'           gamma = c(gamma, rep(0, length(q)))
#'           alpha = c(alpha, rep(0, length(q)))
#'           beta = c(beta, beta_obs)
#'
#'           order = rep(q, each=length(level) + 1)[under_max_alpha]
#'
#'           beta = data.frame(Estimate = beta, order)
#'
#'           C = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
#'           U = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
#'           V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
#'           S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
#'
#'           beta = beta$Estimate
#'
#'           cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
#'
#'           # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
#'         }) %>% abind(along = 3) %>% apply(1:2, sd)
#'         # end = Sys.time()
#'         # end - start
#'
#'         # stopCluster(cl)
#'         # plan(sequential)
#'
#'       } else {
#'
#'         se = matrix(0, ncol = 7, nrow = nrow(beta))
#'         colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
#'         se = as.data.frame(se)
#'
#'       }
#'
#'     }
#'
#'
#'     se = as.data.frame(se)
#'
#'     gamma = gamma %>% mutate(s.e. = se$gamma[1:(nrow(se) - length(q))],
#'                              LCL = Estimate - tmp * se$gamma[1:(nrow(se) - length(q))],
#'                              UCL = Estimate + tmp * se$gamma[1:(nrow(se) - length(q))],
#'                              Region = region_name)
#'
#'     alpha = alpha %>% mutate(s.e. = se$alpha[1:(nrow(se) - length(q))],
#'                              LCL = Estimate - tmp * se$alpha[1:(nrow(se) - length(q))],
#'                              UCL = Estimate + tmp * se$alpha[1:(nrow(se) - length(q))],
#'                              Region = region_name)
#'
#'     beta = beta %>% mutate(  s.e. = se$beta,
#'                              LCL = Estimate - tmp * se$beta,
#'                              UCL = Estimate + tmp * se$beta,
#'                              Region = region_name)
#'
#'     C = C %>% mutate(        s.e. = se$C,
#'                              LCL = Estimate - tmp * se$C,
#'                              UCL = Estimate + tmp * se$C,
#'                              Region = region_name)
#'
#'
#'     U = U %>% mutate(        s.e. = se$U,
#'                              LCL = Estimate - tmp * se$U,
#'                              UCL = Estimate + tmp * se$U,
#'                              Region = region_name)
#'
#'     V = V %>% mutate(        s.e. = se$V,
#'                              LCL = Estimate - tmp * se$V,
#'                              UCL = Estimate + tmp * se$V,
#'                              Region = region_name)
#'
#'     S = S %>% mutate(        s.e. = se$S,
#'                              LCL = Estimate - tmp * se$S,
#'                              UCL = Estimate + tmp * se$S,
#'                              Region = region_name)
#'
#'     list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)
#'
#'   }
#'
#'   output = lapply(1:length(data_list), function(i) for_each_region(data = data_list[[i]], region_name = region_names[i], N = Ns[i]))
#'   names(output) = region_names
#'
#'   return(output)
#'
#' }
#'
#'
#'
#' #' ggplot for Beta diversity
#' #'
#' #' \code{ggiNEXTbeta}: ggplot for Interpolation and extrapolation of Beta diversity with order q
#' #'
#' #' @param output the output from iNEXTbeta
#' #' @param type selection of plot type : \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ; \code{type = 'D'} for plotting 4 turnover dissimilarities.
#' #' @param measurement character indicating the label of y-axis.
#' #' @param scale Are scales shared across all facets (\code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})? Default is \code{"free"}.
#' #' @param main The title of the plot.
#' #' @param transp a value between 0 and 1 controlling transparency. \code{transp = 0} is completely transparent, default is 0.4.
#' #'
#' #' @return a figure for Beta diversity or dissimilarity diversity.
#' #'
#' #' @examples
#' #' ## Taxonomic diversity for abundance data
#' #' data(ant)
#' #' output1 = iNEXTbeta(data = ant, diversity = 'TD', datatype = 'abundance', level = seq(0.5, 1, 0.05),
#' #'                     nboot = 20)
#' #'
#' #' ggiNEXTbeta(output1, type = 'B', measurement = 'TD', scale = 'free', main = NULL, transp = 0.4)
#' #' ggiNEXTbeta(output1, type = 'D', measurement = 'TD', scale = 'free', main = NULL, transp = 0.4)
#' #'
#' #'
#' #' ## Taxonomic diversity for incidence data
#' #' data(fish)
#' #' output2 = iNEXTbeta(data = fish, diversity = 'TD', datatype = 'incidence_raw', level = seq(0.7, 0.95, 0.05),
#' #'                     nboot = 20, conf = 0.95)
#' #'
#' #' ggiNEXTbeta(output2, type = 'B', measurement = 'TD', scale = 'free', main = NULL, transp = 0.4)
#' #' ggiNEXTbeta(output2, type = 'D', measurement = 'TD', scale = 'free', main = NULL, transp = 0.4)
#' #'
#' #'
#' #' ## Phylogenetic Hill numbers for abundance data
#' #' data(ant)
#' #' data(ant_tree)
#' #' output3 = iNEXTbeta(data = ant, diversity = 'PD', datatype = 'abundance',  level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 20, conf = 0.95, PDtree = ant_tree, PDreftime = NULL, PDtype = 'meanPD')
#' #'
#' #' ggiNEXTbeta(output3, type = 'B', measurement = 'meanPD', scale = 'free', main = NULL, transp = 0.4)
#' #' ggiNEXTbeta(output3, type = 'D', measurement = 'meanPD', scale = 'free', main = NULL, transp = 0.4)
#' #'
#' #'
#' #' ## Phylogenetic Hill numbers for incidence data
#' #' data(fish)
#' #' data(fish_tree)
#' #' output4 = iNEXTbeta(data = fish, diversity = 'PD', datatype = 'incidence_raw', level = seq(0.7, 0.95, 0.05),
#' #'                     nboot = 0, conf = 0.95, PDtree = fish_tree, PDreftime = NULL, PDtype = 'meanPD')
#' #'
#' #' ggiNEXTbeta(output4, type = 'B', measurement = 'meanPD', scale = 'free', main = NULL, transp=0.4)
#' #' ggiNEXTbeta(output4, type = 'D', measurement = 'meanPD', scale = 'free', main = NULL, transp= 0.4)
#' #'
#' #'
#' #' ## Functional diversity for abundance data under single threshold
#' #' data(ant)
#' #' data(ant_distM)
#' #' output5 = iNEXTbeta(data = ant, diversity = 'FD', datatype = 'abundance', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 20, conf = 0.95, FDdistM = ant_distM, FDtype = 'single', FDtau = NULL)
#' #'
#' #' ggiNEXTbeta(output5, type = 'B', measurement = 'FD_tau', scale = 'free', main = NULL, transp = 0.4)
#' #' ggiNEXTbeta(output5, type = 'D', measurement = 'FD_tau', scale = 'free', main = NULL, transp = 0.4)
#' #'
#' #'
#' #' ## Functional diversity for incidence data under single threshold
#' #' data(fish)
#' #' data(fish_distM)
#' #' output6 = iNEXTbeta(data = fish, diversity = 'FD', datatype = 'incidence_raw', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 20, conf = 0.95, FDdistM = fish_distM, FDtype = 'single', FDtau = NULL)
#' #'
#' #' ggiNEXTbeta(output6, type = 'B', measurement = 'FD_tau', scale = 'free', main = NULL, transp = 0.4)
#' #' ggiNEXTbeta(output6, type = 'D', measurement = 'FD_tau', scale = 'free', main = NULL, transp = 0.4)
#' #'
#' #'
#' #' #' ## Functional diversity for abundance data with thresholds integrating from 0 to 1
#' #' data(ant)
#' #' data(ant_distM)
#' #' output7 = iNEXTbeta(data = ant, diversity = 'FD', datatype = 'abundance', level = seq(0.5, 0.9, 0.05),
#' #'                     nboot = 10, conf = 0.95, FDdistM = ant_distM, FDtype = 'AUC', FDcut_number = 30)
#' #'
#' #' ggiNEXTbeta(output7, type = 'B', measurement = 'FD_AUC', scale = 'free', main = NULL, transp = 0.4)
#' #' ggiNEXTbeta(output7, type = 'D', measurement = 'FD_AUC', scale = 'free', main = NULL, transp = 0.4)
#' #'
#' #'
#' #' ## Functional diversity for incidence data with thresholds integrating from 0 to 1
#' #' data(fish)
#' #' data(fish_distM)
#' #' output8 = iNEXTbeta(data = fish, diversity = 'FD', datatype = 'incidence_raw', level = seq(0.5, 0.9, 0.05),
#' #'                    nboot = 10, conf = 0.95, FDdistM = fish_distM, FDtype = 'AUC', FDcut_number = 30)
#' #'
#' #' ggiNEXTbeta(output8, type = 'B', measurement = 'FD_AUC', scale = 'free', main = NULL, transp = 0.4)
#' #' ggiNEXTbeta(output8, type = 'D', measurement = 'FD_AUC', scale = 'free', main = NULL, transp = 0.4)
#' #'
#' #' @export
#' ggiNEXTbeta = function(output, type = c('B', 'D'),
#'                        measurement = c('TD', 'PD','meanPD' ,'FD_tau', 'FD_AUC'),
#'                        scale = 'free', main = NULL, transp = 0.4){
#'
#'   if (type == 'B'){
#'
#'     gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
#'     alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
#'     beta =  lapply(output, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()
#'     beta = beta %>% filter(Method != 'Observed')
#'     beta[beta == 'Observed_alpha'] = 'Observed'
#'
#'     # Dropping out the points extrapolated over double reference size
#'     gamma1 = data.frame() ; alpha1 = data.frame() ; beta1 = data.frame()
#'
#'     for(i in 1:length(unique(gamma$Region))){
#'
#'       Gamma <- gamma %>% filter(Region==unique(gamma$Region)[i]) ; ref_size = unique(Gamma[Gamma$Method=="Observed",]$Size)
#'       Gamma = Gamma %>% filter(!(Order==0 & round(Size)>2*ref_size))
#'
#'       Alpha <- alpha %>% filter(Region==unique(gamma$Region)[i]) ; Alpha = Alpha %>% filter(!(Order==0 & round(Size)>2*ref_size))
#'       Beta <- beta %>% filter(Region==unique(gamma$Region)[i]) ; Beta = Beta %>% filter(!(Order==0 & round(Size)>2*ref_size))
#'
#'       gamma1 = rbind(gamma1,Gamma) ; alpha1 = rbind(alpha1,Alpha) ; beta1 = rbind(beta1,Beta)
#'
#'     }
#'
#'     gamma = gamma1 ; alpha = alpha1 ; beta= beta1
#'
#'     df = rbind(gamma, alpha, beta)
#'     for (i in unique(gamma$Order)) df$Order[df$Order == i] = paste0('q = ', i)
#'     df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))
#'
#'     id_obs = which(df$Method == 'Observed')
#'
#'     for (i in 1:length(id_obs)) {
#'
#'       new = df[id_obs[i],]
#'       new$level = new$level - 0.0001
#'       new$Method = 'Interpolated'
#'
#'       newe = df[id_obs[i],]
#'       newe$level = newe$level + 0.0001
#'       newe$Method = 'Extrapolated'
#'
#'       df = rbind(df, new, newe)
#'
#'     }
#'
#'     if (measurement == 'TD') { ylab = "Taxonomic diversity" }
#'     if (measurement == 'PD') { ylab = "Phylogenetic diversity" }
#'     if (measurement == 'meanPD') { ylab = "Phylogenetic Hill number" }
#'     if (measurement == 'FD_tau') { ylab = "Functional diversity (given tau)" }
#'     if (measurement == 'FD_AUC') { ylab = "Functional diversity (AUC)" }
#'
#'   }
#'
#'   if (type == 'D'){
#'
#'     C = lapply(output, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
#'     U = lapply(output, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
#'     V = lapply(output, function(y) y[["V"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
#'     S = lapply(output, function(y) y[["S"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()
#'     C = C %>% filter(Method != 'Observed')
#'     U = U %>% filter(Method != 'Observed')
#'     V = V %>% filter(Method != 'Observed')
#'     S = S %>% filter(Method != 'Observed')
#'     C[C == 'Observed_alpha'] = U[U == 'Observed_alpha'] = V[V == 'Observed_alpha'] = S[S == 'Observed_alpha'] = 'Observed'
#'
#'     # Dropping out the points extrapolated over double reference size
#'     c1 = data.frame() ; u1 = data.frame() ; v1 = data.frame() ; s1 = data.frame()
#'
#'     for(i in 1:length(unique(C$Region))){
#'
#'       CC <- C %>% filter(Region==unique(C$Region)[i]) ; ref_size = unique(CC[CC$Method=="Observed",]$Size)
#'       CC = CC %>% filter(!(Order==0 & round(Size)>2*ref_size))
#'
#'       UU <- U %>% filter(Region==unique(C$Region)[i]) ; UU = UU %>% filter(!(Order==0 & round(Size)>2*ref_size))
#'       VV <- V %>% filter(Region==unique(C$Region)[i]) ; VV = VV %>% filter(!(Order==0 & round(Size)>2*ref_size))
#'       SS <- S %>% filter(Region==unique(C$Region)[i]) ; SS = SS %>% filter(!(Order==0 & round(Size)>2*ref_size))
#'
#'       c1 = rbind(c1,CC) ; u1 = rbind(u1,UU) ; v1 = rbind(v1,VV) ; s1 = rbind(s1,SS)
#'
#'     }
#'
#'     C = c1 ; U = u1 ; V = v1 ; S = s1
#'
#'
#'     df = rbind(C, U, V, S)
#'     for (i in unique(C$Order)) df$Order[df$Order == i] = paste0('q = ', i)
#'     df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))
#'
#'     id_obs = which(df$Method == 'Observed')
#'
#'     for (i in 1:length(id_obs)) {
#'
#'       new = df[id_obs[i],]
#'       new$level = new$level - 0.0001
#'       new$Method = 'Interpolated'
#'
#'       newe = df[id_obs[i],]
#'       newe$level = newe$level + 0.0001
#'       newe$Method = 'Extrapolated'
#'
#'       df = rbind(df, new, newe)
#'
#'     }
#'
#'     if (measurement == 'TD') { ylab = "Taxonomic dissimilarity" }
#'     if (measurement == 'PD') { ylab = "Phylogenetic dissimilarity" }
#'     if (measurement == 'meanPD') { ylab = "Phylogenetic dissimilarity" }
#'     if (measurement == 'FD_tau') { ylab = "Functional dissimilarity (given tau)" }
#'     if (measurement == 'FD_AUC') { ylab = "Functional dissimilarity (AUC)" }
#'
#'   }
#'
#'   lty = c(Interpolated = "solid", Extrapolated = "dashed")
#'   df$Method = factor(df$Method, levels = c('Interpolated', 'Extrapolated', 'Observed'))
#'
#'   double_size = unique(df[df$Method == "Observed",]$Size)*2
#'   double_extrapolation = df %>% filter(Method == "Extrapolated" & round(Size) %in% double_size)
#'
#'   ggplot(data = df, aes(x = level, y = Estimate, col = Region)) +
#'     geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=transp) +
#'     geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) + scale_linetype_manual(values = lty) +
#'     # geom_line(lty=2) +
#'     geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 3) +
#'     geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 3, stroke = 1.5)+
#'     geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 3) +
#'     geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 3, stroke = 1.5) +
#'     facet_grid(div_type ~ Order, scales = scale) +
#'     theme_bw() +
#'     theme(legend.position = "bottom", legend.title = element_blank()) +
#'     labs(x = 'Sample coverage', y = ylab, title = main)
#' }
#'
#'
#' coverage_to_size = function(x, C, datatype = 'abundance'){
#'
#'   if (datatype == 'abundance'){
#'
#'     n <- sum(x)
#'     refC <- iNEXT.3D:::Chat.Ind(x, n)
#'     f <- function(m, C) abs(iNEXT.3D:::Chat.Ind(x, m) - C)
#'     if (refC == C) {
#'       mm = n
#'     } else if (refC > C) {
#'       opt <- optimize(f, C = C, lower = 0, upper = sum(x))
#'       mm <- opt$minimum
#'       mm <- round(mm)
#'     } else if (refC < C) {
#'       f1 <- sum(x == 1)
#'       f2 <- sum(x == 2)
#'       if (f1 > 0 & f2 > 0) {
#'         A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
#'       }
#'       if (f1 > 1 & f2 == 0) {
#'         A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
#'       }
#'       if (f1 == 1 & f2 == 0) {
#'         A <- 1
#'       }
#'       if (f1 == 0 & f2 == 0) {
#'         A <- 1
#'       }
#'       if (f1 == 0 & f2 > 0) {
#'         A <- 1
#'       }
#'       mm <- (log(n/f1) + log(1 - C))/log(A) - 1
#'       if (is.nan(mm) == TRUE) mm = Inf
#'       mm <- n + mm
#'       mm <- round(mm)
#'     }
#'
#'   } else {
#'
#'     m <- NULL
#'     n <- max(x)
#'     refC <- iNEXT.3D:::Chat.Sam(x, n)
#'     f <- function(m, C) abs(iNEXT.3D:::Chat.Sam(x, m) - C)
#'     if (refC == C) {
#'       mm = n
#'     } else if (refC > C) {
#'       opt <- optimize(f, C = C, lower = 0, upper = max(x))
#'       mm <- opt$minimum
#'       mm <- round(mm)
#'     } else if (refC < C) {
#'       f1 <- sum(x == 1)
#'       f2 <- sum(x == 2)
#'       U <- sum(x) - max(x)
#'       if (f1 > 0 & f2 > 0) {
#'         A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
#'       }
#'       if (f1 > 1 & f2 == 0) {
#'         A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
#'       }
#'       if (f1 == 1 & f2 == 0) {
#'         A <- 1
#'       }
#'       if (f1 == 0 & f2 == 0) {
#'         A <- 1
#'       }
#'       if (f1 == 0 & f2 > 0) {
#'         A <- 1
#'       }
#'       mm <- (log(U/f1) + log(1 - C))/log(A) - 1
#'       if (is.nan(mm) == TRUE) mm = Inf
#'       mm <- n + mm
#'       mm <- round(mm)
#'     }
#'
#'   }
#'
#'   return(mm)
#' }
#'
#' bootstrap_population_multiple_assemblage = function(data, data_gamma, datatype){
#'
#'   if (datatype == 'abundance'){
#'
#'     S_obs = sum(data_gamma > 0)
#'     n = sum(data_gamma)
#'     f1 = sum(data_gamma == 1)
#'     f2 = sum(data_gamma == 2)
#'     f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()
#'
#'     output = apply(data, 2, function(x){
#'
#'       p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)
#'
#'       if(length(p_i_hat) != length(x)){
#'
#'         p_i_hat_unobs = p_i_hat[(length(x)+1):length(p_i_hat)]
#'         p_i_hat_obs = p_i_hat[1:length(x)]
#'         p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
#'         candidate = which(p_i_hat==0)
#'
#'         chosen = sample(x = candidate, size = min(length(p_i_hat_unobs), length(candidate)), replace = F)
#'         p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)
#'
#'         p_i_hat
#'
#'       } else {
#'
#'         p_i_hat = c(p_i_hat, rep(0, f0_hat))
#'         p_i_hat
#'
#'       }
#'     })
#'
#'   }
#'
#'   if (datatype == 'incidence'){
#'
#'     S_obs = sum(data_gamma > 0)
#'     t = data_gamma[1]
#'     Q1 = sum(data_gamma == 1)
#'     Q2 = sum(data_gamma == 2)
#'     Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling
#'
#'     output = apply(data, 2, function(x){
#'
#'       pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)
#'
#'       if(length(pi_i_hat) != (length(x) - 1)){
#'
#'         pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
#'         pi_i_hat_obs = pi_i_hat[1:(length(x)-1)]
#'         pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
#'         candidate = which(pi_i_hat == 0)
#'         chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs), length(candidate)), replace = F)
#'         pi_i_hat[chosen] = pi_i_hat_unobs
#'
#'         pi_i_hat
#'
#'       } else {
#'
#'         pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
#'         pi_i_hat
#'
#'       }
#'     })
#'
#'   }
#'
#'   return(output)
#'
#' }
#'
#' Bootstrap_distance_matrix = function(data, distance_matrix, f0.hat, datatype){
#'
#'   if (datatype == "incidence_freq") {
#'     n = data[1]
#'     X = data[-1]
#'     u = sum(data)
#'   } else if (datatype == "abundance") {
#'     n = sum(data)
#'     X = data
#'   }
#'
#'   # n = sum(data)
#'   distance = as.matrix(distance_matrix)
#'   dij = distance
#'   # X = data
#'
#'   F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
#'   F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
#'
#'   if (datatype == "abundance") {
#'     F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
#'     F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
#'   } else if (datatype == "incidence_freq") {
#'     F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
#'     F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
#'   }
#'
#'   # F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
#'   # F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
#'
#'   if (f0.hat == 0) {
#'     d = dij
#'   } else if (f0.hat == 1) {
#'     random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
#'     d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
#'     d00 = matrix(0, f0.hat, f0.hat)
#'     d <- cbind(dij, d.0bar )
#'     aa <- cbind(t(d.0bar), d00 )
#'     d <- rbind(d, aa)
#'     diag(d) = 0
#'   } else {
#'     random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
#'     d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
#'
#'     fo.num = (f0.hat * (f0.hat-1) )/2
#'     random_d00 = as.vector(rmultinom(1, 1000, rep(1/fo.num, fo.num) ) )/1000
#'     d00 = matrix(0, f0.hat, f0.hat)
#'     d00[upper.tri(d00)] = (F00hat/2)*random_d00
#'     d00 <- pmax(d00, t(d00))###signmatrix
#'     d <- cbind(dij, d.0bar )
#'     aa <- cbind(t(d.0bar), d00 )
#'     d <- rbind(d, aa)
#'     diag(d) = 0
#'   }
#'
#'   return(d)
#'
#' }
#'
#' FD.m.est_0 = function (ai_vi, m, q, nT) {
#'   EFD = function(m, qs, obs, asy, beta, av) {
#'     m = m - nT
#'     out <- sapply(1:length(qs), function(i) {
#'       if (qs[i] != 2) {
#'         obs[i] + (asy[i] - obs[i]) * (1 - (1 - beta[i])^m)
#'       }
#'       else if (qs[i] == 2) {
#'         V_bar^2/sum((av[, 2]) * ((1/(nT + m)) * (av[, 1]/nT) + ((nT + m - 1)/(nT + m)) * (av[, 1] * (av[, 1] - 1)/(nT * (nT - 1)))))
#'       }
#'     })
#'     return(out)
#'   }
#'   V_bar <- sum(ai_vi$ai[, 1] * ai_vi$vi[, 1])/nT
#'   asy <- iNEXT.3D:::FD_est(ai_vi, q, nT)$est
#'   obs <- iNEXT.3D:::FD_mle(ai_vi, q)
#'   out <- sapply(1:ncol(ai_vi$ai), function(i) {
#'     ai <- ai_vi$ai[, i]
#'     ai[ai < 1] <- 1
#'     av = cbind(ai = round(ai), vi = ai_vi$vi[, i])
#'     RFD_m = iNEXT.3D:::RFD(av, nT, nT - 1, q, V_bar)
#'     beta <- rep(0, length(q))
#'     asy_i <- asy[, i]
#'     obs_i <- obs[, i]
#'     asy_i <- sapply(1:length(q), function(j) {
#'       max(asy_i[j], obs_i[j], RFD_m[j])
#'     })
#'
#'     obs_i <- sapply(1:length(q), function(j) {
#'       max(RFD_m[j], obs_i[j])
#'     })
#'
#'     beta0plus <- which(asy_i != obs_i)
#'     beta[beta0plus] <- (obs_i[beta0plus] - RFD_m[beta0plus])/(asy_i[beta0plus] - RFD_m[beta0plus])
#'     sapply(m, function(mm) {
#'       if (mm < nT) {
#'         iNEXT.3D:::RFD(av, nT, mm, q, V_bar)
#'       }
#'       else if (mm == nT) {
#'         obs_i
#'       }
#'       else if (mm == Inf) {
#'         asy_i
#'       }
#'       else {
#'         EFD(m = mm, qs = q, obs = obs_i, asy = asy_i,
#'             beta = beta, av = av)
#'       }
#'     }) %>% t %>% as.numeric
#'   })
#'   matrix(out, ncol = ncol(ai_vi$ai))
#' }
