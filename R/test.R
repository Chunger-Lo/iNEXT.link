# library(ggplot2)
# library(knitr)
# library(tibble)
# library(iNEXT.3D)
# library(iNEXT.beta)
# library(iNEXT.link)
# library(ade4)
# library(sets)
# library(tidyr)
# library(tidyverse)
# library(ape)
# library(dplyr)
# library(phytools)
# library(magrittr)
# library(chaoUtility)
# library(future.apply)
# library(abind)
#
# # source('R//iNEXTlink.R')
# # source('R//myfunc.R')
# # source('R//subfunction.R')
# # source('R//crazy_min.R')
#
# data("puerto.rico")
# network_list = puerto.rico$data
# coltree = puerto.rico$col.tree
# rowtree = puerto.rico$row.tree
#
# data("beetles_pool")
# network_list = beetles_pool$data
# coltree = beetles_pool$col.tree
# rowtree = NULL
# #
# #
# datainf1 = DataInfo.link(data = network_list,
#                          diversity = 'TD', datatype = "abundance")
# ## taxonomic ----
# asy1 = Asy.link(network_list, diversity = 'TD', q = c(0,1,2),nboot = 10)
# obs1 = Obs.link(network_list, diversity = 'TD', q = c(0,1,2),nboot = 10)
# sc1 <- SC.link(x = network_list, datatype = "abundance", nboot = 10)
# spec1 = Spec.link(network_list, q = seq(0,2,1), datatype = "abundance", nboot =10)
#
# dissimilarity1 = iNEXTbeta.link(networks = network_list, level = seq(0.5,1,0.05),
#                                 datatype='abundance',
#                                 q = c(0, 1, 2),diversity='TD',
#                                 nboot = 10, conf = 0.95)
# # ## phylogenetic ----
# b = 5
# sample2 <- iNEXT.link::iNEXT.link(data = network_list, diversity = 'PD', datatype = "abundance",
#                       nboot = b,col.tree = coltree, row.tree = rowtree)
#
#
# asy2 = iNEXT.link::Asy.link(data = network_list, diversity = 'PD', q = seq(0,2,0.25),nboot = b,
#                 col.tree = coltree, row.tree = rowtree)
# obs2 = iNEXT.link::Obs.link(data = network_list, diversity = 'PD', q = seq(0,2,0.25),nboot = b,
#                 col.tree = coltree, row.tree = rowtree)
# dissimilarity2 = iNEXT.link::iNEXTbeta.link(networks = network_list, level = seq(0.5,1,0.2),
#                                  datatype='abundance',
#                                  q = c(0, 1, 2),diversity='PD',
#                                  nboot = 5, conf = 0.95,
#                                  col.tree = coltree, row.tree = rowtree)
