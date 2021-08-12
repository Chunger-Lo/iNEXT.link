# # # devtools::install_github('crazymin2266/iNEXT.beta')
# # # devtools::install_github('KaiHsiangHu/iNEXT.3D')
# # # devtools::install_github('KaiHsiangHu/iNEXTPD2')
# # # devtools::install_github('KaiHsiangHu/iNEXT.4steps')
# # # devtools::install_github('chaolab2019/chaoUtility')
# # # devtools::install_github('Chunger-Lo/iNEXT.link')
# # # devtools::install_github('YanHanChen/PhD')
# # # #
# library(ggplot2)
# library(knitr)
# library(tibble)
# library(iNEXT.3D)
# library(iNEXT.beta)
# library(iNEXT.4steps)
# library(ade4)
# library(tidyr)
# library(tidyverse)
# library(ape)
# library(dplyr)
# library(phytools)
# library(PhD)
# library(iNEXTPD2)
#
# library(magrittr)
# library(chaoUtility)
# library(future.apply)
# library(abind)
#
# library(iNEXT.link)
# # library(sets)
#
# source("R//iNEXTlink.R")
# source("R//subfunction.R")
# source("R//myfunc.R")
#
#
# data(Norfolk)
# data(puerto.rico)
# load(file = "data//beetles_treat_pool.rda")
#
# #DataInfo.link ----
# DataInfo.link(puerto.rico$data, diversity = 'TD', datatype="abundance")
# DataInfo.link(puerto.rico$data, diversity = 'PD', datatype="abundance",
#               row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
#
# #SC.link ----
# sc1 = SC.link(Norfolk, datatype = "abundance", nboot = 10)
# ggSC.link(sc1)
#
# #iNEXT.link ----
# out1 <- iNEXT.link(Norfolk, diversity = 'TD',datatype = "abundance", nboot = 10)
# out2 <- iNEXT.link(puerto.rico$data, diversity = 'PD', datatype="abundance",
#                    row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree, nboot = 3, PDtype = 'PD')
#
# #ggiNEXT.link ----
# ggiNEXT.link(outcome = out1, diversity = 'TD', type = 1)
# ggiNEXT.link(outcome = out1, diversity = 'TD', type = 2)
# ggiNEXT.link(outcome = out1, diversity = 'TD', type = 3)
#
# ggiNEXT.link(outcome = out2, diversity = 'PD', type = 1)
# # ggiNEXT.link(outcome = out2, diversity = 'PD', type = 2)
# ggiNEXT.link(outcome = out2, diversity = 'PD', type = 3)
#
# #Asy.link ----
# asy1 <- Asy.link(Norfolk, diversity = 'TD', datatype = "abundance", nboot = 10)
# asy2 <- Asy.link(puerto.rico$data, diversity = 'PD', datatype = "abundance", nboot = 10,
#                  row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
# #ggAsy.link ----
# ggAsy.link(asy1)
# ggAsy.link(asy2, diversity = 'PD')
#
# #Obs.link ----
# obs1 <- Obs.link(Norfolk, diversity = 'TD', datatype = "abundance", nboot = 10)
# obs2 <- Obs.link(data = puerto.rico$data, diversity = 'PD', datatype = "abundance",nboot = 10,
#                  row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
# #ggObs.link ----
# ggObs.link(obs1)
# ggObs.link(obs2, diversity = 'PD')
#
# #estimateD.link ----
# est1 <- estimateD.link(Norfolk, diversity = 'TD', datatype="abundance", base="coverage", level=0.7, nboot = 10)
# est2 <- estimateD.link(Norfolk, diversity = 'TD', datatype="abundance", base="size", level=0.7, nboot = 10)
# est3 <- estimateD.link(data = puerto.rico$data, diversity = 'PD', datatype="abundance", base="coverage",
#                        level=0.7, nboot = 10,
#                        row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
# est1
# est2
# est3
#
#
# #iNEXTbeta.link ----
# beta1 = iNEXTbeta.link(data = puerto.rico$data%>%lapply(function(x) round(x/10)),
#                        level = seq(0.5, 0.9, 0.4), datatype='abundance',q = c(0, 1, 2),
#                        diversity = 'TD', nboot = 10, conf = 0.95)
# beta2 = iNEXTbeta.link(data = puerto.rico$data%>%lapply(function(x) round(x/10)),
#                        level = seq(0.5, 0.9, 0.4), datatype='abundance',q = c(0, 1, 2),
#                        diversity = 'PD', nboot = 10, conf = 0.95,
#                        row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
# #ggiNEXTbeta.link ----
# ggiNEXTbeta.link(beta1,diversity = 'TD', type = 'B')
# ggiNEXTbeta.link(beta1,diversity = 'TD', type = 'D')
# ggiNEXTbeta.link(beta2,diversity = 'PD', type = 'B')
# ggiNEXTbeta.link(beta2,diversity = 'PD', type = 'D')
#
#
# #Spec.link ----
# Est <- Spec.link(data = Norfolk, datatype = "abundance", q = c(0,1,2),
#                  nboot = 30, method = "Estimated")
# Emp <- Spec.link(data = Norfolk, datatype = "abundance", q = c(0,1,2),
#                  nboot = 30, method = "Empirical")
# ggSpec.link(Est)
# ggSpec.link(Emp)
#
