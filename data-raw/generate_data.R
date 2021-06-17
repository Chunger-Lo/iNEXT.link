library(phytools)
library(dplyr)
library(stringr)
library(tidyverse)
# Norfolk ----
Shelfanger <- read.table("data-raw//Norfolk//Shelfanger.txt")
Hickling <- read.table("data-raw//Norfolk//Hickling.txt")
Norfolk <- list(Shelfanger = Shelfanger,Hickling = Hickling)
usethis::use_data(Norfolk, overwrite = T)

# puerto.rico ----
Caguana <- read.table("data-raw//puerto.rico//Caguana.txt")
Cialitos <- read.table("data-raw//puerto.rico//Cialitos.txt")
Cordillera <- read.table("data-raw//puerto.rico//Cordillera.txt")
Fronton <- read.table("data-raw//puerto.rico//Fronton.txt")
data <- list(Cordillera = Cordillera,Caguana = Caguana,Fronton = Fronton,Cialitos = Cialitos)
puerto.rico <- list(Cordillera = Cordillera,Caguana = Caguana,Fronton = Fronton,Cialitos = Cialitos)


rowtree <- read.newick("data-raw//puerto.rico//rowtree.txt")
coltree <- read.newick("data-raw//puerto.rico//coltree.txt")

puerto.rico = list()
puerto.rico[["data"]] = data
puerto.rico[["col.tree"]] = coltree
puerto.rico[["row.tree"]] = rowtree
usethis::use_data(puerto.rico, overwrite = T)

# Beetles interaction ----
# load(file = "data-raw//beetles//beetles_tree.rda")
load(file = "data-raw//beetles//for_anne.rda")
tab_selected = tab[,-c(1,2,3,8,9)]
#Tilia_cordata  $Ernoporus_tiliae

### temp remove
tab_selected[tab_selected$tree_species=="Tilia_cordata","Ernoporus.tiliae"] = 0
###
colnames(tab_selected) = gsub('\\.', '_',colnames(tab_selected))
tab_selected$tree_species= gsub('\\.', '_',tab_selected$tree_species)

# speices_name = tab_selected%>%colnames()%>%
#   .[5:length(.)]
#
# ind = sapply(speices_name, function(name){
#   name%in%tr$tip.label
# })

names(tab_selected) <- gsub("Anobium_fulvicorne", "Hemicoelus_fulvicornis", names(tab_selected))
names(tab_selected) <- gsub("Xyleborus_saxesenii", "Xyleborinus_saxesenii", names(tab_selected))
names(tab_selected) <- gsub("Xyleborus_dispar", "Anisandrus_dispar", names(tab_selected))
names(tab_selected) <- gsub("Xyleborus_germanus", "Xylosandrus_germanus", names(tab_selected))
names(tab_selected) <- gsub("Phymatodes_alni", "Poecilium_alni", names(tab_selected))
names(tab_selected) <- gsub("Phymatodes_pusillus", "Poecilium_pusillum", names(tab_selected))
names(tab_selected) <- gsub("Xyloterus_signatus", "Trypodendron_signatum", names(tab_selected))
names(tab_selected) <- gsub("Rhopalopus_femoratus", "Ropalopus_femoratus", names(tab_selected))
names(tab_selected) <- gsub("Molorchus_umbellatarum", "Glaphyra_umbellatarum", names(tab_selected))
names(tab_selected) <- gsub("Nemosoma_elongatum", "Nemozoma_elongatum", names(tab_selected))
names(tab_selected) <- gsub("Grammoptera_variegata", "Grammoptera_abdominalis", names(tab_selected))
names(tab_selected) <- gsub("Dasytes_cyaneus", "Dasytes_caeruleus", names(tab_selected))
names(tab_selected) <- gsub("Cryphalus_abietis", "Cryphalus_asperatus", names(tab_selected))
names(tab_selected) <- gsub("Phyllodrepa_ioptera", "Dropephylla_ioptera", names(tab_selected))
names(tab_selected) <- gsub("Rhaphitropis_machicus", "Rhaphitropis_marchica", names(tab_selected))
treatG = tab_selected%>%filter(treatment == 'G')%>%.[,-c(2,3)]
treatO = tab_selected%>%filter(treatment == 'O')%>%.[,-c(2,3)]



net_treat_plot = list()
net_treat_plot[['G']] = list()
net_treat_plot[['G']][['A']] = treatG%>%filter(plot == 'A')%>%.[,-1]%>%column_to_rownames('tree_species')
net_treat_plot[['G']][['B']] = treatG%>%filter(plot == 'B')%>%.[,-1]%>%column_to_rownames('tree_species')
net_treat_plot[['G']][['C']] = treatG%>%filter(plot == 'C')%>%.[,-1]%>%column_to_rownames('tree_species')
net_treat_plot[['O']] = list()
net_treat_plot[['O']][['A']] = treatO%>%filter(plot == 'A')%>%.[,-1]%>%column_to_rownames('tree_species')
net_treat_plot[['O']][['B']] = treatO%>%filter(plot == 'B')%>%.[,-1]%>%column_to_rownames('tree_species')
net_treat_plot[['O']][['C']] = treatO%>%filter(plot == 'C')%>%.[,-1]%>%column_to_rownames('tree_species')

beetles_treat_plot = list()
beetles_treat_plot[["data"]] = net_treat_plot
beetles_treat_plot[["col.tree"]] = tr
usethis::use_data(beetles_treat_plot, overwrite = T)

beetles_treat_plot_remove_dominance = list()
beetles_treat_plot_remove_dominance[["data"]] = net_treat_plot
beetles_treat_plot_remove_dominance[["col.tree"]] = tr
# usethis::use_data(beetles_treat_plot, overwrite = T)
usethis::use_data(beetles_treat_plot_remove_dominance, overwrite = T)
### treat pool
net_treat_pool = list()
net_treat_pool[['G']] = treatG%>%.[,-1]%>%group_by(tree_species)%>%summarise_all(sum)%>%column_to_rownames('tree_species')
net_treat_pool[['O']] = treatO%>%.[,-1]%>%group_by(tree_species)%>%summarise_all(sum)%>%column_to_rownames('tree_species')

beetles_treat_pool = list()
beetles_treat_pool[["data"]] = net_treat_pool
beetles_treat_pool[["col.tree"]] = tr
usethis::use_data(beetles_treat_pool, overwrite = T)

beetles_treat_pool_remove = list()
beetles_treat_pool_remove[["data"]] = net_treat_pool
beetles_treat_pool_remove[["col.tree"]] = tr
usethis::use_data(beetles_treat_pool_remove, overwrite = T)



## by plot
plotA = tab_selected%>%filter(plot == 'A')%>%.[,-c(1,3)]
plotB = tab_selected%>%filter(plot == 'B')%>%.[,-c(1,3)]
plotC = tab_selected%>%filter(plot == 'C')%>%.[,-c(1,3)]

interactions = list()
interactions[['A']] = list()
interactions[['A']][['G']] = plotA%>%filter(treatment == 'G')%>%.[,-1]%>%column_to_rownames('tree_species')
interactions[['A']][['O']] = plotA%>%filter(treatment == 'O')%>%.[,-1]%>%column_to_rownames('tree_species')

interactions[['B']] = list()
interactions[['B']][['G']] = plotB%>%filter(treatment == 'G')%>%.[,-1]%>%column_to_rownames('tree_species')
interactions[['B']][['O']] = plotB%>%filter(treatment == 'O')%>%.[,-1]%>%column_to_rownames('tree_species')

interactions[['C']] = list()
interactions[['C']][['G']] = plotC%>%filter(treatment == 'G')%>%.[,-1]%>%column_to_rownames('tree_species')
interactions[['C']][['O']] = plotC%>%filter(treatment == 'O')%>%.[,-1]%>%column_to_rownames('tree_species')

beetles_plot_treat = list()
beetles_plot_treat[["data"]] = interactions
beetles_plot_treat[["col.tree"]] = tr

beetles_plot_treat_remove_dominance = list()
beetles_plot_treat_remove_dominance[["data"]] = interactions
beetles_plot_treat_remove_dominance[["col.tree"]] = tr
usethis::use_data(beetles_plot_treat, overwrite = T)
usethis::use_data(beetles_plot_treat_remove_dominance, overwrite = T)

### treat pool
net_plot_pool = list()
net_plot_pool[['A']] = plotA%>%.[,-1]%>%group_by(tree_species)%>%summarise_all(sum)%>%column_to_rownames('tree_species')
net_plot_pool[['B']] = plotB%>%.[,-1]%>%group_by(tree_species)%>%summarise_all(sum)%>%column_to_rownames('tree_species')
net_plot_pool[['C']] = plotC%>%.[,-1]%>%group_by(tree_species)%>%summarise_all(sum)%>%column_to_rownames('tree_species')

beetles_plot_pool = list()
beetles_plot_pool[["data"]] = net_plot_pool
beetles_plot_pool[["col.tree"]] = tr
usethis::use_data(beetles_plot_pool, overwrite = T)

beetles_plot_pool_remove = list()
beetles_plot_pool_remove[["data"]] = net_plot_pool
beetles_plot_pool_remove[["col.tree"]] = tr
usethis::use_data(beetles_plot_pool_remove, overwrite = T)



#
# # beetles(pool) ----
# treatG = tab_selected%>%filter(treatment == 'G')%>%.[,-c(1,2,3)]
# treatO = tab_selected%>%filter(treatment == 'O')%>%.[,-c(1,2,3)]
#
# treatG_pool = treatG%>%group_by(tree_species)%>%summarize_all(sum)%>%column_to_rownames('tree_species')
# treatO_pool = treatO%>%group_by(tree_species)%>%summarize_all(sum)%>%column_to_rownames('tree_species')
# interactions = list()
# interactions[['G']] = treatG_pool
# interactions[['O']] = treatO_pool
#
# beetles_pool = list()
# beetles_pool[["data"]] = interactions
# beetles_pool[["col.tree"]] = tr
# usethis::use_data(beetles_pool, overwrite = T)



