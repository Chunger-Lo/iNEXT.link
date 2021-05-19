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
load(file = "data-raw//beetles//beetles_tree.rda")
tab_selected = tab[,-c(1,2,3,8,9)]
colnames(tab_selected) = colnames(tab_selected)%>%str_replace_all('\\.','_')

speices_name = tab_selected%>%colnames()%>%
  .[5:length(.)]

ind = sapply(speices_name, function(name){
  name%in%tr$tip.label
})
speices_name[!ind]

net_treat_plot = list()

treatG = tab_selected%>%filter(treatment == 'G')%>%.[,-c(2,3)]
treatO = tab_selected%>%filter(treatment == 'O')%>%.[,-c(2,3)]

net_treat_plot[['G']] = list()
net_treat_plot[['G']][['A']] = treatG%>%filter(plot == 'A')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
net_treat_plot[['G']][['B']] = treatG%>%filter(plot == 'B')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
net_treat_plot[['G']][['C']] = treatG%>%filter(plot == 'C')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
net_treat_plot[['O']] = list()
net_treat_plot[['O']][['A']] = treatO%>%filter(plot == 'A')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
net_treat_plot[['O']][['B']] = treatO%>%filter(plot == 'B')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
net_treat_plot[['O']][['C']] = treatO%>%filter(plot == 'C')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]

beetles_treat_plot = list()
beetles_treat_plot[["data"]] = net_treat_plot
beetles_treat_plot[["col.tree"]] = tr
usethis::use_data(beetles_treat_plot, overwrite = T)


## by plot
plotA = tab_selected%>%filter(plot == 'A')%>%.[,-c(1,3)]
plotB = tab_selected%>%filter(plot == 'B')%>%.[,-c(1,3)]
plotC = tab_selected%>%filter(plot == 'C')%>%.[,-c(1,3)]

interactions = list()
interactions[['A']] = list()
interactions[['A']][['G']] = plotA%>%filter(treatment == 'G')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
interactions[['A']][['O']] = plotA%>%filter(treatment == 'O')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]

interactions[['B']] = list()
interactions[['B']][['G']] = plotB%>%filter(treatment == 'G')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
interactions[['B']][['O']] = plotB%>%filter(treatment == 'O')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]

interactions[['C']] = list()
interactions[['C']][['G']] = plotC%>%filter(treatment == 'G')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
interactions[['C']][['O']] = plotC%>%filter(treatment == 'O')%>%.[,-1]%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]

beetles_plot_treat = list()
beetles_plot_treat[["data"]] = interactions
beetles_plot_treat[["col.tree"]] = tr
usethis::use_data(beetles_plot_treat, overwrite = T)

# beetles(pool) ----
treatG = tab_selected%>%filter(treatment == 'G')%>%.[,-c(1,2,3)]
treatO = tab_selected%>%filter(treatment == 'O')%>%.[,-c(1,2,3)]

treatG_pool = treatG%>%group_by(tree_species)%>%summarize_all(sum)%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
treatO_pool = treatO%>%group_by(tree_species)%>%summarize_all(sum)%>%column_to_rownames('tree_species')%>%.[,speices_name[ind]]
interactions = list()
interactions[['G']] = treatG_pool
interactions[['O']] = treatO_pool

beetles_pool = list()
beetles_pool[["data"]] = interactions
beetles_pool[["col.tree"]] = tr
usethis::use_data(beetles_pool, overwrite = T)


