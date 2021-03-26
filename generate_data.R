Shelfanger <- read.table("F://Chao//Network diversity//shinyapp//heir.entropy//data//Shelfanger.txt")
Hickling <- read.table("F://Chao//Network diversity//shinyapp//heir.entropy//data/Hickling.txt")
Norfolk <- list(Shelfanger = Shelfanger,Hickling = Hickling)
usethis::use_data(Norfolk, overwrite = T)

Caguana <- read.table("F://Chao//Network diversity//shinyapp//heir.entropy//data/Caguana.txt")
Cialitos <- read.table("F://Chao//Network diversity//shinyapp//heir.entropy//data/Cialitos.txt")
Cordillera <- read.table("F://Chao//Network diversity//shinyapp//heir.entropy//data/Cordillera.txt")
Fronton <- read.table("F://Chao//Network diversity//shinyapp//heir.entropy//data/Fronton.txt")
data <- list(Cordillera = Cordillera,Caguana = Caguana,Fronton = Fronton,Cialitos = Cialitos)
puerto.rico <- list(Cordillera = Cordillera,Caguana = Caguana,Fronton = Fronton,Cialitos = Cialitos)

library(phytools)
rowtree <- read.newick("F://Chao//Network diversity//shinyapp//heir.entropy//tree/rowtree.txt")
coltree <- read.newick("F://Chao//Network diversity//shinyapp//heir.entropy//tree/coltree.txt")

puerto.rico = list()
puerto.rico[["data"]] = data
puerto.rico[["col.tree"]] = coltree
puerto.rico[["row.tree"]] = rowtree
usethis::use_data(puerto.rico, overwrite = T)



