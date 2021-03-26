library(plyr)
library(sets)
library(tidyverse)
library(dplyr)
library(tibble)
source("iNEXT3D_0316.r")
source("iNEXT_beta.r")
source("iNEXTlink.r")
source("myfunc.r")
source("subfunction.r")

x = Norfolk
data_type = "abundance"
q= c(0,1,2);level = "taxonomic"
nboot = 30; conf=0.95;by = "coverage"
coverage_expected = seq(0.5,1,0.05);max_alpha_coverage= F


