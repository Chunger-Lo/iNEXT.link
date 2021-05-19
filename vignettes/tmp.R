length(beetles[[1]]%>%rowSums()>0)
length(beetles[[2]]%>%rowSums()>0)

length(beetles[[1]]%>%colSums()>0)
length(beetles[[2]]%>%colSums()>0)

x = beetles

data_long <- lapply(x, function(tab){
  tab = rownames_to_column(tab, "row.name")
  long = gather(data = tab,key = "col.name", value= "abundance", -row.name)
  # long = mutate(long,int_name = paste0(col.name, "x", row.name))
  long = unite(long, "int_name",row.name:col.name)
  long$abundance = as.numeric(long$abundance)
  return(long)
})

names_tab = lapply(seq_along(x), function(i){
  as.set(data_long[[i]]$int_name)
})
# res_set = set()
# for(i in 1:(length(names_tab)-1)){
#   res_set = set_union(names_tab[[i]], names_tab[[i+1]])
# }
res_set = set()
for(i in 1:length(names_tab)){
  res_set = set_union(names_tab[[i]], res_set)
}

combined = data.frame(sp = sapply(res_set, as.character))
for(i in 1:length(names_tab)){
  combined = combined%>%
    left_join( data_long[[i]], by = c("sp" = "int_name"))
}

colnames(combined)[2:(1+length(x))] = paste0("abundance", 1:length(x))
combined[is.na(combined)] = 0
combined = column_to_rownames(combined, "sp")

sum(combined[,1] == combined[,2]) / nrow(combined)
