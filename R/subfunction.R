#########
datainfphy <- function(data, row.tree = NULL,col.tree = NULL, datatype){
  if (class(data)!="dataframe") data <- as.data.frame(data)
  # if(datatype == "abundance"){
  #   rownames(res) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "Coverage","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10")
  #   rownames(tmp) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")
  #
  # }
  # if(datatype == "incidence_freq"){
  #   rownames(tmp) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")
  #   rownames(res) <- c("U", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "Coverage","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")
  # }
  res <- matrix(0,12,1,dimnames=list(1:12, "value"))
  rownames(res) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")

  res[1,1] <- as.integer(sum(data))
  res[2,1] <-  sum(ncol(data),nrow(data))
  res[3,1] <-  sum(data>0)
  res[4,1] <-  nrow(data)
  res[5,1] <-  ncol(data)
  res[6,1] <-  round(sum(data>0)/ncol(data)/nrow(data),4)
  res[7,1] <-  sum(data == 1)
  res[8,1] <-  sum(data == 2)
  phy <- create.aili(data,row.tree = row.tree,col.tree = col.tree)
  res[9,1] <- sum(phy$branch.length[phy$branch.abun==1])
  res[10,1] <- sum(phy$branch.length[phy$branch.abun==2])
  res[11,1] <- nrow(phy)
  res[12,1] <- sum(phy$branch.length*phy$branch.abun)/sum(data)

  res = res%>%t()%>%as.data.frame()
  return(res)
}
datainf <- function(data, datatype){
  if (class(data)!="dataframe") data <- as.data.frame(data)
  res <- matrix(0,17,1,dimnames=list(1:17, "value"))

  if(datatype == "abundance"){
    rownames(res) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "Coverage","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10")
  }
  if(datatype == "incidence_freq"){
    rownames(res) <- c("U", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "Coverage","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")
  }

  res[1,1] <- as.integer(sum(data))
  res[2,1] <-  sum(ncol(data),nrow(data))
  res[3,1] <-  sum(data>0)
  res[4,1] <-  nrow(data)
  res[5,1] <-  ncol(data)
  res[6,1] <-  round(sum(data>0)/ncol(data)/nrow(data),4)
  res[8:17,1] <- c(sum(data==1),sum(data==2),sum(data==3),sum(data==4),sum(data==5),sum(data==6),sum(data==7),sum(data==8),sum(data==9),sum(data==10))
  f1 = sum(data==1)
  f2 = sum(data==2)
  n = sum(data)
  res[7,1] <- round(1 - (f1/n)*((n-1)*f1/((n-1)*f1+2*f2)),4) #C.hat

  res = res%>%t()%>%as.data.frame()
  return(res)
}
#plot structure
plot.tree2 <- function(mat){
  #number of lower level must be large than or equal to the number of higher level
  t <- apply(mat,MARGIN = 1, function(x) length(unique(x)))
  if(sum((t[-1]-t[-length(t)])<0)>0) stop("number of lower level must be large than or equal to the number of higher level, please renew your structure matrix.")
  rownames(mat) <- paste0("level", 1:nrow(mat))
  colnames(mat) <- paste0("community", 1:ncol(mat))
  mat <- data.frame(t(mat), stringsAsFactors = F)
  m <- ncol(mat)
  mat$pathString <- apply(mat,1,paste,collapse="/")
  population <- as.Node(mat)
  useRtreeList <- ToListExplicit(population, unname = TRUE)

  # radialNetwork(useRtreeList, fontSize = 10, opacity = 0.9)
  diagonalNetwork(useRtreeList,fontSize = 27, opacity = 10, linkColour = "#828282", nodeStroke = "#6495ED")
}

Diversity_profile <- function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 1:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum((1-p1)^r/r)))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      #ifelse(A==0,NA,A^(1/(1-q)))
      A^(1/(1-q))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      r <- 0:(n-1)
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}
Diversity_profile_MLE <- function(x,q){
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}
Diversity_Tsallis <- function(x,q){
  qD = Diversity_profile(x, q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}
Diversity_Tsallis_MLE <- function(x,q){
  qD = Diversity_profile_MLE(x,q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}
bootstrap_forq = function(data,B,q,conf,FUNNAME){
  data <- data[data!=0]
  n <- sum(data)
  f1 = sum(data==1); f2 = sum(data==2)
  f0 = ceiling(ifelse( f2>0, (n-1)*f1^2/n/2/f2, (n-1)*f1*(f1-1)/2/n ))
  C_hat = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
  lamda_hat = (1-C_hat)/sum((data/n)*(1-data/n)^n)
  pi_hat = (data/n)*(1-lamda_hat*((1-data/n)^n))
  p_hat = c( pi_hat, rep( (1-C_hat)/f0, f0 ))
  random = rmultinom( B, n, p_hat )
  #Bt_estimate <- sapply(c(1:B),function(i) FUNNAME(random[,i],q))
  Bt_estimate <- apply(random,MARGIN = 2,function(i) FUNNAME(i,q))
  estimate <- FUNNAME(data,q)
  #Interval_mean = apply( Bt_estimate, 1, mean)
  Interval_mean = rowMeans(Bt_estimate)
  Interval_sd = apply(Bt_estimate, 1, sd)
  Interval_quantileL = apply( Bt_estimate, 1, quantile, p=(1-conf)/2)
  Interval_quantileU = apply( Bt_estimate, 1, quantile, p=1-(1-conf)/2)
  Upper_bound = estimate+Interval_quantileU-Interval_mean
  Lower_bound = estimate+Interval_quantileL-Interval_mean
  result <- cbind("estimate"=estimate,"sd"=Interval_sd,"LCL"=Lower_bound,"UCL"=Upper_bound)
  result
}
# data = tmp;
# B = 10;
# q=c(0,1,2);conf= 0.95
MakeTable_Proposeprofile = function(data, B, q, conf){
  Diversity = bootstrap_forq(data, B, q, conf, Diversity_profile)
  Entropy = bootstrap_forq(data, B, q, conf, Diversity_Tsallis)
  # tmp <- Diversity_Tsallis(Diversity[,1],q)
  # Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
  output = rbind(data.frame("Order.q" = q,"Target"="Diversity","Estimate"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4]),
                 data.frame("Order.q" = q,"Target"="Entropy","Estimate"=Entropy[,1],"s.e."=Entropy[,2],"LCL"=Entropy[,3],"UCL"=Entropy[,4]))
  output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)

  return(output)
}
MakeTable_Empericalprofile = function(data, B, q, conf){
  Diversity = bootstrap_forq( data, B,q,conf,Diversity_profile_MLE)
  Entropy = bootstrap_forq( data, B,q,conf,Diversity_Tsallis_MLE)
  # tmp <- Diversity_Tsallis(Diversity[,1],q)
  # Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
  output = rbind(data.frame("Order.q" = q,"Target"="Diversity","Emperical"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4]),
                 data.frame("Order.q" = q,"Target"="Entropy","Emperical"=Entropy[,1],"s.e."=Entropy[,2],"LCL"=Entropy[,3],"UCL"=Entropy[,4]))
  output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)

  return(output)
}


### Evenness
# x = Norfolk
# x = data
# x = puerto.rico$data

Evenness.profile <- function(x, q, datatype = c("abundance","incidence_freq"), method, E.class, C = NULL) {
  # data_long <- lapply(x, function(tab){
  #   as.matrix(tab)%>%c()}
  #   ## 拉長
  # )

  if (method == "Estimated") {
    estqD = estimate3D(x, class = 'TD', q, datatype, base = "coverage", level = C, nboot = 0)
    estS = estimate3D(x, class = 'TD', 0, datatype, base = "coverage", level = C, nboot = 0)

    estqD = estimate3D(x, class = 'TD', q, datatype, base = "coverage", level = C, nboot = 0)
    estS = estimate3D(x, class = 'TD', 0, datatype, base = "coverage", level = C, nboot = 0)

    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q, empqD[empqD$Assemblage == names(x)[k], "qD"], empS[empS$Assemblage == names(x)[k], "qD"], i, x[[k]]/sum(x[[k]])))
      if (class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow = 1))}
      rownames(tmp) = q
      tmp
    })
  } else if (method == "Empirical") {

    empqD = ObsND(x, q = q, datatype = datatype, nboot = 30)
    empS = empqD%>%filter(Order.q == 0)

    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q,
                                                       qD =empqD[empqD$Assemblage == names(x)[k], "qD"],
                                                       S = empS[empS$Assemblage == names(x)[k], "qD"],
                                                       E.class = i))
                                                       # x[[k]]/sum(x[[k]])))
      if (class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow = 1))}
      rownames(tmp) = q
      tmp
    })
    # empqD = Obs3D(x, class = 'TD', q = q, datatype = datatype, nboot = 0)
    # empS = Obs3D(x, class = 'TD', q = 0, datatype = datatype, nboot = 0)
    #
    # out = lapply(E.class, function(i) {
    #   tmp = sapply(1:length(x), function(k) even.class(q, empqD[empqD$Assemblage == names(x)[k], "qD"], empS[empS$Assemblage == names(x)[k], "qD"], i, x[[k]]/sum(x[[k]])))
    #   if (class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow = 1))}
    #   rownames(tmp) = q
    #   tmp
    # })
  }

  names(out) = paste("E", E.class, sep="")
  return(out)
}

###
ready4beta <- function(x){
  ### 拉長所有interaction, 視為不同subplot, 取聯集補0
  ## transform 2d matrix to vector
  ## expand to the union of all networks
  ## replace na to zero
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
}
#########
long_to_wide = function(data_long = data_gamma){
  n <- 2
  reg <- paste0('((?:[^_]*_){',n,'})(.*)')

  temp = data_long%>%as.data.frame()%>%rownames_to_column("sp")

  splitted  = sapply(temp$sp, function(str){
    gsubfn::strapply(str, reg, c, simplify =~ sub('_$', '', x))
  })%>%t()%>%as.data.frame()%>% remove_rownames()%>%set_colnames(c('row_sp', 'col_sp'))

  mat = splitted%>%cbind(abundance = temp$.)%>%
    tidyr::spread(key = "col_sp", value = "abundance")%>%column_to_rownames("row_sp")
  mat[is.na(mat)] = 0

  return(mat)
}

create.aili <- function(data,row.tree = NULL,col.tree = NULL) {
  if (class(data)[1] != "matrix") { data <- as.matrix(data) }

  if ((is.null(row.tree)) == 0 & (is.null(col.tree) == 1)){
    tip <- row.tree$tip.label[-match(rownames(data),row.tree$tip.label)]
    mytree <- drop.tip(row.tree,tip)
    mytree <- phylo2phytree(mytree)
    tmp <- apply(data, 2, function(abun){
      phyExpandData(x=abun, labels=rownames(data), phy=mytree, datatype="abundance")
    })
    tmp <- lapply(1:length(tmp), function(x){
      tmp1 <- as.data.frame(tmp[[x]])
      tmp1$spe.c <- colnames(data)[x]
      tmp1$interaction <- paste(tmp1$label,tmp1$spe.c,sep = "_")
      tmp1
    })
    tmp <- do.call(rbind,tmp)
    out <- data.frame(branch.abun = tmp$branch.abun, branch.length = tmp$branch.length,
                      tgroup = tmp$tgroup,interaction = tmp$interaction)
  }
  ## beetles (only coltree)
  if ((is.null(row.tree) == 1) & (is.null(col.tree) == 0)){
    ind = is.na(match(colnames(data),col.tree$tip.label))
    tip_notin_data <- col.tree$tip.label[ind]
    mytree <- drop.tip(col.tree,tip_notin_data)
    mytree <- phylo2phytree(mytree)

    tmp <- apply(data, 1, function(abun){
      # phyBranchAL_Abu(phylo = mytree, data = abun, rootExtend = T, refT = NULL)
      phyExpandData(x=abun, labels=colnames(data), phy=mytree, datatype="abundance")
    })
    tmp <- lapply(1:length(tmp), function(x){
      tmp1 <- tmp[[x]]%>%as.data.frame()
      tmp1$spe.r <- rownames(data)[x]
      tmp1$interaction <- paste(tmp1$spe.r,tmp1$label,sep = "_")
      tmp1
    })
    tmp <- do.call(rbind,tmp)
    out <- data.frame(branch.abun = tmp$branch.abun, branch.length = tmp$branch.length,tgroup = tmp$tgroup,interaction = tmp$interaction)
    out$branch.abun[is.na(out$branch.abun)] = 0
  }

  if ((is.null(row.tree) == 0) & (is.null(col.tree) == 0)){
    col.tip <- col.tree$tip.label[-match(colnames(data),col.tree$tip.label)]
    mytree.col <- drop.tip(col.tree,col.tip)
    mytree.col <- phylo2phytree(mytree.col)
    row.tip <- row.tree$tip.label[-match(rownames(data),row.tree$tip.label)]
    mytree.row <- drop.tip(row.tree,row.tip)
    mytree.row <- phylo2phytree(mytree.row)

    tmp <- apply(data, 1, function(abun){
      phyExpandData(x=abun, labels=colnames(data), phy=mytree.col, datatype="abundance")
    })
    for (i in 1:length(tmp)) {
      tmp[[i]]$Species <- names(tmp)[i]
    }
    tmp <- do.call(rbind,tmp)
    #tmp[is.na(tmp)] <- 0
    t1 <- tmp%>%as.data.frame()%>%dplyr::select(label, branch.abun, Species)
    # t1 <- as.data.frame(tmp[,c(4,7,8)])
    t1 <- spread(t1,Species,branch.abun)%>%column_to_rownames("label")%>%as.matrix()
    # rownames(t1) <- t1[,1]
    mat <- as.matrix(t1)
    label <- unique(tmp$label)
    species <- unique(tmp$Species)

    # mat <- matrix(0,length(label),length(species))
    # rownames(mat) <- label
    # colnames(mat) <- species
    # for (i in seq_len(nrow(mat))) {
    #   for (j in seq_len(ncol(mat))) {
    #     if (sum(tmp$label == label[i] & tmp$Species == species[j])>0) {
    #       mat[i,j] <- tmp[tmp$label == label[i] & tmp$Species == species[j],]$branch.abun
    #     }
    #   }
    # }
    ### col by col
    tmp1 <- apply(mat, 1, function(abun){
      phyExpandData(x=abun, labels=rownames(data), phy=mytree.row, datatype="abundance")
    })
    for (i in 1:length(tmp1)) {
      tmp1[[i]]$Species <- names(tmp1)[i]
    }
    tmp1 <- do.call(rbind,tmp1)
    # t2 <- as.data.frame(tmp1[,c(3,4,5,7,8)])
    t2 <- as.data.frame(tmp1[, c("branch.length", 'label','tgroup','node.age','branch.abun','Species')])

    tmp1[is.na(tmp1)] <- 0
    t2$r.length <- 0
    t2$r.group <- 0
    t2$group <- 0
    for (i in seq_len(length(label))) {
      t2[t2$Species == label[i],]$r.length <-  tmp[tmp$label == label[i],]$branch.length[1]
      t2[t2$Species == label[i],]$r.group <-  tmp[tmp$label == label[i],]$tgroup[1]
    }
    t2$group <- ifelse(t2$tgroup == t2$r.group,t2$tgroup,"Inode")
    t2$interaction <- paste(t2$label,t2$Species,sep = "_")

    out <- data.frame(branch.abun = t2$branch.abun, branch.length = t2$branch.length*t2$r.length,
                      tgroup = t2$group, interaction = t2$interaction)
  }
  # out <- out[out$branch.abun > 0,]
  # out <- na.omit(out)
  # out <- out[order(out$tgroup,decreasing = T),]
  rownames(out) <- NULL
  return(out)
}
# datainfphy <- function(data, row.tree = NULL,col.tree = NULL){
#   if (class(data)!="dataframe") data <- as.data.frame(data)
#   tmp <- matrix(0,12,1,dimnames=list(1:12, "value"))
#   rownames(tmp) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")
#   tmp[1,1] <- as.integer(sum(data))
#   tmp[2,1] <-  sum(ncol(data),nrow(data))
#   tmp[3,1] <-  sum(data>0)
#   tmp[4,1] <-  nrow(data)
#   tmp[5,1] <-  ncol(data)
#   tmp[6,1] <-  round(sum(data>0)/ncol(data)/nrow(data),4)
#   tmp[7,1] <-  sum(data == 1)
#   tmp[8,1] <-  sum(data == 2)
#   phy <- create.aili(data,row.tree = row.tree,col.tree = col.tree)
#   tmp[9,1] <- sum(phy$branch.length[phy$branch.abun==1])
#   tmp[10,1] <- sum(phy$branch.length[phy$branch.abun==2])
#   tmp[11,1] <- nrow(phy)
#   tmp[12,1] <- sum(phy$branch.length*phy$branch.abun)/sum(data)
#   out <- data.frame(tmp)
#   return(out)
# }
expanddata <- function(data){
  out <- lapply(data,function(x){
    if(ncol(x) != 1){
      tmp <- as.data.frame(x)
    tmp$spe.r <- rownames(tmp)
    tmp1 <- gather(tmp,spe.c,abun,-spe.r)
    tmp2 <- data.frame(link = paste(tmp1$spe.r,tmp1$spe.c,sep = "_"),abun = tmp1$abun)
    out <- tmp2[tmp2$abun >0,]
    }
    if(ncol(x) == 1){
      tmp <- as.data.frame(x)
    # tmp$spe.r <- rownames(tmp)
    # tmp1 <- gather(tmp,spe.c,abun,-spe.r)
    tmp <- data.frame(link = rownames(x),abun = tmp$`x[[1]]`)
    # tmp2[tmp2$abun >0,]
    out <- tmp[tmp$abun >0,]
    }
    out
  })
  tmp <- do.call(rbind,out)
  link <- unique(tmp$link)
  dama <- matrix(0,length(link),length(data),dimnames = list(link,names(data)))
  for (i in 1:length(data)) {
    dama[match(out[[i]]$link,link),i] <- out[[i]]$abun
  }
  as.data.frame(dama)
}


sample.boot <- function(exdata,B) {
  gamma <- rowSums(exdata)
  n <- sum(sapply(unique(gamma), function(x){x*sum(gamma == x)}))
  f1 <- sum(gamma == 1)
  f2 <- sum(gamma == 2)
  f0.hat <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
  adjust.pi <- function(data) {
    f1 <- sum(data == 1)
    f2 <- sum(data == 2)
    n <- sum(sapply(unique(data), function(x){x*sum(data == x)}))
    C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
    lambda <- ifelse(C.hat != 1,(1-C.hat)/sum(data/n*(1-data/n)^n),0)
    f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
    ###
    f0 <- ifelse( sum(data>0) + f0 > f0.hat +nrow(exdata),nrow(exdata) + f0.hat - sum(data>0) ,f0)
    ###
    p.seen <- data/n*(1-lambda*(1-data/n)^n)
    p.seen <- p.seen[p.seen>0]
    p.unseen <- rep((1-C.hat)/f0,f0)
    list(seen = p.seen,unseen = p.unseen)
  }
  pi <- lapply(1:ncol(exdata), function(x){adjust.pi(exdata[,x])})
  xi <- function(data){
    out <- lapply(1:ncol(data),function(x) {
      tmp <- data[,x]
      ni <- sum(sapply(unique(tmp), function(y){y*sum(tmp == y)}))
      tmp[tmp>0] <- pi[[x]]$seen
      tmp[(length(tmp)+1):(length(tmp)+f0.hat)] <- 0
      tmp[sample(which(tmp == 0),length(pi[[x]]$unseen))] <- pi[[x]]$unseen
      rmultinom(1,ni,tmp)
    })
    as.data.frame(do.call(cbind,out))
  }
  lapply(1:B, function(x){xi(exdata)})
}

hiera <- function(mat){
  if(class(mat) != "matrix") mat <- as.matrix(mat)
  nsite <- ncol(mat)
  H <- nrow(mat)
  newmat <- matrix(0, H, nsite)
  newmat[1, ] <- rep(1, nsite)
  newmat[H, ] <- seq_len(nsite)
  for(i in 2:H){
    newmat[i, ] <- match(mat[i, ], unique(mat[i, ]))
    up <- newmat[i-1, ]; down <- newmat[i, ]
    if(length(unique(down)) != length(unique) && length(unique(up))!=1) {
      for(j in unique(down)){
        id1 <- which(down == j)[1]
        id <- which(down == j)[-1]
        idnew <- id[up[id] != up[id1]]
        new <- max(down)+seq_along(idnew)
        down[idnew] <- new[match(up[idnew], unique(up[idnew]))]
      }
      newmat[i, ] <- down
    }
  }
  newmat
}
hierachical <- function(exdata,mat) {
  hi <- hiera(mat)
  out <- lapply(1:nrow(hi), function(x){
    sapply(unique(hi[x,]), function(y) {
      rowSums(exdata[which(hi[x,]==y)])
    })
  })
  names(out) <- c("Gamma",paste("Alpha",(nrow(hi)-1):1))
  lapply(out, function(x){
    c(unlist(x))
  })
}


sample.boot.PD <- function(data,phydata,B,row.tree = NULL,col.tree = NULL) {
  gamma <- rowSums(expanddata(data))
  n <- sum(sapply(unique(gamma), function(x){x*sum(gamma == x)}))
  f1 <- sum(gamma == 1)
  f2 <- sum(gamma == 2)
  f0.hat <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
  adjust.pi <- function(data) {
    f1 <- sum(data == 1)
    f2 <- sum(data == 2)
    n <- sum(sapply(unique(data), function(x){x*sum(data == x)}))
    C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
    lambda <- ifelse(C.hat != 1,(1-C.hat)/sum(data/n*(1-data/n)^n),0)
    f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
    p.seen <- data/n*(1-lambda*(1-data/n)^n)
    p.seen <- p.seen[p.seen>0]
    p.unseen <- (1-C.hat)/f0
    list(seen = p.seen,unseen = p.unseen)
  }
  pi <- lapply(data, function(x){adjust.pi(x)})
  xi <- function(data,phydata,row.tree,col.tree){
    out <- lapply(1:length(data),function(x) {
      tmp <- data[[x]]
      tmp <- as.matrix(tmp)
      f1 <- sum(tmp == 1)
      f2 <- sum(tmp == 2)
      f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
      phy <- phydata[[x]]
      g1 <- sum(phy$branch.length[phy$branch.abun==1])
      g2 <- sum(phy$branch.length[phy$branch.abun==2])
      g0 <- ifelse(g2>g1*f2/2/f1,(n-1)/n*g1^2/2/g2,(n-1)/n*g1*(f1-1)/2/(f2+1))/f0
      ni <- sum(sapply(unique(c(tmp)), function(y){y*sum(tmp == y)}))
      tmp[tmp>0] <- pi[[x]]$seen
      #tmp[(length(tmp)+1):(length(tmp)+f0.hat)] <- 0
      la <- sample(c(which(tmp == 0),(length(tmp)+1):(length(tmp)+f0.hat)),f0)
      tmp[la[la <= length(tmp)]] <- pi[[x]]$unseen
      p <- c(tmp,rep(pi[[x]]$unseen,sum(la> length(tmp))))
      ai <- rmultinom(1,ni,p)
      tmp <- matrix(ai[1:length(tmp)],nrow(tmp),ncol(tmp),dimnames = list(rownames(tmp),colnames(tmp)))
      out <- create.aili(tmp,row.tree = row.tree,col.tree = col.tree)
      out <- rbind(out,data.frame(branch.abun = ifelse(length(ai)>length(tmp),ai[-c(1:length(tmp))],0),branch.length=g0, tgroup="Tip",interaction = paste("Unseen species",la[la >length(tmp)]-length(tmp))))
      out
    })
    names(out) <- names(data)
    out
  }
  lapply(1:B, function(x){xi(data,phydata,row.tree,col.tree)})
}


hierachical.PD <- function(phydata,mat) {
  hi <- hiera(mat)
  out <- lapply(1:nrow(hi), function(x){
    lapply(unique(hi[x,]), function(y) {
      sub <- which(hi[x,]==y)
      sub <- lapply(sub,function(x){
        phydata[[x]]
      })
      sub <- do.call(rbind,sub)
      sub %>% group_by(interaction)  %>% mutate(branch.length = sum(branch.abun*branch.length)/sum(branch.abun)) %>%
        mutate(branch.abun = sum(branch.abun)) %>% ungroup %>% distinct(interaction,.keep_all = T) %>% as.data.frame
    })
  })
  names(out) <- c("Gamma",paste("Alpha",seq_len(nrow(hi)-1)))
  lapply(1:nrow(hi), function(x){
    tmp <- do.call(rbind,out[[x]])
    tmp <- tmp[tmp$branch.abun >0,]
  })
}

sample.boot.phy <- function(data,B,row.tree = NULL,col.tree = NULL) {
  data <- as.matrix(data)
  # n <- sum(sapply(unique(c(data)), function(x){x*sum(data == x)}))
  n <- sum(data)
  f1 <- sum(data == 1)
  f2 <- sum(data == 2)
  f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))

  adjust.pi <- function(data) {
    data_straight <- c(data)
    # f1 <- sum(data_straight == 1)
    # f2 <- sum(data_straight == 2)
    ## n <- sum(sapply(unique(data_straight), function(x){x*sum(data_straight == x)}))
    # n <- sum(data_straight)
    C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
    lambda <- (1-C.hat)/sum(data_straight/n*(1-data_straight/n)^n)
    # f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
    p.seen <- data_straight/n*(1-lambda*(1-data_straight/n)^n)
    p.seen <- p.seen[p.seen>0]
    p.unseen <- (1-C.hat)/f0
    list(seen = p.seen,unseen = p.unseen)
  }
  phy <- create.aili(data,row.tree,col.tree)
  g1 <- sum(phy$branch.length[phy$branch.abun==1])
  g2 <- sum(phy$branch.length[phy$branch.abun==2])
  g0 <- ifelse(g2>g1*f2/2/f1,(n-1)/n*g1^2/2/g2,(n-1)/n*g1*(f1-1)/2/(f2+1))/f0
  pi <- adjust.pi(data)
  # a_i -> p_i
  pi_matrix = data
  pi_matrix[pi_matrix>0] <- pi$seen

  # transfrom pi from S1xS1 to B1xB2
  seen_interaction_aili = create.aili(pi_matrix,row.tree,col.tree)
  unseen_interaction_aili = data.frame(branch.abun = rep(pi$unseen,f0), branch.length = g0 / f0,
                                       tgroup = "Tip",interaction = "unseen")
  p <- rbind(seen_interaction_aili,unseen_interaction_aili)%>%
    mutate(branch.abun = ifelse(tgroup == 'Root', 1, branch.abun))

  total_nodes_num = nrow(p)

  lapply(1:B, function(x){
    ai <- rbinom(total_nodes_num,n,p$branch.abun)
    # rmultinom(n = 1, size = total_nodes_num, prob =  p$branch.abun)%>%c()
    out <- cbind(ai, p[,c("branch.length", "tgroup", "interaction")])
    colnames(out)[1] = 'branch.abun'
    # out <- data.frame(branch.abun = ai,branch.length = p$branch.length,tgroup = p$tgroup,interaction = p$interaction)
    out <- out[out$branch.abun>0,]
  })
}

get.netphydiv <- function(data,q,B,row.tree = NULL,col.tree = NULL,conf) {
  plan(multisession)
  phydata <- future_lapply(data,function(x){
    create.aili(x,row.tree = row.tree,col.tree = col.tree)%>%
      filter(branch.abun > 0)
  }, future.seed=NULL)

  mle <- lapply(phydata, function(x){
      PhD:::PD.qprofile(aL = x, q = q, cal = "PD", nt = sum(x[x$tgroup == "Tip","branch.abun"]))/
      sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)
  })
  est <- lapply(phydata,function(x){
    # (ai, Lis, q, nt, cal)
    # ai = x$branch.abun;
    # Lis = x$branch.length
    # q = c(0,1,2)
    # nt = nt = sum(x[x$tgroup == "Tip",]$branch.abun)
    # cal = 'PD'
      my_PhD.q.est(ai = x$branch.abun,Lis = x$branch.length,q,nt = sum(x[x$tgroup == "Tip","branch.abun"]), cal = 'PD')/
        sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)
  })
  boot.sam <- lapply(data, function(x){
    sample.boot.phy(x,B =B,row.tree = row.tree,col.tree = col.tree)
  })
  plan(multisession)
  mle.boot <- future_lapply(boot.sam, function(x){
    lapply(x, function(y){
          PhD:::PD.qprofile(y,q,cal = "PD",nt = sum(y[y$tgroup == "Tip",]$branch.abun))/sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)
    })
  }, future.seed=NULL)
  mle.sd <- lapply(mle.boot, function(x){
    tmp <- do.call(rbind,x)
    sapply(1:length(q), function(x){
      sd(tmp[,x])
    })
  })
  est.boot <- lapply(boot.sam, function(x){
    lapply(x, function(y){
      ## NEW
      my_PhD.q.est(ai = y$branch.abun,Lis = y$branch.length,q,nt = sum(y[y$tgroup == "Tip","branch.abun"]), cal = 'PD')/
        sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)
      ## OLD
      # PhD:::PhD.q.est(y,q,datatype = "abundace",nt = sum(y[y$tgroup == "Tip",]$branch.abun))/sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)
    })
  })
  est.sd <- lapply(est.boot, function(x){
    tmp <- do.call(rbind,x)
    sapply(1:length(q), function(x){
      sd(tmp[,x])
    })
  })
  ci <- qnorm(conf/2+0.5)
  out <- c()
  for (i in 1:length(data)) {
    tmp <- rbind(data.frame(Order.q = q,Estimate = mle[[i]],method = "Empirical",
                            UCL = mle[[i]] + mle.sd[[i]]*ci,
                            LCL = mle[[i]] - mle.sd[[i]]*ci),
                 data.frame(Order.q = q,Estimate = est[[i]],method = "Estimate",
                            UCL = est[[i]] + est.sd[[i]]*ci,
                            LCL = est[[i]] - est.sd[[i]]*ci))
    tmp$Region <- names(mle)[i]
    out <- rbind(out,tmp)
  }
  return(as.data.frame(out))
}

get.netphydiv_iNE <- function(data,q,B,row.tree = NULL,col.tree = NULL,conf, knots = 40) {
  q <- unique(ceiling(q))
  ci <- qnorm(conf/2+0.5)
  inex <- function(data,q,B,row.tree = NULL,col.tree = NULL) {
    data <- as.matrix(data)
    n <- sum(data)
    m <- sort(unique(ceiling(c(seq(1,2*n,length.out = knots),n))),decreasing = F)
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
      plan(multisession)
      PD <- future_lapply(m,function(y){
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
      PD.table <- data.frame(m=m,method = ifelse(m<n,"interpolated",ifelse(n == m,"observed","extrapolated")),Order.q = x,PD = PD, PD.UCL = PD+ci * PD.sd,PD.LCL = PD - ci * PD.sd)
      out <- left_join(PD.table,sc.table)
      out
    })
    do.call(rbind,out)
  }



  # out <- future_lapply(data, FUN = function(x){
  #   inex(data = x,q,B,row.tree,col.tree)
  # })
  out <- lapply(data, function(x){
    inex(data = x,q,B,row.tree,col.tree)
  })
  out1 <- c()
  for (i in 1:length(out)) {
    tmp <- out[[i]]
    tmp$Region <- names(out)[i]
    out1 <- rbind(out1,tmp)
  }
  return(out1)
}


get.netphyhiera <- function(data,mat,q,B,row.tree = NULL,col.tree = NULL,conf){
  phydata <- lapply(data, function(x){create.aili(x,row.tree = row.tree,col.tree = col.tree)})
  boot <- sample.boot.PD(data,phydata,B,row.tree,col.tree)
  out.hi <- lapply(boot, function(x){
    hierachical.PD(x,mat)
  })
  mle.boot <- lapply(out.hi, function(x){
    lapply(x, function(y){
          PhD:::PD.qprofile(y,q,cal = "PD",nt = sum(y[y$tgroup == "Tip",]$branch.abun))/sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)
    })
  })
  est.boot <- lapply(out.hi, function(x){
    lapply(x, function(y){
          PhD:::PhD.q.est(y,q,datatype = "abundace",nt = sum(y[y$tgroup == "Tip",]$branch.abun))/sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)
    })
  })
  mle.tmp <- do.call(cbind,mle.boot)
  est.tmp <- do.call(cbind,est.boot)
  mle.boot <- lapply(1:nrow(mat), function(x){
    do.call(cbind,mle.tmp[x,])/length(unique(hiera(mat)[x,]))
  })
  est.boot <- lapply(1:nrow(mat), function(x){
    do.call(cbind,est.tmp[x,])/length(unique(hiera(mat)[x,]))
  })

  beta.mle.boot <- lapply(2:nrow(mat), function(x){
    mle.boot[[x-1]]/mle.boot[[x]]
  })

  beta.est.boot <- lapply(2:nrow(mat), function(x){
    est.boot[[x-1]]/est.boot[[x]]
  })


  hidata <- hierachical.PD(phydata,mat)
  mle <- lapply(hidata, function(x){
        PhD:::PD.qprofile(x,q,cal = "PD",nt = sum(x[x$tgroup == "Tip",]$branch.abun))/sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)
  })
  est <- lapply(hidata,function(x){
        PhD:::PhD.q.est(x,q,datatype = "abundace",nt = sum(x[x$tgroup == "Tip",]$branch.abun))/sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)
  })
  mle <- lapply(1:nrow(mat), function(x){
    mle[[x]]/length(unique(hiera(mat)[x,]))
  })
  est <- lapply(1:nrow(mat), function(x){
    est[[x]]/length(unique(hiera(mat)[x,]))
  })
  beta.mle <- lapply(2:nrow(mat), function(x){
    mle[[x-1]]/mle[[x]]
  })
  beta.est <- lapply(2:nrow(mat), function(x){
    est[[x-1]]/est[[x]]
  })
  ci <- qnorm(conf/2+0.5)
  gamma <- rbind(data.frame(Order.q=q,Estimate=mle[[1]],Method = "Empirical",
                            UCL = mle[[1]]+ci*sapply(1:length(q),function(x){sd(mle.boot[[1]][x,])}),
                            LCL = mle[[1]]-ci*sapply(1:length(q),function(x){sd(mle.boot[[1]][x,])})),
                 data.frame(Order.q=q,Estimate=est[[1]],Method = "Estimate",
                            UCL = est[[1]]+ci*sapply(1:length(q),function(x){sd(est.boot[[1]][x,])}),
                            LCL = est[[1]]-ci*sapply(1:length(q),function(x){sd(est.boot[[1]][x,])})))
  gamma <- cbind(gamma,data.frame(Type = "Gamma"))

  alpha <- lapply(2:nrow(mat), function(x){
    rbind(data.frame(Order.q=q,Estimate=mle[[x]],Method = "Empirical",
                     UCL = mle[[x]]+ci*sapply(1:length(q),function(y){sd(mle.boot[[x]][y,])}),
                     LCL = mle[[x]]-ci*sapply(1:length(q),function(y){sd(mle.boot[[x]][y,])})),
          data.frame(Order.q=q,Estimate=est[[x]],Method = "Estimate",
                     UCL = est[[x]]+ci*sapply(1:length(q),function(y){sd(est.boot[[x]][y,])}),
                     LCL = est[[x]]-ci*sapply(1:length(q),function(y){sd(est.boot[[x]][y,])})))
  })
  names(alpha) <- paste("Alpha",(nrow(mat)-1):1)
  alpha <- lapply(1:length(alpha), function(x){
    alpha[[x]] <- cbind(alpha[[x]],data.frame(Type = names(alpha)[x]))
  })

  beta <- lapply(1:(nrow(mat)-1), function(x){
    rbind(data.frame(Order.q=q,Estimate=beta.mle[[x]],Method = "Empirical",
                     UCL = beta.mle[[x]]+ci*sapply(1:length(q),function(y){sd(beta.mle.boot[[x]][y,])}),
                     LCL = beta.mle[[x]]-ci*sapply(1:length(q),function(y){sd(beta.mle.boot[[x]][y,])})),
          data.frame(Order.q=q,Estimate=beta.est[[x]],Method = "Estimate",
                     UCL = beta.est[[x]]+ci*sapply(1:length(q),function(y){sd(beta.est.boot[[x]][y,])}),
                     LCL = beta.est[[x]]-ci*sapply(1:length(q),function(y){sd(beta.est.boot[[x]][y,])})))
  })
  names(beta) <- paste("Beta",(nrow(mat)-1):1)
  beta <- lapply(1:length(beta), function(x){
    beta[[x]] <- cbind(beta[[x]],data.frame(Type = names(beta)[x]))
  })
  out <- rbind(gamma,do.call(rbind,alpha),do.call(rbind,beta))
  out
}

bootstrap_population_multiple_network = function(data, data_gamma, data_type){

  if (data_type == 'abundance'){

    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()

    output = apply(data, 2, function(x){

      p_i_hat = chaoUtility:::bootp_one_abu(Spec = x, zero = T)

      ## with unseen species (interaction)
      if(length(p_i_hat) != length(x)){

        p_i_hat_unobs = p_i_hat[length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)

        chosen = sample(x = candidate, size = min(length(p_i_hat) - length(x), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)

        p_i_hat

      } else {

        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat

      }
    })
  }

  if (data_type == 'incidence'){

    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling

    output = apply(data, 2, function(x){

      pi_i_hat = chaoUtility:::bootp_one_inc(Spec = x, zero = T)

      if(length(pi_i_hat) != length(x)){

        pi_i_hat_unobs = pi_i_hat[length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:length(x)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat==0)
        chosen = sample(x = candidate, size = min(length(pi_i_hat) - length(x), length(candidate)), replace = F)
        pi_i_hat[chosen] = pi_i_hat_unobs

        pi_i_hat

      } else {

        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat

      }
    })
  }
  return(output)
}

bootstrap_population_multiple_assemblage = function(data, data_gamma, data_type){

  if (data_type == 'abundance'){

    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()

    output = apply(data, 2, function(x){

      p_i_hat = chaoUtility:::bootp_one_abu(Spec = x, zero = T)

      if(length(p_i_hat) != length(x)){

        p_i_hat_unobs = p_i_hat[length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)

        chosen = sample(x = candidate, size = min(length(p_i_hat) - length(x), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)

        p_i_hat

      } else {

        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat

      }
    })
  }

  if (data_type == 'incidence'){

    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling

    output = apply(data, 2, function(x){

      pi_i_hat = chaoUtility:::bootp_one_inc(Spec = x, zero = T)

      if(length(pi_i_hat) != length(x)){

        pi_i_hat_unobs = pi_i_hat[length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:length(x)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat==0)
        chosen = sample(x = candidate, size = min(length(pi_i_hat) - length(x), length(candidate)), replace = F)
        pi_i_hat[chosen] = pi_i_hat_unobs

        pi_i_hat

      } else {

        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat

      }
    })
  }
  return(output)
}

