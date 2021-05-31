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
