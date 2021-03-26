
# package -----------------------------------------------------------------


library(iNEXT)
library(iNEXTPD2)
library(FunD)
library(dplyr)
library(stringr)
library(ggplot2)
library(xlsx)
library(ape)

# iNEXT3D -----------------------------------------------------------------


iNEXT3D = function(data,
                   q=c(0,1,2),
                   endpoint = NULL,
                   data_type=c('abundance', 'incidence_raw', 'incidence_freq'),
                   diversity_type=c('taxonomic', 'phylogenetic', 'functional'),
                   nboot = 0,
                   conf = 0.95,
                   PD_level=c('PD', 'meanPD'),
                   phy_tree=NULL,
                   reftime = NULL,
                   distance_matrix=NULL,
                   tau_type = c('single', 'AUC'),
                   tau=NULL)
{

  # warning -----------------------------------------------------------------


  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(data_type, TYPE))){
    stop("invalid datatype")
  }
  if(pmatch(data_type, TYPE) == -1){
    stop("ambiguous datatype")
  }
  data_type <- match.arg(data_type, TYPE)

  if(data_type == "incidence"){
    stop(' please try data_type="incidence_freq" or data_type= "incidence_raw".')
  }
  class_x = class(data)[1]

  # TD ----------------------------------------------------------------------


  if("taxonomic" %in% diversity_type){

    TD_func <- function(){

      if(data_type != "incidence_raw"){

        x <- lapply(1:ncol(data), function(i) data[, i])

        if(class_x == "matrix"){
          names(x) <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

      }

      if(data_type == "incidence_raw"){

        x <- lapply(data, iNEXT:::as.incfreq)

        if(class_x %in% c("matrix", "data.frame")){
          names(x) <- "Region_1"
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }


      }

      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      out <- iNEXT(x,
                   q = q,
                   datatype = data_type,
                   endpoint = endpoint,
                   nboot = nboot,
                   conf = conf)

      out$iNextEst$size_based$Type <- "TD"
      out$iNextEst$coverage_based$Type <- "TD"
      out$AsyEst$Type <- "TD"
      return(out)

    }
    TD_out = TD_func()

  }


  # PD ----------------------------------------------------------------------

  if("phylogenetic" %in% diversity_type){

    PD_func <- function(){
      if(data_type == "incidecne_freq"){stop("incidence_raw data needed for phylogenetic diversity")}

      if(data_type =="abundance"){

        if(class_x %in% c("matrix", "data.frame")){

          x <- as.matrix(data)
          colnames(x) = if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)

        }

        if(class_x == "list" & length(data) == 1){

          x <- as.matrix(data[[1]])

          if(is.null(colnames(x))){

            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:ncol(x)) else names(data)

          }
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) i %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          if(is.null(colnames(x))){
            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
          }
        }

        infor = PDInfo(x, tree = phy_tree, reftime = reftime)


        ### Can't put single Order.q in "meanPD"
        iNextEst = iNEXTPD(x,
                           datatype = "abundance",
                           tree = phy_tree,
                           q = c(q,0,1),
                           reftime = reftime,
                           type = PD_level,
                           endpoint = endpoint,
                           nboot = nboot,
                           conf = conf)
        iNextEst$size_based = filter(iNextEst$size_based, Order.q %in% q)
        iNextEst$coverage_based = filter(iNextEst$coverage_based, Order.q %in% q)


        asy = PhdAsy(data = x,
                     datatype = "abundance",
                     tree = phy_tree,
                     q = c(0,1,2),
                     reftime = reftime,
                     type = PD_level,
                     nboot = nboot,
                     conf = conf)

        colnames(asy)[3:5] = c("Estimator", "LCL", "UCL")

        obs = PhdObs(x,
                     datatype = "abundance",
                     tree = phy_tree,
                     q = c(0,1,2),
                     reftime = reftime,
                     type = PD_level,
                     nboot = nboot,
                     conf = conf)
        colnames(obs)[3] = "Observed"

        AsyEst=full_join(select(obs, Assemblage, Order.q, Observed),
                         select(asy, Assemblage, Order.q, Estimator, LCL, UCL, Reference.time, Type),
                         by = c("Assemblage", "Order.q"))

        out=list(DataInfo = infor,
                 iNextEst = iNextEst,
                 AsyEst = AsyEst)


      }

      if(data_type == "incidence_raw"){

        if(class_x != "list"){
          x <- as.matrix(data)
          nT <- ncol(x)
          names(nT) = "Region_1"
        }

        if(class_x == "list" & length(data) == 1){
          x <- as.matrix(data[[1]])
          rownames(x) <- rownames(data[[1]])
          nT <- ncol(x)
          names(nT) <- ifelse(is.null(names(data)), "Region_1", names(data))
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) data.frame(i) %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          nT <- sapply(data, ncol)
          names(nT) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

        }

        infor = PDInfo(x,
                       nT = nT,
                       datatype = "incidence_raw",
                       tree = phy_tree,
                       reftime = reftime)

        iNextEst = iNEXTPD(x,
                           nT = nT,
                           datatype = "incidence_raw",
                           tree = phy_tree,
                           q = q,
                           reftime = reftime,
                           type = PD_level,
                           endpoint = endpoint,
                           nboot = nboot,
                           conf = conf)

        colnames(iNextEst$size_based)[colnames(iNextEst$size_based)=="nt"] = "t"
        colnames(iNextEst$coverage_based)[colnames(iNextEst$coverage_based)=="nt"] = "t"

        asy = PhdAsy(data = x,
                     nT = nT,
                     datatype = "incidence_raw",
                     tree = phy_tree,
                     q = c(0,1,2),
                     reftime = reftime,
                     type = PD_level,
                     nboot = nboot,
                     conf = conf)

        colnames(asy)[3:5] = c("Estimator", "LCL", "UCL")

        obs = PhdObs(x,
                     nT = nT,
                     datatype = "incidence_raw",
                     tree = phy_tree,
                     q = c(0,1,2),
                     reftime = reftime,
                     type = PD_level,
                     nboot = nboot,
                     conf = conf)
        colnames(obs)[3] = "Observed"

        AsyEst=full_join(select(obs, Assemblage, Order.q, Observed),
                         select(asy, Assemblage, Order.q, Estimator, LCL, UCL, Reference.time, Type),
                         by = c("Assemblage", "Order.q"))

        out=list(DataInfo = infor,
                 iNextEst = iNextEst,
                 AsyEst = AsyEst)
      }
      return(out)
    }

    PD_out <- PD_func()

  }

  # FD ----------------------------------------------------------------------


  if("functional" %in% diversity_type){

    FD_func <- function(){

      if(data_type != "incidence_raw"){

        if(class_x == "matrix"){
          names_data <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names_data <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

        x <- lapply(1:ncol(data), function(i) {
          out = data.frame(data[,i])
          rownames(out) = rownames(data)
          names(out) <- names_data[i]
          return(out)
        })



      }

      if(data_type == "incidence_raw"){



        if(class_x %in% c("matrix", "data.frame")){

          x <- lapply(data, function(i){
            out <- data.frame(iNEXT:::as.incfreq(i))
            names(out) <- "Region_1"
            return(out)
          })


        }else{

          names_x <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

          x <- lapply(1:length(data), function(i){
            out <- data.frame(iNEXT:::as.incfreq(data[[i]]))
            names(out) <- names_x[i]
            return(out)
          } )
        }


      }



      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      if(tau_type == "single"){

        coverage_based <- lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]

          if(data_type == "abundance"){
            index <- match(rownames(distance_matrix),rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[c(index)]
          }else{
            #index <- c(1, (sapply(rownames(i)[-1], function(k) which(k == rownames(distance_matrix)))+1)) %>% as.vector()
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[c(index),])
            rownames(i) = name[c(index)]
          }

          names(i) = ass

          if(data_type == "abundance"){
            cmin = iNEXT:::Chat.Ind(as.matrix(i), sum(as.matrix(i)))
            cmax = iNEXT:::Chat.Ind(as.matrix(i), 2*sum(as.matrix(i)))

          }else{
            cmin = iNEXT:::Chat.Sam(as.matrix(i), as.matrix(i)[1,])
            cmax = iNEXT:::Chat.Sam(as.matrix(i), 2*as.matrix(i)[1,])

          }

          out = FunD::EstimateFD(i,
                                 datatype = data_type,
                                 distM = as.matrix(distance_matrix),
                                 q = q,
                                 nboot = nboot,
                                 conf = conf,
                                 threshold = tau,
                                 level = c(seq(0.5, 0.95, 0.05), seq(cmin,cmax,0.1),cmax))

          colnames(out)[1] = "Assemblage"

          if(data_type == "incidence_freq") colnames(out)[colnames(out) == "m"] = "t"
          colnames(out)[colnames(out) == "method"] = "Method"
          out$Method[out$Method == "Interpolated"] = "Rarefaction"
          out$Method[out$Method == "Extrapolated"] = "Extrapolation"
          out$Method[out$goalSC==cmin] = "Observed"

          return(out)
        } ) %>% do.call(rbind,.)

        size_based=lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- match(rownames(distance_matrix),rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[c(index)]
          }else{
            #index <- c(1, (sapply(rownames(i)[-1], function(k) which(k == rownames(distance_matrix)))+1))
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }

          names(i) = ass


          out = (FunD::iNEXTFD(data = i,
                               distM = as.matrix(distance_matrix),
                               datatype = data_type,
                               nboot = nboot,
                               conf = conf,
                               threshold = tau,
                               q = q,
                               endpoint = endpoint))$inext

          colnames(out)[colnames(out) == "Community"] = "Assemblage"

          if(data_type == "incidence_freq") colnames(out)[colnames(out) == "m"] = "t"
          colnames(out)[colnames(out) == "method"] = "Method"
          out$Method[out$Method == "Interpolated"] = "Rarefaction"
          out$Method[out$Method == "Extrapolated"] = "Extrapolation"

          return(out)
        })%>% do.call(rbind,.)

        AsyEst=lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau = if(is.null(tau)) threhold_func(i, distance_matrix)[2] else tau
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            #index <- c(1, (sapply(rownames(i)[-1], function(k) which(k == rownames(distance_matrix)))+1))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau = if(is.null(tau)) threhold_func(i[-1,], distance_matrix)[2] else tau

          }

          names(i) = ass

          obs = FunD:::FDtable_mle(i,
                                   datatype = data_type,
                                   dij = as.matrix(distance_matrix),
                                   q = c(0,1,2),
                                   nboot = 0,
                                   conf = conf,
                                   tau = tau)
          colnames(obs)[colnames(obs) == "Empirical"] = "Observed"
          asy_info=FunD:::FDtable_est(i,
                                      datatype = data_type,
                                      dij = as.matrix(distance_matrix),
                                      q = c(0,1,2),
                                      nboot = nboot,
                                      conf = conf,
                                      tau = tau)
          asy = asy_info$Estoutput
          #info = asy_info$info

          out = full_join(select(obs, "Order.q", "Observed", "Community"),
                          asy,
                          by = c("Order.q", "Community"))

          colnames(out)[3:4] = c("Assemblage", "Estimator")
          out = out[c("Assemblage", "Order.q", "Observed", "Estimator", "LCL", "UCL", "tau")]

          list(out = out, info = asy_info$info)
        })

        DataInfo = do.call(rbind, lapply(AsyEst, function(i) i$info))
        colnames(DataInfo)[colnames(DataInfo) == "Community"] = "Assemblage"


        AsyEst = do.call(rbind, lapply(AsyEst, function(i) i$out))

        colnames(AsyEst)[colnames(AsyEst) == "tau"] = "threshold"
        coverage_based$Type = "FD"
        AsyEst$Type = "FD"

        out = list(DataInfo = DataInfo,
                   iNextEst=list(size_based = size_based,
                                 coverage_based = coverage_based),
                   AsyEst = AsyEst)

      }
      if(tau_type == "AUC"){

        DataInfo=lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau = if(is.null(tau)) threhold_func(i, distance_matrix)[2] else tau

          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau = if(is.null(tau)) threhold_func(i[-1,], distance_matrix)[2] else tau

          }

          names(i) = ass


          FunD:::FDtable_est(i,
                             datatype = data_type,
                             dij = as.matrix(distance_matrix),
                             q = c(0,1,2),
                             nboot = 0,
                             conf = conf,
                             tau = tau)$info
        }) %>% do.call(rbind,.)

        colnames(DataInfo)[1:2] = c("Assemblage", "Threshole")


        coverage_based = lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }

          names(i) = ass

          if(data_type == "abundance"){
            cmin = iNEXT:::Chat.Ind(as.matrix(i), sum(as.matrix(i)))
            cmax = iNEXT:::Chat.Ind(as.matrix(i), 2*sum(as.matrix(i)))

          }else{
            cmin = iNEXT:::Chat.Sam(as.matrix(i), as.matrix(i)[1,])
            cmax = iNEXT:::Chat.Sam(as.matrix(i), 2*as.matrix(i)[1,])

          }

          out = FunD:::AUCtable_invFD(i,
                                      dij = as.matrix(distance_matrix),
                                      q = q,
                                      datatype = data_type,
                                      level = c(seq(0.5,0.9,0.05),seq(cmin, cmax, 0.1),cmax),
                                      nboot = nboot,
                                      conf = conf)
          out$method[out$goalSC==cmin]="Observed"
          return(out)

        }) %>% do.call(rbind,.)


        coverage_based = select(coverage_based, "Community":"SC")

        colnames(coverage_based)[colnames(coverage_based) == "m"] = ifelse(data_type == "abundance", "m", "t")
        colnames(coverage_based)[colnames(coverage_based) == "Community"] = "Assemblage"
        colnames(coverage_based)[colnames(coverage_based) == "method"] = "Method"
        coverage_based$Method[coverage_based$Method == "Interpolated"] = "Rarefaction"
        coverage_based$Method[coverage_based$Method == "Extrapolated"] = "Extrapolation"

        size_based=lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }

          names(i) = ass

          if(data_type == "incidence_freq"){
            m = sort(unique(c(1, c(i[[1]])[1], 2*c(i[[1]])[1], round(seq(1, 2*c(i[[1]])[1], (2*c(i[[1]])[1]-1)/10)))))
          }else{
            m = sort(unique(c(1, sum(c(i[[1]])), 2*c(i[[1]])[1], round(seq(1, 2*sum(c(i[[1]])), (2*sum(c(i[[1]]))-1)/10)))))
          }

          out = lapply(m, function(size) FunD:::AUCtable_iNextFD(i,
                                                                 dij=as.matrix(distance_matrix),
                                                                 q = q,
                                                                 nboot= nboot,
                                                                 datatype = data_type,
                                                                 m = size)) %>% do.call(rbind,.)
          colnames(out)[colnames(out) == "Community"] = "Assemblage"

          if(data_type == "incidence_freq") colnames(out)[colnames(out) == "m"] = "t"
          colnames(out)[colnames(out) == "method"] = "Method"

          return(out)
        })%>% do.call(rbind,.)


        asy = lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }

          names(i) = ass

          out = FunD:::AUCtable_est(i,
                                    dij = as.matrix(distance_matrix),
                                    q = c(0,1,2),
                                    datatype = data_type,
                                    nboot = nboot,
                                    conf = conf)

        }) %>% do.call(rbind,.)

        colnames(asy)[colnames(asy) == "AUC"] = "Estimator"


        obs = lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }

          names(i) = ass

          out = FunD:::AUCtable_mle(i,
                                    dij = as.matrix(distance_matrix),
                                    q = c(0,1,2),
                                    datatype = data_type,
                                    nboot = nboot,
                                    conf = conf)

        }) %>% do.call(rbind,.)

        colnames(obs)[colnames(obs) == "AUC"] = "Observed"


        AsyEst = full_join(select(obs, "Community", "Order.q", "Observed"),
                           asy,
                           by = c("Community", "Order.q"))

        colnames(AsyEst)[colnames(AsyEst) == "Community"] = "Assemblage"

        coverage_based$Type = "FD"
        AsyEst$Type = "FD"

        out = list(DataInfo = DataInfo,
                   iNextEst = list(size_based = size_based,
                                   coverage_based = coverage_based),
                   AsyEst = AsyEst)


      }
      return(out)

    }
    FD_out <- FD_func()

  }

  match_index <- pmatch(diversity_type, c("taxonomic","phylogenetic","functional"))

  if(length(match_index) == 1){
    if(match_index == 1){
      out <- list(TD = TD_out)
    }
    if(match_index == 2){
      out <- list(PD = PD_out)
    }
    if(match_index == 3){
      out <- list(FD = FD_out)
    }

  }else if(length(match_index) == 2){
    if(sum(sort(match_index) == c(1,2)) == 2){
      out <- list(TD = TD_out,
                  PD = PD_out)
    }
    if(sum(sort(match_index) %in% c(1,3)) == 2){
      out <- list(TD = TD_out,
                  FD = FD_out)
    }
    if(sum(sort(match_index) == c(2,3)) == 2){
      out <- list(PD = PD_out,
                  FD = FD_out)
    }
  }else{
    if(sum(sort(match_index) == c(1,2,3)) == 3){
      out <- list(TD = TD_out,
                  PD = PD_out,
                  FD = FD_out)
    }


  }



  return(out)

}


# estimate3D --------------------------------------------------------------

estimate3D <- function(data,
                       q = c(0, 1, 2),
                       data_type = "abundance",
                       diversity_type = c("taxonomic", "phylogenetic", "functional"),
                       base = "coverage",
                       level = NULL,
                       nboot = 0,
                       conf = 0.95,
                       phy_tree = NULL,
                       reftime = NULL,
                       PD_level=c('PD', 'meanPD'),
                       distance_matrix = NULL,
                       tau_type = c("single", "AUC"),
                       tau = NULL)
{

  # warning -----------------------------------------------------------------


  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(data_type, TYPE))){
    stop("invalid datatype")
  }
  if(pmatch(data_type, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(data_type, TYPE)
  class_x <- class(data)[1]

  if(data_type == "incidence"){
    stop(' please try data_type="incidence_freq" or data_type= "incidence_raw".')
  }
  class_x = class(data)[1]

  # TD ----------------------------------------------------------------------


  if("taxonomic" %in% diversity_type){


    TD_func <- function(){

      if(data_type != "incidence_raw"){
        data = as.data.frame(data)
        x <- lapply(1:ncol(data), function(i) data[, i])

        if(class_x == "matrix"){
          names(x) <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

      }

      if(data_type == "incidence_raw"){

        x <- lapply(data, iNEXT:::as.incfreq)

        if(class_x %in% c("matrix", "data.frame")){
          names(x) <- "Region_1"
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }


      }

      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      out <- estimateD(x,
                       q = q,
                       base = base,
                       level = level,
                       nboot = nboot,
                       conf = conf,
                       datatype = data_type)
      rownames(out) = 1:nrow(out)
      out$Type <- "TD"
      return(out)

    }
    TD_out <- TD_func()


  }


  # PD ----------------------------------------------------------------------

  if("phylogenetic" %in% diversity_type){

    if(data_type == "incidecne_freq"){stop("incidence_raw data needed for phylogenetic diversity")}

    PD_func <- function(){

      if(data_type =="abundance"){

        if(class_x %in% c("matrix", "data.frame")){

          x <- as.matrix(data)
          colnames(x) = if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)

        }

        if(class_x == "list" & length(data) == 1){

          x <- as.matrix(data[[1]])

          if(is.null(colnames(x))){

            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:ncol(x)) else names(data)

          }
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) i %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          if(is.null(colnames(x))){
            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
          }
        }

        out <- estimatePD(data = x,
                          tree = phy_tree,
                          q = q,
                          type = PD_level,
                          reftime = NULL,
                          nboot = nboot,
                          conf = conf,
                          level = level
        )

      }

      if(data_type == "incidence_raw"){

        if(class_x != "list"){
          x <- as.matrix(data)
          nT <- ncol(x)
          names(nT) = "Region_1"
        }

        if(class_x == "list" & length(data) == 1){
          x <- as.matrix(data[[1]])
          rownames(x) <- rownames(data[[1]])
          nT <- ncol(x)
          names(nT) <- ifelse(is.null(names(data)), "Region_1", names(data))
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) data.frame(i) %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          nT <- sapply(data, ncol)
          names(nT) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

        }

        out <- estimatePD(data = x,
                          tree = phy_tree,
                          q = q,
                          type = PD_level,
                          reftime = NULL,
                          nboot = nboot,
                          conf = conf,
                          datatype = "incidence_raw",
                          nT = nT,
                          level = level)
        colnames(out)[2] =ifelse(data_type == "abundance", "m", "t")
      }
      return(data.frame(out))
    }
    PD_out <- PD_func()


  }

  # FD ----------------------------------------------------------------------


  if("functional" %in% diversity_type){

    FD_func <- function(){

      if(data_type != "incidence_raw"){

        if(class_x == "matrix"){
          names_data <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names_data <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

        x <- lapply(1:ncol(data), function(i) {
          out = data.frame(data[,i])
          rownames(out) = rownames(data)
          names(out) <- names_data[i]
          return(out)
        })



      }

      if(data_type == "incidence_raw"){



        if(class_x %in% c("matrix", "data.frame")){

          x <- lapply(data, function(i){
            out <- data.frame(iNEXT:::as.incfreq(i))
            names(out) <- "Region_1"
            return(out)
          })


        }else{

          names_x <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

          x <- lapply(1:length(data), function(i){
            out <- data.frame(iNEXT:::as.incfreq(data[[i]]))
            names(out) <- names_x[i]
            return(out)
          } )
        }


      }



      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      if(tau_type == "single"){
        out <- lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }

          names(i) = ass

          out <- data.frame(FunD:::EstimateFD(data = i,
                                              q = q,
                                              distM = as.matrix(distance_matrix),
                                              datatype = data_type,
                                              level = level,
                                              nboot = nboot,
                                              conf = conf,
                                              threshold = tau))
          colnames(out)[c(1, 2, 3)] <- c("Assemblage", ifelse(data_type == "abundance", "m", "t"), "Method")
          return(out)
        } ) %>% do.call(rbind,.)
        out <- data.frame(out)
        rownames(out) = 1:nrow(out)

      }
      if(tau_type == "AUC"){
        out <- lapply(x, function(i){
          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){

            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))

            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }
          names(i) = ass

          out <- FunD:::AUCtable_invFD(i,
                                       q = q,
                                       dij = as.matrix(distance_matrix),
                                       level = level,
                                       nboot = nboot,
                                       conf = conf,
                                       datatype = data_type)

          colnames(out)[c(1, 4, 5)] <- c("Assemblage", ifelse(data_type == "abundance", "m", "t"), "Method")
          out <- out[!(colnames(out) %in% c("SC.UCL", "SC.LCL"))]
          return(out)
        } ) %>% do.call(rbind,.)
        out <- data.frame(out)
        out <- out[c(1,4,5,2,9,6:8,3)]
        rownames(out) = 1:nrow(out)
      }
      out$Type = ifelse(tau_type == "AUC", "FD (AUC)", "FD")
      out$Method[out$Method=="Interpolated"] <- "Rarefaction"
      out$Method[out$Method=="Extrapolated"] <- "Extrapolation"

      return(out)

    }

    FD_out <- FD_func()

  }
  match_index <- pmatch(diversity_type, c("taxonomic","phylogenetic","functional"))

  if(length(match_index) == 1){
    if(match_index == 1){
      out <- list(TD = TD_out)
    }
    if(match_index == 2){
      out <- list(PD = PD_out)
    }
    if(match_index == 3){
      out <- list(FD = FD_out)
    }

  }else if(length(match_index) == 2){
    if(sum(sort(match_index) == c(1,2)) == 2){
      out <- list(TD = TD_out,
                  PD = PD_out)
    }
    if(sum(sort(match_index) %in% c(1,3)) == 2){
      out <- list(TD = TD_out,
                  FD = FD_out)
    }
    if(sum(sort(match_index) == c(2,3)) == 2){
      out <- list(PD = PD_out,
                  FD = FD_out)
    }
  }else{
    if(sum(sort(match_index) == c(1,2,3)) == 3){
      out <- list(TD = TD_out,
                  PD = PD_out,
                  FD = FD_out)
    }


  }
  return(out)
}


# Asy3D -------------------------------------------------------------------

Asy3D <- function(data,
                  q = seq(0, 2, 0.2),
                  data_type = "abundance",
                  diversity_type = c("taxonomic", "phylogenetic", "functional"),
                  nboot = 0,
                  conf = 0.95,
                  phy_tree = NULL,
                  reftime = NULL,
                  PD_level=c('PD', 'meanPD'),
                  distance_matrix = NULL,
                  tau_type = c("single", "AUC"),
                  tau = NULL)
{

  # warning -----------------------------------------------------------------


  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(data_type, TYPE))){
    stop("invalid datatype")
  }
  if(pmatch(data_type, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(data_type, TYPE)
  class_x <- class(data)[1]

  if(data_type == "incidence"){
    stop(' please try data_type="incidence_freq" or data_type= "incidence_raw".')
  }
  class_x = class(data)[1]

  # TD ----------------------------------------------------------------------


  if("taxonomic" %in% diversity_type){


    TD_func <- function(){

      if(data_type != "incidence_raw"){

        x <- lapply(1:ncol(data), function(i) data[, i])

        if(class_x == "matrix"){
          names(x) <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

      }

      if(data_type == "incidence_raw"){

        x <- lapply(data, iNEXT:::as.incfreq)

        if(class_x %in% c("matrix", "data.frame")){
          names(x) <- "Region_1"
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }


      }

      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      out <- AsyD(x,
                  q = q,
                  nboot = nboot,
                  conf = conf,
                  datatype = data_type)
      rownames(out) = 1:nrow(out)
      out$Type <- "TD"
      colnames(out)[colnames(out) == "method"] <- "Method"
      return(out)

    }
    TD_out <- TD_func()
  }


  # PD ----------------------------------------------------------------------

  if("phylogenetic" %in% diversity_type){

    if(data_type == "incidecne_freq"){stop("incidence_raw data needed for phylogenetic diversity")}

    PD_func <- function(){

      if(data_type =="abundance"){

        if(class_x %in% c("matrix", "data.frame")){

          x <- as.matrix(data)
          colnames(x) = if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)

        }

        if(class_x == "list" & length(data) == 1){

          x <- as.matrix(data[[1]])

          if(is.null(colnames(x))){

            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:ncol(x)) else names(data)

          }
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) i %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          if(is.null(colnames(x))){
            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
          }
        }

        asy <- PhdAsy(data = x,
                      tree = phy_tree,
                      q = q,
                      type = PD_level,
                      reftime = NULL,
                      nboot = nboot,
                      conf = conf,
                      datatype = "abundance")

        obs <- PhdObs(data = x,
                      tree = phy_tree,
                      q = q,
                      type = PD_level,
                      reftime = NULL,
                      nboot = nboot,
                      conf = conf,
                      datatype = "abundance")
        out <- bind_rows(asy, obs)

      }

      if(data_type == "incidence_raw"){

        if(class_x != "list"){
          x <- as.matrix(data)
          nT <- ncol(x)
          names(nT) = "Region_1"
        }

        if(class_x == "list" & length(data) == 1){
          x <- as.matrix(data[[1]])
          rownames(x) <- rownames(data[[1]])
          nT <- ncol(x)
          names(nT) <- ifelse(is.null(names(data)), "Region_1", names(data))
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) data.frame(i) %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          nT <- sapply(data, ncol)
          names(nT) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

        }

        asy <- PhdAsy(data = x,
                      tree = phy_tree,
                      q = q,
                      type = PD_level,
                      reftime = NULL,
                      nboot = nboot,
                      conf = conf,
                      datatype = "incidence_raw",
                      nT = nT)
        obs <- PhdObs(data = x,
                      tree = phy_tree,
                      q = q,
                      type = PD_level,
                      reftime = NULL,
                      nboot = nboot,
                      conf = conf,
                      datatype = "incidence_raw",
                      nT = nT)
        out <- bind_rows(asy, obs)
      }

      out <- out[c(2:5, 1, 6:8)]
      out$Method <- ifelse(out$Method == "Asymptotic", "Estimated", "Empirical")
      return(data.frame(out))
    }

    PD_out <- PD_func()

  }

  # FD ----------------------------------------------------------------------


  if("functional" %in% diversity_type){

    FD_func <- function(){

      if(data_type != "incidence_raw"){

        if(class_x == "matrix"){
          names_data <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names_data <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

        x <- lapply(1:ncol(data), function(i) {
          out = data.frame(data[,i])
          rownames(out) = rownames(data)
          names(out) <- names_data[i]
          return(out)
        })



      }

      if(data_type == "incidence_raw"){



        if(class_x %in% c("matrix", "data.frame")){

          x <- lapply(data, function(i){
            out <- data.frame(iNEXT:::as.incfreq(i))
            names(out) <- "Region_1"
            return(out)
          })


        }else{

          names_x <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

          x <- lapply(1:length(data), function(i){
            out <- data.frame(iNEXT:::as.incfreq(data[[i]]))
            names(out) <- names_x[i]
            return(out)
          } )
        }


      }



      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      if(tau_type == "single"){

        out <- lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))
            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau <- if(is.null(tau)) threhold_func(as.matrix(i),distance_matrix)[2] else tau

          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau <- if(is.null(tau)) threhold_func((i)[-1,],distance_matrix)[2] else tau

          }

          names(i) = ass


          asy <- select(data.frame(FunD:::FDtable_est(data = i,
                                                      q = q,
                                                      dij = as.matrix(distance_matrix),
                                                      datatype = data_type,
                                                      nboot = nboot,
                                                      conf = conf,
                                                      tau = tau)), starts_with("Estoutput"))
          obs <- (data.frame(FunD:::FDtable_mle(data = i,
                                                q = q,
                                                dij = as.matrix(distance_matrix),
                                                datatype = data_type,
                                                nboot = nboot,
                                                conf = conf,
                                                tau = tau)))

          colnames(asy) <- c("Order.q", "qFD", "qFD.LCL", "qFD.UCL", "threshold", "Assemblage")
          colnames(obs) <- c("Order.q", "qFD", "qFD.LCL", "qFD.UCL", "threshold", "Assemblage")

          asy$Method <- "Estimated"
          obs$Method <- "Empirical"


          out <- bind_rows(asy, obs)
          out[c("qFD", "qFD.LCL", "qFD.UCL")]= apply(out[c("qFD", "qFD.LCL", "qFD.UCL")], 2, as.vector)

          return(out)
        } ) %>% do.call(rbind,.)
        out <- data.frame(out)
        rownames(out) = 1:nrow(out)

      }
      if(tau_type == "AUC"){
        out <- lapply(x, function(i){
          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }
          names(i) = ass

          asy <- FunD:::AUCtable_est(i,
                                     q = q,
                                     dij = as.matrix(distance_matrix),
                                     nboot = nboot,
                                     conf = conf,
                                     datatype = data_type)

          obs <- FunD:::AUCtable_mle(i,
                                     q = q,
                                     dij = as.matrix(distance_matrix),
                                     nboot = nboot,
                                     conf = conf,
                                     datatype = data_type)
          obs$Method <- "Empirical"
          asy$Method <- "Estimated"
          out <- bind_rows(asy, obs)
          colnames(out) <- c("Assemblage", "Order.q", "AUC", "AUC.LCL", "AUC.UCL","Method")
          out <- out[,c(2:5,1,6)]
          return(out)
        } ) %>% do.call(rbind,.)
        out <- data.frame(out)
        rownames(out) = 1:nrow(out)
      }
      out$Type = ifelse(tau_type == "AUC", "FD (AUC)", "FD")
      return(out)

    }
    FD_out <- FD_func()

  }

  match_index <- pmatch(diversity_type, c("taxonomic","phylogenetic","functional"))

  if(length(match_index) == 1){
    if(match_index == 1){
      out <- list(TD = TD_out)
    }
    if(match_index == 2){
      out <- list(PD = PD_out)
    }
    if(match_index == 3){
      out <- list(FD = FD_out)
    }

  }else if(length(match_index) == 2){
    if(sum(sort(match_index) == c(1,2)) == 2){
      out <- list(TD = TD_out,
                  PD = PD_out)
    }
    if(sum(sort(match_index) %in% c(1,3)) == 2){
      out <- list(TD = TD_out,
                  FD = FD_out)
    }
    if(sum(sort(match_index) == c(2,3)) == 2){
      out <- list(PD = PD_out,
                  FD = FD_out)
    }
  }else{
    if(sum(sort(match_index) == c(1,2,3)) == 3){
      out <- list(TD = TD_out,
                  PD = PD_out,
                  FD = FD_out)
    }


  }
  return(out)
}

# Obs3D -------------------------------------------------------------------

# Obs3D -------------------------------------------------------------------
# data = as.numeric(bootstrap_data_gamma)
Obs3D <- function(data,
                  q = seq(0, 2, 0.2),
                  data_type = "abundance",
                  diversity_type = c("taxonomic", "phylogenetic", "functional"),
                  nboot = 0,
                  conf = 0.95,
                  phy_tree = NULL,
                  reftime = NULL,
                  PD_level=c('PD', 'meanPD'),
                  distance_matrix = NULL,
                  tau_type = c("single", "AUC"),
                  tau = NULL)
{

  # warning -----------------------------------------------------------------


  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(data_type, TYPE))){
    stop("invalid datatype")
  }
  if(pmatch(data_type, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(data_type, TYPE)
  class_x <- class(data)[1]

  if(data_type == "incidence"){
    stop(' please try data_type="incidence_freq" or data_type= "incidence_raw".')
  }
  class_x = class(data)[1]
  data = as.data.frame(data)
  # TD ----------------------------------------------------------------------


  if("taxonomic" %in% diversity_type){


    TD_func <- function(){

      if(data_type != "incidence_raw"){

        x <- lapply(1:ncol(data), function(i) data[, i])

        if(class_x == "matrix"){
          names(x) <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

      }

      if(data_type == "incidence_raw"){

        x <- lapply(data, iNEXT:::as.incfreq)

        if(class_x %in% c("matrix", "data.frame")){
          names(x) <- "Region_1"
        }else{
          names(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }


      }

      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      out <- AsyD(x,
                  q = q,
                  nboot = nboot,
                  conf = conf,
                  datatype = data_type)
      rownames(out) = 1:nrow(out)
      out$Type <- "TD"
      colnames(out)[colnames(out) == "method"] <- "Method"
      return(out %>% filter(Method == "Empirical"))

    }
    TD_out <- TD_func()
  }


  # PD ----------------------------------------------------------------------

  if("phylogenetic" %in% diversity_type){

    if(data_type == "incidecne_freq"){stop("incidence_raw data needed for phylogenetic diversity")}

    PD_func <- function(){

      if(data_type =="abundance"){

        if(class_x %in% c("matrix", "data.frame")){

          x <- as.matrix(data)
          colnames(x) = if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)

        }

        if(class_x == "list" & length(data) == 1){

          x <- as.matrix(data[[1]])

          if(is.null(colnames(x))){

            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:ncol(x)) else names(data)

          }
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) i %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          if(is.null(colnames(x))){
            colnames(x) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
          }
        }



        obs <- PhdObs(data = x,
                      tree = phy_tree,
                      q = q,
                      type = PD_level,
                      reftime = NULL,
                      nboot = nboot,
                      conf = conf,
                      datatype = "abundance")
        out <- obs

      }

      if(data_type == "incidence_raw"){

        if(class_x != "list"){
          x <- as.matrix(data)
          nT <- ncol(x)
          names(nT) = "Region_1"
        }

        if(class_x == "list" & length(data) == 1){
          x <- as.matrix(data[[1]])
          rownames(x) <- rownames(data[[1]])
          nT <- ncol(x)
          names(nT) <- ifelse(is.null(names(data)), "Region_1", names(data))
        }

        if(class_x == "list" & length(data) > 1){

          x_list <- lapply(data, function(i) data.frame(i) %>% mutate(species = rownames(i)))
          x <- x_list[[1]]
          for(i in 2:length(x_list)){

            x <- full_join(x, x_list[[i]], by = "species")

          }

          rownames(x) <- x$species
          x <- x[colnames(x)!="species"]
          x[is.na(x)] <- 0

          nT <- sapply(data, ncol)
          names(nT) <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

        }


        obs <- PhdObs(data = x,
                      tree = phy_tree,
                      q = q,
                      type = PD_level,
                      reftime = NULL,
                      nboot = nboot,
                      conf = conf,
                      datatype = "incidence_raw",
                      nT = nT)
        out <- obs
      }

      out <- out[c(2:5, 1, 6:8)]
      out$Method <- ifelse(out$Method == "Asymptotic", "Estimated", "Empirical")
      return(data.frame(out))
    }

    PD_out <- PD_func()

  }

  # FD ----------------------------------------------------------------------


  if("functional" %in% diversity_type){

    FD_func <- function(){

      if(data_type != "incidence_raw"){

        if(class_x == "matrix"){
          names_data <- if(is.null(colnames(data))) paste0("Region_", 1:ncol(data)) else colnames(data)
        }else{
          names_data <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)
        }

        x <- lapply(1:ncol(data), function(i) {
          out = data.frame(data[,i])
          rownames(out) = rownames(data)
          names(out) <- names_data[i]
          return(out)
        })



      }

      if(data_type == "incidence_raw"){



        if(class_x %in% c("matrix", "data.frame")){

          x <- lapply(data, function(i){
            out <- data.frame(iNEXT:::as.incfreq(i))
            names(out) <- "Region_1"
            return(out)
          })


        }else{

          names_x <- if(is.null(names(data))) paste0("Region_", 1:length(data)) else names(data)

          x <- lapply(1:length(data), function(i){
            out <- data.frame(iNEXT:::as.incfreq(data[[i]]))
            names(out) <- names_x[i]
            return(out)
          } )
        }


      }



      data_type <- ifelse(data_type=="abundance", "abundance", "incidence_freq")

      if(tau_type == "single"){

        out <- lapply(x, function(i){

          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))
            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau <- if(is.null(tau)) threhold_func(as.matrix(i),distance_matrix)[2] else tau

          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
            tau <- if(is.null(tau)) threhold_func((i)[-1,],distance_matrix)[2] else tau

          }

          names(i) = ass



          obs <- (data.frame(FunD:::FDtable_mle(data = i,
                                                q = q,
                                                dij = as.matrix(distance_matrix),
                                                datatype = data_type,
                                                nboot = nboot,
                                                conf = conf,
                                                tau = tau)))

          colnames(obs) <- c("Order.q", "qFD", "qFD.LCL", "qFD.UCL", "threshold", "Assemblage")

          obs$Method <- "Empirical"


          out <- obs
          out[c("qFD", "qFD.LCL", "qFD.UCL")]= apply(out[c("qFD", "qFD.LCL", "qFD.UCL")], 2, as.vector)

          return(out)
        } ) %>% do.call(rbind,.)
        out <- data.frame(out)
        rownames(out) = 1:nrow(out)

      }
      if(tau_type == "AUC"){
        out <- lapply(x, function(i){
          name = rownames(i)
          ass = names(i)
          distance_matrix <- distance_matrix[rownames(distance_matrix) %in% rownames(i),rownames(distance_matrix) %in% rownames(i)]
          if(data_type == "abundance"){
            index <- order(rownames(distance_matrix), rownames(i))
            #index <- sapply(rownames(i), function(k) which(k == rownames(distance_matrix)))
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }else{
            index <- c(1,(order(rownames(i)[-1],rownames(distance_matrix)))+1)
            i <- data.frame(i[index,])
            rownames(i) = name[index]
          }
          names(i) = ass



          obs <- FunD:::AUCtable_mle(i,
                                     q = q,
                                     dij = as.matrix(distance_matrix),
                                     nboot = nboot,
                                     conf = conf,
                                     datatype = data_type)
          obs$Method <- "Empirical"
          out <- obs
          colnames(out) <- c("Assemblage", "Order.q", "AUC", "AUC.LCL", "AUC.UCL","Method")
          out <- out[,c(2:5,1,6)]
          return(out)
        } ) %>% do.call(rbind,.)
        out <- data.frame(out)
        rownames(out) = 1:nrow(out)
      }
      out$Type = ifelse(tau_type == "AUC", "FD (AUC)", "FD")
      return(out)

    }
    FD_out <- FD_func()

  }

  match_index <- pmatch(diversity_type, c("taxonomic","phylogenetic","functional"))

  if(length(match_index) == 1){
    if(match_index == 1){
      out <- list(TD = TD_out)
    }
    if(match_index == 2){
      out <- list(PD = PD_out)
    }
    if(match_index == 3){
      out <- list(FD = FD_out)
    }

  }else if(length(match_index) == 2){
    if(sum(sort(match_index) == c(1,2)) == 2){
      out <- list(TD = TD_out,
                  PD = PD_out)
    }
    if(sum(sort(match_index) %in% c(1,3)) == 2){
      out <- list(TD = TD_out,
                  FD = FD_out)
    }
    if(sum(sort(match_index) == c(2,3)) == 2){
      out <- list(PD = PD_out,
                  FD = FD_out)
    }
  }else{
    if(sum(sort(match_index) == c(1,2,3)) == 3){
      out <- list(TD = TD_out,
                  PD = PD_out,
                  FD = FD_out)
    }


  }
  return(out)
}

# subfunction -------------------------------------------------------------

threhold_func=function(data,disM){
  tmp <- apply(apply(as.matrix(data), 2, function(i) i/sum(i)), 1, mean)
  dmean <- sum ( (tmp %*% t(tmp) ) * disM)
  dmin <- min(disM[disM>0])
  dmax <- max(disM)
  taus <- dmean
  c(dmin,taus,dmax)
}


