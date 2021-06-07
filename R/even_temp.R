even.class = function (q, qD, S, E.class)
{
  tmp = c()
  if (E.class == 1)
    tmp = ifelse(q != 1, (1 - qD^(1 - q))/(1 - S^(1 - q)),
                 log(qD)/log(S))
  if (E.class == 2)
    tmp = ifelse(q != 1, (1 - qD^(q - 1))/(1 - S^(q - 1)),
                 log(qD)/log(S))
  if (E.class == 3)
    tmp = (qD - 1)/(S - 1)
  if (E.class == 4)
    tmp = (1 - 1/qD)/(1 - 1/S)
  if (E.class == 5)
    tmp = log(qD)/log(S)
  return(tmp)
}

Evenness.profile <- function (x, q, datatype = c("abundance", "incidence_freq"),
          method, E.class, C = NULL)
{
  if (method == "Estimated") {
    x = long[[1]]

    C = 0.7
    estqD = estimateD.link(x, diversity = "PD", q, datatype,
                       base = "coverage", level = C, nboot = 0,
                       col.tree = col.tree, row.tree = row.tree)
    estS = estqD%>%filter(Order.q == 0 )


    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q = q,
                                                       qD = estqD%>%pull(qPD),
                                                       S = estS%>%pull(qPD),
                                                       i))
      if (class(tmp)[1] %in% c("numeric", "integer")) {
        tmp = t(as.matrix(tmp, nrow = 1))
      }
      rownames(tmp) = q
      tmp
    })
  }
  else if (method == "Empirical") {
    empqD = obs3D(x, diversity = "TD", q = q, datatype = datatype,
                  nboot = 0)
    empS = obs3D(x, diversity = "TD", q = 0, datatype = datatype,
                 nboot = 0)
    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q,
                                                       empqD[empqD$Assemblage == names(x)[k], "qD"],
                                                       empS[empS$Assemblage == names(x)[k], "qD"],
                                                       i, x[[k]]/sum(x[[k]])))
      if (class(tmp)[1] %in% c("numeric", "integer")) {
        tmp = t(as.matrix(tmp, nrow = 1))
      }
      rownames(tmp) = q
      tmp
    })
  }
  names(out) = paste("E", E.class, sep = "")
  return(out)
}
