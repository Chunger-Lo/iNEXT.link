## PHD
my_PhD.q.est <- function (ai, Lis, q, nt, cal)
{
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  I1 <- which(ai == 1)
  I2 <- which(ai == 2)
  f1 <- length(I1)
  f2 <- length(I2)
  if (f2 > 0) {
    A = 2 * f2/((nt - 1) * f1 + 2 * f2)
  }
  else if (f2 == 0 & f1 > 0) {
    A = 2/((nt - 1) * (f1 - 1) + 2)
  }
  else {
    A = 1
  }
  S <- length(ai)
  if (1 %in% q) {
    ai_h1_I <- ai <= (nt - 1)
    h1_pt2 <- rep(0, S)
    ai_h1 <- ai[ai_h1_I]
    h1_pt2[ai_h1_I] <- tibble(ai = ai) %>% .[ai_h1_I, ] %>%
      mutate(diga = digamma(nt) - digamma(ai)) %>% apply(.,
                                                         1, prod)/nt
  }
  if (2 %in% q) {
    q2_pt2 <- unlist(ai * (ai - 1)/nt/(nt - 1))
  }
  if (sum(abs(q - round(q)) != 0) > 0 | max(q) > 2) {
    deltas_pt2 <- sapply(0:(nt - 1), function(k) {
      ai_delt_I <- ai <= (nt - k)
      deltas_pt2 <- rep(0, S)
      deltas_pt2[ai_delt_I] <- delta_part2(ai = ai[ai_delt_I],
                                           k = k, n = nt)
      deltas_pt2
    }) %>% t()
  }
  Sub <- function(q, g1, g2, PD_obs, t_bar, Li) {
    if (q == 0) {
      ans <- PD_obs + iNEXTPD2:::Dq0(nt, f1, f2, g1, g2)
    }
    else if (q == 1) {
      h2 <- iNEXTPD2:::Dq1_2(nt, g1, A)
      h1 <- sum(Li * h1_pt2)
      h <- h1 + h2
      ans <- t_bar * exp(h/t_bar)
    }
    else if (q == 2) {
      ans <- t_bar^2/sum(Li * q2_pt2)
    }
    else {
      k <- 0:(nt - 1)
      deltas <- as.numeric(deltas_pt2 %*% Li)
      a <- (choose(q - 1, k) * (-1)^k * deltas) %>% sum
      b <- ifelse(g1 == 0 | A == 1, 0, (g1 * ((1 - A)^(1 -
                                                         nt))/nt) * (A^(q - 1) - sum(choose(q - 1, k) *
                                                                                       (A - 1)^k)))
      ans <- ((a + b)/(t_bar^q))^(1/(1 - q))
    }
    return(ans)
  }
  Lis = as.data.frame(Lis)
  est <- sapply(1:ncol(Lis), function(i) {
    Li = Lis[, i]
    t_bar <- t_bars[i]
    PD_obs <- sum(Li)
    g1 <- sum(Li[I1])
    g2 <- sum(Li[I2])
    est <- sapply(q, function(q_) Sub(q = q_, g1 = g1, g2 = g2,
                                      PD_obs = PD_obs, t_bar = t_bar, Li = Li))
  })
  if (cal == "PD") {
    est <- as.numeric(est)
  }
  else if (cal == "meanPD") {
    est <- as.numeric(sapply(1:length(t_bars), function(i) {
      est[, i]/t_bars[i]
    }))
  }
  return(est)
}


my_PhD.m.est <- function (ai, Lis, m, q, nt, cal)
{
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  if (sum(m > nt) > 0) {
    EPD = function(m, obs, asy) {
      m = m - nt
      Lis = as.data.frame(Lis)
      out <- sapply(1:ncol(Lis), function(i) {
        asy_i <- asy[, i]
        obs_i <- obs[, i]
        RPD_m_i <- RPD_m[, i]
        Li <- Lis[, i]
        t_bar <- t_bars[i]
        asy_i <- sapply(1:length(q), function(j) {
          max(asy_i[j], obs_i[j])
        })
        beta <- rep(0, length(q))
        beta0plus <- which(asy_i != obs_i)
        beta[beta0plus] <- (obs_i[beta0plus] - RPD_m_i[beta0plus])/(asy_i[beta0plus] -
                                                                      RPD_m_i[beta0plus])
        outq <- sapply(1:length(q), function(i) {
          if (q[i] != 2) {
            obs_i[i] + (asy_i[i] - obs_i[i]) * (1 -
                                                  (1 - beta[i])^m)
          }
          else if (q[i] == 2) {
            1/sum((Li/(t_bar)^2) * ((1/(nt + m)) * (ai/nt) +
                                      ((nt + m - 1)/(nt + m)) * (ai * (ai -
                                                                         1)/(nt * (nt - 1)))))
          }
        })
        outq
      })
      return(out)
    }
    RPD_m <- PhD:::RPD(ai%>%as.matrix(), Lis%>%as.matrix(), nt, nt - 1, q)
    obs <- PhD:::RPD(ai%>%as.matrix(), Lis%>%as.matrix(), nt, nt, q)
    asy <- matrix(my_PhD.q.est(ai = ai, Lis = Lis, q = q, nt = nt,
                            cal = cal), nrow = length(q), ncol = length(t_bars))
  }
  else if (sum(m == nt) > 0) {
    obs <- PhD:::RPD(ai%>%as.matrix(), Lis%>%as.matrix(), nt, nt, q)
  }
  if (cal == "PD") {
    out <- sapply(m, function(mm) {
      if (mm < nt) {
        ans <- PhD:::RPD(ai = ai%>%as.matrix(), Lis = Lis%>%as.matrix(), n = nt, m = mm,
                   q = q)
      }
      else if (mm == nt) {
        ans <- obs
      }
      else {
        ans <- EPD(m = mm, obs = obs, asy = asy)
      }
      return(as.numeric(ans))
    })
  }
  else if (cal == "meanPD") {
    out <- sapply(m, function(mm) {
      if (mm < nt) {
        ans <- PhD:::RPD(ai = ai%>%as.matrix(), Lis = Lis%>%as.matrix(), n = nt, m = mm,
                   q = q)
      }
      else if (mm == nt) {
        ans <- obs
      }
      else {
        ans <- EPD(m = mm, obs = obs, asy = asy)
      }
      ans <- sapply(1:length(t_bars), function(i) {
        ans[, i]/t_bars[i]
      })
      as.numeric(ans)
    })
  }
  out <- matrix(out, ncol = length(m))
  return(out)
}

# coverage_to_size <- function(x, C){
#   x <- x[x > 0]
#   m <- NULL
#   n <- sum(x)
#   refC <- iNEXT:::Chat.Ind(x, n)
#   f <- function(m, C) abs(iNEXT:::Chat.Ind(x, m) - C)
#   mm <- sapply(C, function(cvrg) {
#     if (refC == cvrg) {
#       mm <- n
#     }
#     else if (refC > cvrg) {
#       opt <- optimize(f, C = cvrg, lower = 0, upper = sum(x))
#       mm <- opt$minimum
#     }
#     else if (refC < cvrg) {
#       f1 <- sum(x == 1)
#       f2 <- sum(x == 2)
#       if (f1 > 0 & f2 > 0) {
#         A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
#       }
#       else if (f1 > 1 & f2 == 0) {
#         A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) +
#                                    2)
#       }
#       else if (f1 == 0 & f2 > 0) {
#         A <- 0
#       }
#       else if (f1 == 1 & f2 == 0) {
#         A <- 0
#       }
#       else if (f1 == 0 & f2 == 0) {
#         A <- 0
#       }
#       mm <- ifelse(A == 0, 0, (log(n/f1) + log(1 - cvrg))/log(A) -
#                      1)
#       mm <- n + mm
#     }
#     mm
#   })
#   mm[mm < 1] <- 1
#   return(mm)
# # }
Coverage_to_size <- function(x, C){
  x <- x[x > 0]
  m <- NULL
  n <- sum(x)
  refC <- iNEXT:::Chat.Ind(x, n)
  f <- function(m, C) abs(iNEXT:::Chat.Ind(x, m) - C)
  mm <- sapply(C, function(cvrg) {
    if (refC == cvrg) {
      mm <- n
    }
    else if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = sum(x))
      mm <- opt$minimum
    }
    else if (refC < cvrg) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      else if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) +
                                   2)
      }
      else if (f1 == 0 & f2 > 0) {
        A <- 0
      }
      else if (f1 == 1 & f2 == 0) {
        A <- 0
      }
      else if (f1 == 0 & f2 == 0) {
        A <- 0
      }
      mm <- ifelse(A == 0, 0, (log(n/f1) + log(1 - cvrg))/log(A) -
                     1)
      mm <- n + mm
    }
    mm
  })
  mm[mm < 1] <- 1
  return(mm)
}
