rm(list = ls())
library(mvtnorm)
library(glmnet)
library(survival)
source("coxph_as.R")
library(emulator)
library(MASS)
library(parallel)
numCores <- detectCores()
get_sd_stdnormal <- function(ddd) {
  sigma_seq <- seq(0.1, 10, by = 0.02)
  target_index <- which.min(abs(pnorm(IQR(ddd) / 2, mean = 0, sd = sigma_seq) - 0.75))
  return(sigma_seq[target_index])
}
generate_synthetic_X <- function(X_obs, M) {
  apply(X_obs, 2, function(v) {
    u <- sample(v, size = M, replace = TRUE)
    if (length(unique(v)) <= 2) {
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rbinom(sum(ind), 1, 0.5)
    } else {
      sd_syn <- get_sd_stdnormal(v)
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rnorm(sum(ind), mean = median(v), sd = sd_syn)
    }
    u
  })
}

N_rep <- 100
M <- 1000
n_seq <- c(100)

NUM_TUNING <- 100

h0_indication <- "const"
rho <- 0.5




betahat_mixcat_function <- function(tau_this, X, Y, status, my_syndata) {
  M <- nrow(my_syndata$X.syn)
  X_merge <- rbind(X, my_syndata$X.syn)
  Y_merge <- c(Y, my_syndata$T.syn)
  status_merge <- c(status, my_syndata$status.syn)
  wt <- c(rep(1, nrow(X)), rep(tau_this / M, M))

  have_error22 <- FALSE
  tryCatch(
    {
      betahat_mix <- as.numeric(coef(coxph(Surv(Y_merge, status_merge) ~ X_merge, weights = wt, ties = "breslow")))
    },
    error = function(e) {
      have_error22 <<- TRUE
    },
    warning = function(w) {
      have_error22 <<- TRUE
    }
  )

  if (have_error22 == TRUE) {
    return(rep(100, ncol(X)))
  }
  return(betahat_mix)
}


n <- 100
for (p in c(20, 40, 60)) {
  for (censor_rate in c(0.1, 0.2, 0.4)) {
    betahat_mixcat_p_store <- matrix(0, nr = N_rep, nc = p)

    FILENAME <- paste0("p", p, "censor", censor_rate * 10, "n", n)
    print(FILENAME)
    file_path <- paste0(
      "~/Desktop/OneDrive - National University of Singapore/survivial/Aug7-2023/simulationgood_prediction_synskew_correct/sequce_weibull/cen",
      censor_rate * 10, "/p", p
    )
    setwd(file_path)


    rep_f <- function(rep_index) {
      set.seed(50 + rep_index)
      print(paste0(n, "rep:", rep_index))
      X_obs <- matrix(rnorm(n * p), n, p)
      # X_obs<-scale(X_obs)
      X_obs[, 1] <- rbinom(n, 1, 0.9)

      X_obs[, 3] <- rchisq(n, 1)
      X_obs[, 2] <- rchisq(n, 4)

      X <- X_obs
      beta <- c(4, -4, 3, -3, 1, -1, 1, -1, rep(1, p - 8)) / sqrt(p)

      Xbeta <- as.numeric(X %*% beta)
      TT <- sapply(1:n, function(i) {
        rexp(1, 0.5 * exp(Xbeta[i]))
      })
      if (censor_rate == 0.1) {
        pre_value <- quantile(TT, 0.95)

        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        target_censor_rate <- censor_rate
        while (abs(target_censor_rate - censor_rates) > 0.02) {
          C <- runif(n, 0, pre_value)
          censor_rates <- sum(C < TT) / n
          # print(censor_rates)
        }
      }
      if (censor_rate == 0.2) {
        pre_value <- quantile(TT, 0.90)

        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        target_censor_rate <- censor_rate
        while (abs(target_censor_rate - censor_rates) > 0.02) {
          C <- runif(n, 0, pre_value)
          censor_rates <- sum(C < TT) / n
          # print(censor_rates)
        }
      }
      if (censor_rate == 0.4) {
        pre_value <- quantile(TT, 0.70)

        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        target_censor_rate <- censor_rate
        while (abs(target_censor_rate - censor_rates) > 0.02) {
          C <- runif(n, 0, pre_value)
          censor_rates <- sum(C < TT) / n
          # print(censor_rates)
        }
      }
      Y_status <- t(sapply(1:n, function(aa) {
        c(min(C[aa], TT[aa]), (C[aa] > TT[aa]))
      }))
      Y <- Y_status[, 1]
      status <- Y_status[, 2]



      X.syn <- generate_synthetic_X(X, M)
      rate_exp <- sum(status) / mean(Y) / n
      T.syn <- rexp(M, rate = rate_exp)
      increase_order <- order(T.syn)
      X.syn <- X.syn[increase_order, ]
      T.syn <- T.syn[increase_order]
      status.syn <- rep(1, M)
      syntheticObj <- list()
      syntheticObj$M <- M
      syntheticObj$xstar <- X.syn
      syntheticObj$ystar <- T.syn
      syntheticObj$h_star <- rate_exp
      my_syndata <- list(T.syn = syntheticObj$ystar, X.syn = syntheticObj$xstar, status.syn = rep(1, M), h_star = rate_exp)

      betahat_tau_p_mix <- betahat_mixcat_function(p, X, Y, status, my_syndata)
      return(betahat_tau_p_mix)
    }
    for (rep_index in 1:N_rep) {
      betahat_mixcat_p_store[rep_index, ] <- rep_f(rep_index)
    }
    save(betahat_mixcat_p_store, file = paste0(FILENAME, "betahat_taumixp"))
  }
}
