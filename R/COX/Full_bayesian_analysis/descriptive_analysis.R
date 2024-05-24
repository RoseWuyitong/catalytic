library(survival)
library(mvtnorm)
library(glmnet)
library(survival)
source("coxph_as.R")
library(emulator)
library(MASS)


get_sd_stdnormal <- function(ddd) {
  sigma_seq <- seq(0.1, 10, by = 0.02)
  target_index <- which.min(abs(pnorm(IQR(ddd) / 2, mean = 0, sd = sigma_seq) - 0.75))
  return(sigma_seq[target_index])
}

generate_synthetic_X_pbc <- function(X_obs, M) {
  # column 7 and 18 are edma, need to resample together
  # categorical_index=c(1,3,4,5,6,7,18)
  syntheticMatrix <- matrix(0, nc = ncol(X_obs), nr = M)
  ###### for categorical
  for (i in c(1, 3, 4, 5, 6)) {
    u <- sample(X_obs[, i], size = M, replace = TRUE)
    ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
    u[ind] <- rbinom(sum(ind), 1, 0.5)
    syntheticMatrix[, i] <- u
  }
  #### for categorical edma
  index_resample_edma <- sample(1:nrow(X_obs), size = M, replace = T)
  ind_edma <- sample(c(TRUE, FALSE), M, replace = TRUE)
  syntheticMatrix[, c(7, 18)] <- X_obs[index_resample_edma, c(7, 18)]
  # sample (1,0) (0,1) and (0,0) random
  replace_index <- which(ind_edma)
  for (row_i in replace_index) {
    replace_value <- NULL
    ru <- runif(1)
    if (ru < 1 / 3) {
      replace_value <- c(1, 0)
    } else {
      if (ru > 2 / 3) {
        replace_value <- c(0, 1)
      } else {
        replace_value <- c(0, 0)
      }
    }
    syntheticMatrix[row_i, c(7, 18)] <- replace_value
  }
  ###### for numerical variable
  for (i in c(2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)) {
    u <- sample(X_obs[, i], size = M, replace = TRUE)
    sd_syn <- get_sd_stdnormal(X_obs[, i])
    ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
    u[ind] <- rnorm(sum(ind), mean = median(X_obs[, i]), sd = sd_syn)
  }
  return(syntheticMatrix)
}

library("rstan")
setting.interval <- function(y, delta, s, J) {
  n <- length(y)

  smax <- max(s)

  case0 <- which(delta == 0)
  case1 <- which(delta == 1)

  case0yleq <- which(delta == 0 & y <= smax)
  case0ygeq <- which(delta == 0 & y > smax)
  case1yleq <- which(delta == 1 & y <= smax)
  case1ygeq <- which(delta == 1 & y > smax)


  ind.d <- ind.r <- matrix(0, n, J)

  for (i in case1yleq) {
    d.mat.ind <- min(which(s - y[i] >= 0))
    ind.d[i, d.mat.ind] <- 1
    ind.r[i, 1:d.mat.ind] <- 1
  }

  for (i in case0yleq) {
    cen.j <- min(which(s - y[i] >= 0))
    ind.r[i, 1:cen.j] <- 1
  }

  if (length(union(case1ygeq, case0ygeq)) > 0) {
    ind.r[union(case1ygeq, case0ygeq), ] <- 1
  }

  ind.r_d <- ind.r - ind.d


  d <- colSums(ind.d)

  list(
    ind.r = ind.r,
    ind.d = ind.d,
    d = d,
    ind.r_d = ind.r_d
  )
}

MODEL_cat <- stan_model("bayesian_cox_cat.stan", auto_write = rstan_options("auto_write" = TRUE))


Bayesian_cox_cat_stan <- function(MODEL_cat, X, Y, status, X.syn, T.syn, status.syn, tau = NULL, h_daga = NULL) {
  # Extract dimensions
  p <- ncol(X)
  n <- nrow(X)
  M <- nrow(X.syn)

  # Set default values
  if (is.null(tau)) {
    tau <- p
  }

  if (is.null(h_daga)) {
    h_daga <- function(t) {
      sum(status) / mean(Y) / n
    }
  }

  # Hyperparameters for gamma process prior construction
  eta0 <- 1
  kappa0 <- 1

  weibull_fit <- survreg(Surv(Y, status) ~ 1, dist = "weibull")
  shape_fit <- 1 / weibull_fit$scale
  scale_fit <- exp(coef(weibull_fit))
  c0 <- 2
  s <- c(sort(Y[status == 1]), 2 * max(Y) - max(Y[-which(Y == max(Y))]))
  s <- unique(s)
  J <- length(s)

  # Calculate intervals
  intv <- setting.interval(Y, status, s, J)

  # Compute H.star and alpha0
  H.star <- alpha0 <- numeric(J)
  for (j in 1:J) {
    H.star[j] <- (s[j] / scale_fit)^shape_fit # eta0 * s[j] ^ kappa0
    alpha0[j] <- c0 * H.star[j]
  }

  hPriorSh <- diff(c(0, alpha0))

  # Compute H_daga_Y_star
  H_daga_Y_star <- sapply(1:M, function(index) {
    integrate(Vectorize(h_daga), lower = 0, upper = T.syn[index])$value
  })

  # For grouped likelihood
  R_tilde_minus_D_tilde <- intv$ind.r_d
  D_tilde <- intv$ind.d

  # Stan model and sampling

  fit <- sampling(
    MODEL_cat,
    data = list(
      J = J,
      hPriorSh = hPriorSh,
      c0 = c0,
      eta0 = eta0,
      kappa0 = kappa0,
      n = n,
      p = p,
      X = X,
      M = M,
      X_star = X.syn,
      Y_star = T.syn,
      status_star = status.syn,
      H_daga_Y_star = H_daga_Y_star,
      tau_downweight = tau,
      R_tilde_minus_D_tilde = R_tilde_minus_D_tilde,
      D_tilde = D_tilde
    ),
    chains = 1
  )

  return(fit)
}




set.seed(123)
mydata <- pbc[complete.cases(pbc), ]
reshuffle <- sample(1:nrow(mydata), nrow(mydata))
mydata <- mydata[reshuffle, ]
X <- mydata[, 4:20]
X$trt <- ifelse(X$trt == 1, 1, 0)
X$sex <- ifelse(X$sex == "f", 1, 0)
X$edema2 <- ifelse(X$edema == 0.5, 1, 0)
X$edema <- ifelse(X$edema == 1, 1, 0)
# X$bili=log(X$bili)
# X$albumin=log(X$albumin)
# X$copper=log(X$copper)
Y <- mydata$time
Y <- Y / 365
status <- (mydata$status == 2) * 1

X <- as.matrix(X)
categorical_index <- c(1, 3, 4, 5, 6, 7, 18)
mean_seq <- apply(X[, -categorical_index], 2, mean)
sd_seq <- apply(X[, -categorical_index], 2, sd)
for (i in 1:nrow(X)) {
  X[i, -categorical_index] <- (X[i, -categorical_index] - mean_seq) / sd_seq
}
head(X)
summary(coxph(Surv(Y, status) ~ X))

M <- 1000
n <- nrow(X)
p <- ncol(X)


X.syn <- generate_synthetic_X_pbc(X, M)
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




fit_cat_p <- Bayesian_cox_cat_stan(MODEL_cat, X, Y, status, X.syn, T.syn, status.syn, tau = p)

fit_sampler <- extract(fit_cat_p)
dim(fit_sampler$beta)
hist(fit_sampler$beta[, 1], main = "Sampling distribution for Coefficient of DPCA", xlab = "")

# give posterior mean and credible interval for each beta with tau=p

post_mean_vector <- rep(0, p)
post_CrI_vector <- rep("ss", p)
for (index in 1:p) {
  print(paste0(
    "post mean: ", round(mean(fit_sampler$beta[, index]), 3), ", 95% CrI: [",
    round(quantile(fit_sampler$beta[, index], 0.025), 2), " , ",
    round(quantile(fit_sampler$beta[, index], 0.975), 2), "]"
  ))
  post_mean_vector[index] <- round(mean(fit_sampler$beta[, index]), 3)
  post_CrI_vector[index] <- paste0(
    "[", round(quantile(fit_sampler$beta[, index], 0.025), 2),
    ",", round(quantile(fit_sampler$beta[, index], 0.975), 2), "]"
  )
}
colnames(X)
Full_col_name <- c(
  "Treatment (DPCA)", "Age", "Gender", "Presence of ascites", "Hepatomegaly", "Blood vessel malformations",
  "Edema1 (despite diuretic therapy)", "Serum bilirunbin",
  "Serum cholesterol", "Serum albumin", "Urine copper", "Phosphotase", "Aspartate aminotransferase",
  "Triglycerides", "Platelet count", "Standardised blood clotting time",
  "Histologic stage of disease", "Edema2 (untreated or successfully treated)"
)
summary_dat <- data.frame(
  Variable = Full_col_name, Posterior_mean_coefficient = post_mean_vector,
  Credible_interval = post_CrI_vector
)
noquote(summary_dat)[c(1, 2, 3, 4, 5, 7, 18, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 6), ]
write.csv(noquote(summary_dat)[c(1, 2, 3, 4, 5, 7, 18, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 6), ], file = "taup_simple_analysis.csv")

# MPLE result


MPLE_dataf <- confint(coxph(Surv(Y, status) ~ X))
mean_MPLE <- round(as.numeric(coef(coxph(Surv(Y, status) ~ X))), 3)
CI_MPLE <- paste0(
  "[", round(MPLE_dataf[, 1], 2),
  ",", round(MPLE_dataf[, 2], 2), "]"
)

summary_dat_MPLE <- data.frame(
  Variable = Full_col_name, mean_coefficient = mean_MPLE,
  C_interval = CI_MPLE
)
noquote(summary_dat_MPLE)[c(1, 2, 3, 4, 5, 7, 18, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 6), ]
write.csv(noquote(summary_dat_MPLE)[c(1, 2, 3, 4, 5, 7, 18, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 6), ], file = "MPLE_simple_analysis.csv")
