

#' Title
#'
#' @param data_
#' @param tau
#' @param method
#' @param ...
#' @param initial_values
#' @param hazard_function
#' @param tol
#' @param max_iter
#' @param entry
#' @param eta0 Hyperparameters for gamma process prior construction
#' @param kappa0 Hyperparameters for gamma process prior construction
#'
#' @return
#' @export
#'
#' @examples
cat_cox_bayes <- function(data_,
                          tau = NULL,
                          method = c("CRE", "WME"),
                          ...,
                          initial_values = NULL,
                          hazard_function = NULL,
                          tol = 1e-5,
                          max_iter = 25,
                          entry = NULL,
                          eta0 = 1,
                          kappa0 = 1) {
  .not_true_stop(
    is.data_cox(data_),
    paste(
      "Please enter `data_` that generated from",
      "catalytic::data_initialization_cox."
    )
  )

  if (is.null(tau)) tau <- ncol(data_$obs_data)
  if (is.null(hazard_function)) hazard_function <- function(t) sum(data_$obs_event) / mean(data_$obs_time) / data_$obs_size

  c0 <- 2
  s <- c(sort(Y[status == 1]), 2 * max(Y) - max(Y[-which(Y == max(Y))]))
  s <- unique(s)
  J <- length(s)

  # Calculate intervals
  intv <- setting.interval(Y, status, s, J)

  # Compute H.star and alpha0
  surv_weibull_model <- do.call(survival::survreg,
                                list(formula=survival::Surv(data_$obs_time, data_$obs_event) ~ 1,
                                     dist = "weibull"))
  weibull_shape <- 1 / surv_weibull_model$scale
  weibull_shape <- exp(coef(surv_weibull_model))
  H.star <- alpha0 <- numeric(J)
  for (j in 1:J) {
    H.star[j] <- (s[j] / weibull_shape)^weibull_shape # eta0 * s[j] ^ kappa0
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
