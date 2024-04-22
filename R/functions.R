get_log_density <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          Y,
                          model = NULL,
                          weights = 1,
                          ...,
                          coefficients = NULL) {
      return(log_density_binomial(Y=as.matrix(Y),
                                  mean=multiply_two_matrices(X,
                                                      if (is.null(model)
                                                          ) coefficients else model$coefficients),
                           weights=weights))
    },
    "gaussian" = function(model = NULL,
                          weights = 1,
                          ...,
                          residuals = NULL) {
      return(log_density_gaussian(residuals= if (is.null(model)
                                          ) residuals else model$residuals,
                           weights=weights))
    }
  )
}


get_diag_approx_cov <- function(X, model) {
  return(diag(chol2inv(chol(second_derivative(X=as.matrix(X),
                                              fitted_values=model$fitted.values,
                                              prior_weights=model$prior.weights)))))
}


get_log_density_grad <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          Y,
                          model = NULL,
                          weights=1,
                          ...,
                          coefficients = NULL) {
      return(
        log_density_grad_binomial(X=as.matrix(X),
                                  Y=as.matrix(Y),
                                  mean=get_mean("binomial")(
                                    X= X,
                                    coefficients = if (is.null(model)
                                        ) coefficients else model$coefficients),
                                    weights=weights)
        )}
  )
}

get_mean <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          model=NULL,
                          ...,
                          coefficients=NULL) {
      return(
        mean_of_binomial(X=as.matrix(X),
                          coefficients=if (is.null(model)
                                              ) coefficients else model$coefficients)
      )
    },
    "gaussian" = function(X,
                          model=NULL,
                          ...,
                          coefficients=NULL) {
      return(
        mean_of_gaussian(X=as.matrix(X),
                          coefficients=if (is.null(model)
        ) coefficients else model$coefficients)
      )
    }
  )
}


get_sim_y <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          model = NULL,
                          ...,
                          coefficients=NULL,
                          mean = NULL,
                          sd = NULL) {
      get_data(
        "binomial"
        )(
        n = nrow(X),
        mean = if (is.null(mean)) {
          get_mean(
            "binomial"
          )(
            X=X,
            model=model
          )
        } else {
          mean
        })
    },
    "gaussian" = function(X,
                          model = NULL,
                          ...,
                          mean = NULL,
                          sd = NULL) {
      get_data(
        "gaussian"
      )(
        n = nrow(X),
        mean = if (is.null(mean)) {
          get_mean(
            "gaussian"
          )(
            X = X,
            coefficients = model$coefficients
          )
        } else {
          mean
        },
        sd = if (is.null(sd)) sqrt(sum(model$residuals^2) / model$df.residual) else sd
      )
    }
  )
}

lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, coefficients, fam, ...) 2 * get_mean(fam)(X, coefficients),
    "classification_error" = function(X, coefficients, fam, ...) 2 * get_mean(fam)(X, coefficients),
    "logistic_deviance" = function(X, coefficients, fam, ...) {
      .not_same_stop(fam,
                     "binomial",
                     "only binomial family can use logistic_deviance")
      return(multiply_two_matrices(X, coefficients))
    }
  )
}

get_iqr_diff <- function(x, invl = 0.5) {
  right_quantile <- (2 - invl) / 2
  iqr_value <- 0
  while (iqr_value == 0 & right_quantile < 1) {
    iqr_value <- diff(
      quantile(
        as.numeric(x),
      c(1 - right_quantile, right_quantile),
      na.rm = FALSE,
      names = FALSE,
      type = 7
    ))
    right_quantile <- right_quantile + 0.05
  }

  return(iqr_value)
}

get_data <- function(fam) {
  switch(fam,
    "binomial" = function(n = 10,
                          mean = 0,
                          ...) {
      stats::rbinom(
        n,
        1,
        mean
      )
    },
    "gaussian" = function(n = 10,
                          mean = 0,
                          sd = 1,
                          ...) {
      stats::rnorm(
        n,
        mean,
        sd
      )
    },
    "uniform" = function(n = 10,
                         mean = 0,
                         ...) {
      rep(mean,
          n)
    }
  )
}

get_discrepancy <- function(func) {
  switch(func,
    "logarithmic_error" = function(Y,
                                   est_Y=NULL,
                             ...,
                             fam=NULL,
                             est_model=NULL,
                             est_coefficients=NULL) {
      # cap the value between 0.0001 and 0.9999 for better calculation process
      return(
        discrepancy_of_logarithmic_error(
          Y=pmax(0.0001, pmin(0.9999, Y)),
          est_Y = pmax(0.0001, pmin(0.9999,
                            if (is.null(est_Y)) get_mean(
                              fam
                            )(
                              X = X,
                              coefficients = if (is.null(est_model)) est_coefficients else est_model$coefficients
                            )  else
                              est_Y))
        )
      )
    },
    "square_error" = function(Y,
                              est_Y=NULL,

                              ...,
                              X=NULL,
                              fam=NULL,
                              est_model=NULL,
                              est_coefficients=NULL) {
      return(
        discrepancy_of_square_error(
          Y=Y,
          est_Y=if (is.null(est_Y)) get_mean(
            fam
          )(
            X = X,
            coefficients = if (is.null(est_model)) est_coefficients else est_model$coefficients
          ) else est_Y
        )
      )
    },
    "classification_error" = function(Y,
                                      est_Y=NULL,


                                      ...,
                                      X=NULL,
                                      fam=NULL,
                                      est_model=NULL,
                                      est_coefficients=NULL) {

      return(
        discrepancy_of_classification_error(
          Y=Y,
          est_Y = if (is.null(est_Y)) get_mean(
            fam
          )(
            X = X,
            coefficients = if (is.null(est_model)) est_coefficients else est_model$coefficients
          ) else est_Y
        )
        )
    },
    "logistic_deviance" = function(X,
                                   Y,
                                   fam,
                                   model=NULL,
                                   coefficients=NULL,
                                   ...) {
      - get_log_density(fam)(
          X = X,
          Y = Y,
        model = model,
        coefficients=coefficients)
    }
  )
}



stan_model_selection <- memoise::memoise(function(fam,
                                                  joint,
                                                  gibbs,
                                                  use_theta) {
  return(rstan::stan_model(
    file = paste0(
      "~/catalytic/R/rstan_script/",
      fam,
      if (joint & isFALSE(gibbs)) "_joint",
      "_bayesian",
      if (fam == "binomial" & joint & use_theta & isFALSE(gibbs)) "_theta",
      if (fam == "binomial" & joint & isFALSE(use_theta) & isFALSE(gibbs)) "_tau",
      ".stan"
    )
  ))
})

standardize_data <- function(data, if_missing = NULL) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  drop_redundant_dims(data)

  for (i in 1:ncol(data)) {
    if (!(NA %in% suppressWarnings(as.numeric(data[, i])))) {
      data[, i] <- as.numeric(data[, i])
    } else if (!(NA %in% suppressWarnings(as.logical(data[, i])))) {
      data[, i] <- as.logical(data[, i])
    } else {
      data[, i] <- as.factor(data[, i])
    }
  }

  return(data)
}


drop_redundant_dims <- function(data) {
  drop_dim <- sapply(data, function(v) is.matrix(v) && NCOL(v) == 1)
  data[, drop_dim] <- lapply(data[, drop_dim, drop = FALSE], drop)
  return(data)
}

validate_glm_formula <- function(f) {
  if (!any(grepl("~", f))) {
    stop("please enter ~ for formula",
      call. = FALSE
    )
  }
}



validate_family <- function(f) {
  if (is.character(f)) {
    f <- get(f, mode = "function", envir = parent.frame(2))
  }
  if (is.function(f)) {
    f <- f()
  }
  if (!is(f, "family")) {
    stop("'family' must be a family.", call. = FALSE)
  }

  return(f)
}

get_family <- function(family) {
  fam <- try(family$family, silent = TRUE)
  if (inherits(fam, "try-error")) {
    fam <- family()$family
  }
  return(fam)
}

family_plus_link <- function(fam) {
  if (is.character(fam)) {
    stopifnot(identical(fam, x$method))
    fam <- paste0("ordered [", fam, "]")
  } else {
    fam <- paste0(fam$family, " [", fam$link, "]")
  }
  return(fam)
}

formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) {
    "--MARK--"
  } else {
    " "
  }
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent) {
    return(char)
  }
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

rhs <- function(formula, with_ = FALSE) {
  # Split the formula string by "~" and extract the first element
  if (with_) {
    return(as.formula(paste("~", deparse(rlang::f_rhs(formula)), collapse = " ")))
  }
  return(rlang::f_rhs(formula))
}

lhs <- function(formula) {
  return(rlang::f_lhs(formula))
}

H_matrix <- function(X, X_star, tau, M, n) {
  X <- as.matrix(X)
  X_star <- as.matrix(X_star)
  return(X %*% (
    chol2inv(chol((t(X) %*% X + (tau / M) * (t(X_star) %*% X_star))))) %*% (
    t(X) + tau / (M * n) * t(X_star) %*% rep(1, M) %*% t(rep(1, n))
  ))
}
