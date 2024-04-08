#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.initialization <- function(data_, data_head = FALSE, detail = TRUE, ...) {
  cat("initialization")
  cat("\n family:            ", family_plus_link(data_$family))
  cat("\n formula:           ", formula_string(data_$formula))
  cat("\n size:              ", data_$size)

  cat("\n------\n")
  cat("Observation Data Information")
  cat("\n size:              ", data_$obs_size)
  if (data_head) {
    cat("\n head(data):              \n")
    print(head(data_$obs_data))
  }
  cat("\n------\n")
  cat("Synthetic Data Information")
  cat("\n size:              ", data_$syn_size)
  cat("\n generation process:\n")
  print(
    format(data_$syn_x_generation_steps,
      justify = "left",
      width = max(nchar(as.character(unlist(data_$syn_x_generation_steps)))) + 1
    ),
    quote = FALSE
  )

  if (data_head) {
    cat("\n head(data):             \n")
    print(head(data_$syn_data))
  }

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.initialization\n")
  }

  invisible(data_)
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_tau_glm <- function(cat_glm, digit = 3, detail = TRUE, ...) {
  cat("cat_tau_glm")
  cat("\n family:            ", family_plus_link(cat_glm$family))
  cat("\n formula:           ", formula_string(cat_glm$formula))
  cat("\n size:              ", cat_glm$data_$size)
  cat("\n tau:               ", cat_glm$tau)

  cat("\n------\n")
  print(
    format(round(summary(cat_glm$model)$coefficients, digit),
      justify = "left"
    ),
    quote = FALSE
  )
  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_tau_glm\n")
  }

  invisible(cat_glm)
}


#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_freq_glm <- function(x, ...) {
  cat("\nFormula for frequntist prior:  ")
  print(x$formula)
  cat(
    "\nRisk estimate method: ", x$risk_estimate_method,
    "\n\nOptimal tau: ", x$tau,
    "\n\nOptimal estimated beta: ", round(x$beta, 3)
  )
  invisible(x)
}


#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_bayesian_glm <- function(x, ...) {
  cat(
    "\nMCMC algorithm: ", x$algorithm,
    "\n\nJoint tau: ", x$joint,
    "\n\nMean of estimated argments: ", round(x$beta, 3)
  )
  invisible(x)
}


#' Title
#'
#' @param cat_gibbs
#' @param digit
#' @param detail
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_gibbs_glm <- function(cat_gibbs, digit = 2, detail = TRUE, ...) {
  cat(paste0(
    "Inference for Gibb's sampling.\n",
    "gibbs_iter=", cat_gibbs$gibbs_iter,
    "; gibbs_warmup=", cat_gibbs$gibbs_warmup, "; ",
    "gibb's sampling post-warmup draws=",
    cat_gibbs$gibbs_iter - cat_gibbs$gibbs_warmup,
    ".\n",
    "For each gibb's sampling iterration: ",
    "Inference for Stan model: anon_model.\n",
    "1 chains, each with stan_iter=", cat_gibbs$stan_iter,
    "; stan_warmup=", cat_gibbs$stan_warmup, "; ",
    "stan sampling post-warmup draws=",
    cat_gibbs$stan_iter - cat_gibbs$stan_warmup,
    ".\n\n"
  ))

  print(round(cat_gibbs$inform_df, digit))

  cat(paste0(
    "\nSamples were drawn using Gibb's Sampling at ",
    formatted_time <- format(Sys.time(), "%a %b %d %H:%M:%S %Y"),
    ".\nFor each parameter, n_eff is a crude measure of effective sample size."
  ))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_gibbs_glm\n")
  }

  invisible(cat_gibbs)
}
