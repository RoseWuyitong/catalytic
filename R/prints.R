#' Title
#'
#' @param ...
#' @param cat_init
#' @param show_data
#' @param detail
#'
#' @return
#' @export
#'
#' @examples
print.cat_initialization <- function(cat_init,
                                     show_data = TRUE,
                                     detail = TRUE,
                                     ...) {
  cat(cat_init$function_name)

  cat("\n formula:           ", get_formula_string(cat_init$formula))
  if (cat_init$function_name == "data_initilization_glm") {
    cat("\n family:            ", get_glm_family(cat_init$family, with_link = TRUE))
  }
  cat("\n size:              ", extract_size(cat_init))

  cat("\n------\n")
  cat("Observation Data Information:")
  cat("\n size:              ", cat_init$obs_size)
  if (show_data) {
    cat("\n head(data) :              \n")
    cat(" [only show the head of first 10 columns] \n\n")
    if (ncol(cat_init$obs_data) > 10) {
      print(head(cat_init$obs_data)[, 1:10])
    } else {
      print(head(cat_init$obs_data))
    }
  }
  cat("\n------\n")
  cat("Synthetic Data Information:")
  cat("\n size:              ", cat_init$syn_size)
  if (show_data) {
    cat("\n head(data):              \n")
    cat(" [only show the head of first 10 columns] \n\n")
    if (ncol(cat_init$syn_data) > 10) {
      print(head(cat_init$syn_data)[, 1:10])
    } else {
      print(head(cat_init$syn_data))
    }
  }

  if (!is.null(cat_init$syn_x_resample_inform)) {
    cat("\n data generation process:")
    cat("\n [only show the first 10 columns] \n\n")
    print(
      format(
        if (nrow(cat_init$syn_x_resample_inform) > 10) {
          cat_init$syn_x_resample_inform[1:10, ]
        } else {
          cat_init$syn_x_resample_inform
        },
        justify = "left",
        width = max(nchar(as.character(unlist(cat_init$syn_x_resample_inform)))) + 1
      ),
      quote = FALSE
    )
  }

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_initialization\n")
  }

  invisible(cat_init)
}


#' Title
#'
#' @param cat_glm_
#' @param digit
#' @param detail
#'
#' @return
#' @export
#'
#' @examples
print.cat_glm <- function(cat_glm_,
                          digit = 3,
                          detail = TRUE) {
  cat(cat_glm_$function_name)
  cat("\n family:            ", get_glm_family(cat_glm_$family, with_link = TRUE))
  cat("\n formula:           ", get_formula_string(cat_glm_$formula))
  cat("\n size:              ", extract_size(cat_glm_$cat_init))
  cat("\n tau:               ", cat_glm_$tau)

  cat("\n------\n")
  cat("coefficients' information:\n\n")

  print(extract_coefs(cat_glm_, digit), quote = FALSE)

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_glm\n")
  }

  invisible(cat_glm_)
}

#' Title
#'
#' @param cat_glm_tune_
#' @param digit
#' @param detail
#'
#' @return
#' @export
#'
#' @examples
print.cat_glm_tune <- function(cat_glm_tune_,
                               digit = 3,
                               detail = TRUE) {
  cat(cat_glm_tune_$function_name)
  cat("\n family:                    ", get_glm_family(cat_glm_tune_$family, with_link = TRUE))
  cat("\n formula:                   ", get_formula_string(cat_glm_tune_$formula))
  cat("\n size:                      ", extract_size(cat_glm_tune_$cat_init))
  cat("\n risk estimate method:      ", cat_glm_tune_$risk_estimate_method)
  cat("\n discrepancy function:      ", cat_glm_tune_$discrepancy_method)

  if (cat_glm_tune_$risk_estimate_method == "parametric_bootstrap") {
    cat("\n bootstrap iteration times: ", cat_glm_tune_$bootstrap_iteration_times)
  } else if (cat_glm_tune_$risk_estimate_method == "cross_validation") {
    cat("\n cross validation folds:    ", cat_glm_tune_$cv_fold_num)
  }
  cat("\n tau_seq:                   ", extract_tau_seq(cat_glm_tune_, digit))
  cat("\n")
  cat("\n optimal tau:               ", round(cat_glm_tune_$tau, digit))
  cat("\n minimun risk estimate:     ", round(cat_glm_tune_$min_risk_estimate, digit))

  cat("\n------\n")
  cat("coefficients' information:\n\n")

  print(extract_coefs(cat_glm_tune_, digit), quote = FALSE)

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_glm_tune\n")
  }

  invisible(cat_glm_tune_)
}

#' Title
#'
#' @param digit
#' @param detail
#' @param ...
#' @param bayes
#'
#' @return
#' @export
#'
#' @examples
print.cat_glm_bayes <- function(cat_glm_bayes_,
                                digit = 3,
                                detail = TRUE, ...) {
  cat(cat_glm_bayes_$function_name)
  cat("\n family:            ", get_glm_family(cat_glm_bayes_$family))
  cat("\n formula:           ", get_formula_string(cat_glm_bayes_$formula))
  cat("\n size:              ", extract_size(cat_glm_bayes_$cat_init))
  if (cat_glm_bayes_$function_name == "cat_glm_bayes") {
    cat("\n tau:               ", cat_glm_bayes_$tau)
  }
  cat("\n stan chain:        ", cat_glm_bayes_$stan_chain)
  cat("\n stan iter:         ", cat_glm_bayes_$stan_iter)
  cat("\n stan warmup:       ", cat_glm_bayes_$stan_warmup)
  cat("\n stan algorithm:    ", cat_glm_bayes_$stan_algorithm)

  cat("\n------\n")
  cat("coefficients' information:\n\n")
  print(extract_stan_summary(cat_glm_bayes_, digit))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_glm_bayes\n")
  }
  invisible(cat_glm_bayes_)
}


#' Title
#'
#' @param digit
#' @param detail
#' @param ...
#' @param cat_glm_bayes_joint_gibbs_
#'
#' @return
#' @export
#'
#' @examples
print.cat_glm_bayes_joint_gibbs <- function(cat_glm_bayes_joint_gibbs_, digit = 2, detail = TRUE, ...) {
  cat(cat_glm_bayes_joint_gibbs_$function_name)
  cat("\n family:            ", get_glm_family(cat_glm_bayes_joint_gibbs_$family))
  cat("\n formula:           ", get_formula_string(cat_glm_bayes_joint_gibbs_$formula))
  cat("\n size:              ", extract_size(cat_glm_bayes_joint_gibbs_$cat_init))

  cat("\n------\n")
  cat("coefficients' information:\n\n")
  cat(paste0(
    "Inference for Gibb's sampling.\n",
    "gibbs_iter=", cat_glm_bayes_joint_gibbs_$gibbs_iter,
    "; gibbs_warmup=", cat_glm_bayes_joint_gibbs_$gibbs_warmup, "; ",
    "gibb's sampling post-warmup draws=",
    cat_glm_bayes_joint_gibbs_$gibbs_iter - cat_glm_bayes_joint_gibbs_$gibbs_warmup,
    ".\n\n"
  ))

  print(round(cat_glm_bayes_joint_gibbs_$inform_df, digit))

  cat(paste0(
    "\nSamples were drawn using Gibb's Sampling at ",
    formatted_time <- format(cat_glm_bayes_joint_gibbs_$sys_time, "%a %b %d %H:%M:%S %Y"),
    ".\nFor each parameter, n_eff is a crude measure of effective sample size."
  ))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_glm_bayes_joint_gibbs\n")
  }

  invisible(cat_glm_bayes_joint_gibbs_)
}


#' Title
#'
#' @param digit
#' @param detail
#' @param ...
#' @param cat_cox_model
#'
#' @return
#' @export
#'
#' @examples
print.cat_cox <- function(cat_cox_,
                          digit = 3,
                          detail = TRUE, ...) {
  cat(cat_cox_$function_name)
  cat("\n formula:           ", get_formula_string(cat_cox_$formula))
  cat("\n size:              ", extract_size(cat_cox_$cat_init))
  cat("\n tau:               ", round(cat_cox_$tau, 3))
  cat("\n method:            ", cat_cox_$method)
  if (cat_cox_$method == "CRE") {
    cat(
      "\n iter:              ",
      paste(
        cat_cox_$iter,
        "/",
        cat_cox_$max_iter
      )
    )
  }

  cat("\n------\n")
  cat("coefficients' information:\n\n")

  if (cat_cox_$method == "CRE") {
    last_iter_coef <- cat_cox_$iter_collect[nrow(cat_cox_$iter_collect), ]
    names(last_iter_coef) <- colnames(cat_cox_$cat_init$adj_obs_x)
    print(round(last_iter_coef, digit))
  } else {
    print(extract_coefs(cat_cox_, digit, FALSE))
  }


  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_cox\n")
  }

  invisible(cat_cox_)
}

#' Title
#'
#' @param cat_cox_model
#' @param digit
#' @param detail
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_cox_tune <- function(cat_cox_tune_,
                               digit = 3,
                               detail = TRUE, ...) {
  cat(cat_cox_tune_$function_name)
  cat("\n formula:              ", get_formula_string(cat_cox_tune_$formula))
  cat("\n size:                 ", cat_cox_tune_$cat_init$size)
  cat("\n method:               ", cat_cox_tune_$method)
  cat("\n tau_seq:              ", extract_tau_seq(cat_cox_tune_, digit))
  cat("\n")
  cat("\n optimal tau:          ", round(cat_cox_tune_$tau, digit))
  cat("\n optimal likelihood:   ", round(cat_cox_tune_$max_likelihood, digit))

  cat("\n------\n")
  cat("coefficients' information:\n")
  print(extract_coefs(cat_cox_tune_, digit, FALSE))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_cox_tune\n")
  }

  invisible(cat_cox_tune_)
}

#' Title
#'
#' @param cat_cox_bayes_
#' @param digit
#' @param detail
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_cox_bayes <- function(cat_cox_bayes_,
                                digit = 3,
                                detail = TRUE, ...) {
  cat(cat_cox_bayes_$function_name)

  cat("\n formula:            ", get_formula_string(cat_cox_bayes_$formula))
  cat("\n size:               ", extract_size(cat_cox_bayes_$cat_init))
  if (cat_cox_bayes_$function_name == "cat_cox_bayes") {
    cat("\n tau_seq:            ", extract_tau_seq(cat_cox_bayes_, digit))
  }
  cat("\n")
  cat("\n stan chain:         ", cat_cox_bayes_$stan_chain)
  cat("\n stan iter:          ", cat_cox_bayes_$stan_iter)
  cat("\n stan warmup:        ", cat_cox_bayes_$stan_warmup)

  cat("\n------\n")
  cat("coefficients' information:\n\n")
  print(extract_stan_summary(cat_cox_bayes_, digit, with_intercept = FALSE))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_cox_bayes\n")
  }

  invisible(cat_cox_bayes_)
}


#' Title
#'
#' @param cat_lmm
#' @param digit
#' @param detail
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_lmm <- function(cat_lmm,
                          digit = 3,
                          detail = TRUE, ...) {
  cat(cat_lmm$function_name)
  cat("\n size:              ", cat_lmm$cat_init$size)
  cat("\n tau:               ", cat_lmm$tau)
  cat("\n maximum iteration: ", cat_lmm$max_iter)
  cat("\n total iteration:   ", nrow(cat_lmm$iter_collect) - 1)
  cat("\n tolerance:         ", cat_lmm$tol)
  cat("\n optimize domain:   ", cat_lmm$optimize_domain)
  cat("")

  cat("\n------\n")
  cat("coefficients' information:\n")

  cat("\n initial value:\n")
  print(round(cat_lmm$iter_collect[1, ], digit))


  cat("\n final value:\n")
  print(round(cat_lmm$iter_collect[nrow(cat_lmm$iter_collect), ], digit))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_lmm\n")
  }

  invisible(cat_lmm)
}


#' Title
#'
#' @param digit
#' @param detail
#' @param ...
#' @param cat_lmm_
#'
#' @return
#' @export
#'
#' @examples
print.cat_lmm_tune <- function(cat_lmm_, digit = 3, detail = TRUE, ...) {
  cat(cat_lmm_$function_name)
  cat("\n size:              ", cat_lmm_$cat_init$size)
  cat("\n tau_seq:           ", cat_lmm_$tau_seq)
  cat("\n optimal tau:       ", cat_lmm_$tau)
  cat("\n minimun error:     ", cat_lmm_$min_risk_estimate)

  cat("\n------\n")
  cat("optimal model:\n")

  print(cat_lmm_$model)

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_lmm_tune\n")
  }

  invisible(cat_lmm_)
}

extract_size <- function(cat_init) {
  return(paste0(
    cat_init$obs_size,
    " (Observation) + ",
    cat_init$syn_size,
    " (Synthetic) = ",
    cat_init$size
  ))
}

extract_coefs <- function(cat_, digit = 3, with_intercept = TRUE) {
  print_coef <- coef(cat_$model)

  if (with_intercept) {
    names(print_coef) <- c("(Intercept)", colnames(cat_$cat_init$adj_x))
  } else {
    names(print_coef) <- colnames(cat_$cat_init$adj_x)
  }
  return(round(print_coef, digit))
}

extract_tau_seq <- function(cat_, digit = 3) {
  tau_seq <- round(cat_$tau_seq, digit)
  return(if (length(tau_seq) > 10) {
    paste(
      paste(tau_seq[1:3], collapse = ", "),
      "...",
      paste(
        tau_seq[(length(tau_seq) - 2):length(tau_seq)],
        collapse = ", "
      )
    )
  } else {
    paste(tau_seq, collapse = ", ")
  })
}

extract_stan_summary <- function(cat_, digit = 3, with_intercept = TRUE) {
  print_stan_summary <- rstan::summary(cat_$stan_sample_model)$summary
  if (with_intercept) {
    rownames(print_stan_summary)[1:(ncol(cat_$cat_init$adj_x) + 1)] <- c(
      "(Intercept)",
      colnames(cat_$cat_init$adj_x)
    )
  } else {
    rownames(print_stan_summary)[1:(ncol(cat_$cat_init$adj_x))] <- colnames(cat_$cat_init$adj_x)
  }
  # TODO move to cat function inside

  return(round(print_stan_summary, digit))
}
