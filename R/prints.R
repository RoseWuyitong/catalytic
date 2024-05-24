#' Title
#'
#' @param ...
#' @param data_
#' @param show_data
#' @param detail
#'
#' @return
#' @export
#'
#' @examples
print.data_ <- function(data_,
                        show_data = TRUE,
                        detail = TRUE,
                        ...) {
  cat(data_$function_name)
  if (data_$function_name == "data_initialization") {
    cat(
      "\n family:            ", family_plus_link(data_$family)
    )
  }
  cat("\n formula:           ", formula_string(data_$formula))
  cat("\n size:              ", data_$size)

  cat("\n------\n")
  cat("Observation Data Information:")
  cat("\n size:              ", data_$obs_size)
  if (show_data) {
    cat("\n head(data) :              \n")
    cat(" [only show the head of first 10 columns] \n")
    if (ncol(data_$obs_data) > 10) {
      print(head(data_$obs_data)[, 1:10])
    } else {
      print(head(data_$obs_data))
    }
  }
  cat("\n------\n")
  cat("Synthetic Data Information:")
  cat("\n size:              ", data_$syn_size)
  if (show_data) {
    cat("\n head(data):              \n")
    cat(" [only show the head of first 10 columns] \n")
    if (ncol(data_$syn_data) > 10) {
      print(head(data_$syn_data)[, 1:10])
    } else {
      print(head(data_$syn_data))
    }
  }

  if (data_$function_name == "data_initialization") {
    cat("\n synthetic data generation process:")
    cat("\n [only show the first 10 columns] \n")
    print(
      format(
        if (nrow(data_$syn_x_generation_steps) > 10) {
          data_$syn_x_generation_steps[1:10, ]
        } else {
          data_$syn_x_generation_steps
        },
        justify = "left",
        width = max(nchar(as.character(unlist(data_$syn_x_generation_steps)))) + 1
      ),
      quote = FALSE
    )
  }
  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.data_\n")
  }

  invisible(data_)
}

#' Title
#'
#' @param ...
#' @param cat_glm
#' @param digit
#' @param detail
#'
#' @return
#' @export
#'
#' @examples
print.none_bayes <- function(cat_glm, digit = 3, detail = TRUE, ...) {
  cat(cat_glm$function_name)
  cat("\n family:            ", family_plus_link(cat_glm$family))
  cat("\n formula:           ", formula_string(cat_glm$formula))
  cat("\n size:              ", cat_glm$data_$size)
  cat("\n tau:               ", cat_glm$tau)

  cat("\n------\n")
  cat("coefficients' information:\n")
  print(
    format(round(summary(cat_glm$model)$coefficients, digit),
      justify = "left"
    ),
    quote = FALSE
  )

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.none_bayes\n")
  }

  invisible(cat_glm)
}


#' Title
#'
#' @param ...
#' @param cat_glm
#'
#' @return
#' @export
#'
#' @examples
print.bayes <- function(cat_glm, ...) {
  cat(cat_glm$function_name)
  cat("\n family:            ", family_plus_link(cat_glm$family))
  cat("\n formula:           ", formula_string(cat_glm$formula))
  cat("\n size:              ", cat_glm$data_$size)
  cat("\n tau:               ", cat_glm$tau)

  cat("\n------\n")
  cat("coefficients' information:\n")
  print(cat_glm$sample_model)

  invisible(cat_glm)
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
print.bayes_gibbs <- function(cat_gibbs, digit = 2, detail = TRUE, ...) {
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
  cat("\n size:              ", cat_lmm$data_$size)
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
  cat("\n size:              ", cat_lmm_$data_$size)
  cat("\n taus:              ", cat_lmm_$taus)
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


# TODO:: combined

#' Title
#'
#' @param digit
#' @param detail
#' @param ...
#' @param cat_cox_data_
#'
#' @return
#' @export
#'
#' @examples
print.cat_cox <- function(cat_cox_data_,
                          digit = 3,
                          detail = TRUE, ...) {
  cat(cat_cox_data_$function_name)
  cat("\n size:              ", cat_cox_data_$data_$size)
  cat("\n method:            ", cat_cox_data_$method)
  cat("\n tau:               ", cat_cox_data_$tau)
  if (cat_cox_data_$method == "CRE") {
    cat(
      "\n iter:              ",
      paste(
        cat_cox_data_$iter,
        "/",
        cat_cox_data_$max_iter
      )
    )
  }

  cat("\n------\n")


  cat("coefficients' information:\n")
  if (cat_cox_data_$method == "CRE") {
    cat("\n initial value:\n")
    print(round(cat_cox_data_$iter_collect[1, ], digit))


    cat("\n final value:\n")
    print(round(cat_cox_data_$iter_collect[nrow(cat_cox_data_$iter_collect), ], digit))
  } else {
    cat(round(cat_cox_data_$coefficients, digit))
  }

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_cox\n")
  }


  invisible(cat_cox_data_)
}
