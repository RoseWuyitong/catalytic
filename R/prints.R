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
print.data_initialization <- function(data_,
                                      show_data = TRUE,
                                      detail = TRUE, ...) {
  cat("initialization")
  cat("\n family:            ", family_plus_link(data_$family))
  cat("\n formula:           ", formula_string(data_$formula))
  cat("\n size:              ", data_$size)

  cat("\n------\n")
  cat("Observation Data Information:")
  cat("\n size:              ", data_$obs_size)
  if (show_data) {
    cat("\n head(data):              \n")
    print(head(data_$obs_data))
  }
  cat("\n------\n")
  cat("Synthetic Data Information:")
  cat("\n size:              ", data_$syn_size)
  if (show_data) {
    cat("\n head(data):              \n")
    print(head(data_$syn_data))
  }
  cat("\n synthetic data generation process:\n")
  print(
    format(data_$syn_x_generation_steps,
      justify = "left",
      width = max(nchar(as.character(unlist(data_$syn_x_generation_steps)))) + 1
    ),
    quote = FALSE
  )


  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.data_initialization\n")
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
