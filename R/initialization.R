initialization <- function(formula,
                           family = gaussian(),
                           data,
                           weights = NULL,
                           subset = NULL,
                           na.action = NULL,
                           offset = NULL,
                           ...,
                           syn_size = NULL,
                           x_degree = NULL) {
  validate_glm_formula(formula)
  family <- validate_family(family)
  fam <- family$family
  data <- standardize_data(data)
  syn_size <- as.integer(if (is.null(syn_size)) ncol(data) * 5 else syn_size)

  syn_x_inform <- initial.auto_resampling(
    obs_x = data[, !names(data) %in% c(lhs(formula))],
    syn_size = syn_size,
    x_degree = if (is.null(x_degree)) rep(1, ncol(data)) else x_degree
  )

  syn_model <- do.call(stats::glm, as.list(match.call()[-1]))

  syn_y <- data.frame(y_value(fam)(generate_df_from_formula(
    formula,
    syn_x_inform$syn_x
  ),
  syn_model))
  colnames(syn_y) <- c(lhs(formula))

  syn_data <- cbind(syn_x_inform$syn_x, syn_y)

  initial_data <- list(
    syn_y_generation_formula = formula,
    family = family,
    obs_data = data,
    syn_size = syn_size,
    syn_data = syn_data,
    combined_data = rbind(data, syn_data),
    syn_model = syn_model,
    syn_x_generation_steps = syn_x_inform$recd_syn_x
  )

  class(initial_data) <- c(initial_data$class, "initialization")

  return(initial_data)
}

initial.auto_resampling <- function(obs_x, syn_size, x_degree, ...) {
  recd_syn_x <- data.frame(
    Type = character(),
    Process = character(),
    stringsAsFactors = FALSE
  )

  syn_x <- as.data.frame(apply(obs_x,
    2,
    sample,
    size = syn_size,
    replace = TRUE
  ))

  for (n in 1:ncol(syn_x)) {
    process <- "Coordinate"
    # Check continue or category
    if (is.continuous(obs_x[, n])) {
      data_type <- "Continuous"
      # Check order and diversity
      if (x_degree[n] >= length(unique(obs_x[, n])) - 1) {
        process <- c(process, "Smoothing")
        syn_x[, n] <- do.call(
          initial.resampling_methods,
          list(
            syn_x[, n],
            obs_x[, n],
            "smoothing"
          )
        )
      }
      # Check skewness
      if (abs(e1071::skewness(obs_x[, n])) >= 1) {
        process <- c(process, "IQR")
        syn_x[, n] <- do.call(
          initial.resampling_methods,
          list(
            syn_x[, n],
            obs_x[, n],
            "IQR"
          )
        )
      }
    } else {
      data_type <- "Category"
      tbl <- table(obs_x[, n])
      # Check data balance
      if (any(tbl < nrow(data) / (2 * length(tbl)))) {
        process <- c(process, "Flattening")
        syn_x[, n] <- do.call(
          initial.resampling_methods,
          list(
            syn_x[, n],
            obs_x[, n],
            "flattening"
          )
        )
      }
    }
    recd_syn_x <- rbind(recd_syn_x, data.frame(
      Type = data_type,
      Process = paste(process,
        collapse = "->"
      )
    ))
  }
  return(list(
    syn_x = syn_x,
    recd_syn_x = recd_syn_x
  ))
}



initial.resampling_methods <- function(x, sample_x, method) {
  switch(method,
    "symmetrizing" = mean(sample_x) + sample(c(-1, 1), length(x), replace = TRUE) * (x - mean(sample_x)),
    "smoothing" = x + stats::rnorm(length(x), 0, min(diff(sort(unique(sample_x)))) / 10),
    "flattening" = (function(x) {
      x[sample(length(x), length(x) / 2)] <- sample(unique(sample_x), length(x) / 2, replace = TRUE)
      return(x)
    })(x),
    "IQR" = (function(x) {
      x[sample(length(x), length(x) / 2)] <- rnorm(length(x) / 2, median(sample_x), -get_iqr(sample_x) / (2 * qnorm(1 / 4)))
      return(x)
    })(x),
    stop(paste(method, "Resampling Method Not available"))
  )
}

print.initialization <- function(x, ...) {
  cat("Formula for generating synthetic response:   ")
  print(x$formula)
  cat(
    "\nFamily: ", x$family$family,
    "\n\nSynthetic data size: ", x$syn_size,
    "\n\nSynthetic data generation steps:\n\n"
  )
  print(x$syn_x_generation_steps)
  cat("\n")
  invisible(x)
}
