initialization <- function(observed_covariate,
                           observed_response,
                           synthetic_data_size = 100,
                           observed_covariate_degree = NULL,
                           synthetic_response_generation = NULL) {
  data_ <- rlang::env(
    obs_x = as.data.frame(observed_covariate),
    obs_y = as.data.frame(observed_response),
    syn_n = as.integer(synthetic_data_size),
    syn_x = NA,
    syn_y = NA,
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    )
  )

  .not_same_stop(
    nrow(data_$obs_x),
    nrow(data_$obs_y),
    "Row number of observed_covariate and observed_response not same."
  )

  data_$obs_x <- adjust_data_type(data_$obs_x)
  data_$obs_x_row_n <- nrow(data_$obs_x)
  data_$x_col_n <- ncol(data_$obs_x)

  if (is.null(X_degree)) {
    X_degree <- rep(1, data_$x_col_n)
  } else if (is.double(X_degree)) {
    .not_same_stop(
      data_$x_col_n,
      length(X_degree),
      "length of degree of matrix is not same as the column number of matrix."
    )
  } else {
    stop("Please enter degree in right format, e.g. c(1,2,1,0,...)")
  }

  data_$x_deg <- X_degree

  # Generate synthetic covariate
  message("generating synthetic covariate")
  data_ <- do.call("initial.auto_resampling", list(data_))

  # Generate synthetic response


  return(data_)
}


initial.auto_resampling <- function(data_) {
  data_$syn_x <- apply(data_$obs_x,
    2,
    sample,
    size = data_$syn_n,
    replace = TRUE
  )

  for (n in 1:ncol(X)) {
    message(paste("Column", n, "is processing"))
    process <- "Coordinate"
    # Check continue or category
    if (is_continue(X[, n])) {
      data_type <- "Continue"
      # Check order and diversity
      if (X_degree[n] >= unique(X[, n]) - 1) {
        process <- c(process, "Smoothing")
        X[, n] <- do.call("initial.resampling", list(X[, n], "smoothing"))
      }
      # Check skewness
      if (abs(e1071::skewness(X[, n])) >= 1) {
        process <- c(process, "IQR")
        X[, n] <- do.call("initial.resampling", list(X[, n], "IQR"))
      }
    } else {
      data_type <- "Category"
      tbl <- table(X[, n])
      thrh <- nrow(X) / (2 * length(tbl))
      for (u in names(tbl)) {
        # Check data balance
        if (tbl[u] < thrh) {
          process <- c(process, "Flattening")
          X[, n] <- do.call("initial.resampling", list(X[, n], "flattening"))
          break
        }
      }
    }
    data_$recd_x <- rbind(data_$recd_x, data.frame(Type = data_type, Process = paste(process, collapse = "->")))
  }
  return(data_)
}



initial.resampling <- function(x, method) {
  message(paste(method, "resampling"))
  switch(method,
    "symmetrizing" = mean(x) + sample(c(-1, 1), length(x), replace = TRUE) * (mean(x) - x),
    "smoothing" = x + stats::rnorm(length(x), 0, min(diff(sort(unique(x)))) / 10),
    "flattening" = unlist(lapply(x, function(x) {
      x[sample(length(x), length(x) / 2)] <- sample(unique(x), length(x) / 2, replace = TRUE)
      return(x)
    })),
    "IQR" = sapply(x, function(i) sample(c(i, rnorm(1, median(x), -IQR(x) / (2 * qnorm(1 / 4)))), 1)),
    stop(paste(method, "Resampling Method Not available"))
  )
}

initial.generate_data_from_dist <- function(from_dist = "uniform",
                                            sample_size = 10,
                                            dist_mean = 0,
                                            dist_sd = 1) {
  as.matrix(
    switch(from_dist,
      "bernoulli" = stats::rbinom(
        sample_size,
        1,
        dist_mean
      )
    ),
    "gaussian" = stats::rnorm(
      sample_size,
      dist_mean,
      dist_sd
    ),
    "uniform" = rep(dist_mean, sample_size)
  )
}
