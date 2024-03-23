#' Title
#'
#' @param formula
#' @param family
#' @param data
#' @param weights
#' @param subset
#' @param na.action
#' @param offset
#' @param ...
#' @param syn_size
#' @param x_degree
#'
#' @return
#' @export
#'
#' @examples
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
  fam <- get_family(family)
  data <- standardize_data(data)
  syn_size <- as.integer(if (is.null(syn_size)) ncol(data) * 5 else syn_size)


  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "family", "subset", "weights", "na.action", "offset"),
    table = names(mf), nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf$data <- data

  syn_model <- do.call(stats::glm, as.list(mf[-1]))

  syn_x_inform <- auto_resampling(
    obs_x = data[, -which(names(data) == lhs(formula))],
    syn_size = syn_size,
    x_degree = if (is.null(x_degree)) rep(1, ncol(data)) else x_degree
  )

  syn_y <- data.frame(y_value(fam)(model.frame(rhs(formula, TRUE), syn_x_inform$syn_x),
    syn_model))
  colnames(syn_y) <- c(lhs(formula))

  data_ <- list(
    formula = formula,
    family = family,
    obs_size = nrow(data),
    obs_data = data,
    obs_x = data[, -which(names(data) == lhs(formula))],
    obs_y = as.data.frame(data[, which(names(data) == lhs(formula))]),
    syn_size = syn_size,
    syn_data = cbind(syn_x_inform$syn_x, syn_y),
    syn_x = syn_x_inform$syn_x,
    syn_y = syn_y,
    syn_model = syn_model,
    syn_x_generation_steps = syn_x_inform$recd_syn_x
  )
  data_$data <- rbind(data_$obs_data, data_$syn_data)
  names(data_$obs_y) <- as.character(lhs(formula))

  data_ <- c(
    data_,
    list(
      size = data_$obs_size + data_$syn_size,
      x = rbind(data_$obs_x, data_$syn_x),
      y = rbind(data_$obs_y, data_$syn_y)
    )
  )

  class(data_) <- c(data_$class, "initialization")

  return(data_)
}

auto_resampling <- function(obs_x, syn_size, x_degree, ...) {
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
          resampling_methods,
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
          resampling_methods,
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
          resampling_methods,
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



resampling_methods <- function(x, sample_x, method) {
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
