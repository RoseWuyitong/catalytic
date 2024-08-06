#' Initialize data for catalytic GLM Logistic Regression
#'
#' This function initializes data for fitting GLM Logistic Regression models with synthetic data using the catalytic package.
#' Adjust parameters such as the formula, family type, and syn_size to customize the initialization process.
#'
#' @param formula A formula specifying the simple model structure, which is used for generating synthetic response.
#'                See `stats::glm` for details on formula syntax.
#' @param family The type of GLM family. Defaults to Gaussian.
#' @param data The original dataset in `data.frame` format.
#' @param syn_size The size of synthetic data to generate. Default is 4 times the number of columns in `data`.
#' @param resample_only Logical indicating whether to generate synthetic data only by resampling original data.
#' @param x_degree Degrees of covariates in the synthetic data generation process. Defaults to NULL.
#' @param na_replace Function used to handle missing values in the original dataset. Defaults to `stats::na.omit`.
#'
#' @return An object of class `cat_init` containing:
#'   \describe{
#'     \item{function_name}{Name of the initialization function.}
#'     \item{formula}{Formula used for model initialization.}
#'     \item{family}{GLM family type used.}
#'     \item{obs_size}{Size of the original dataset.}
#'     \item{obs_data}{Original dataset (observations).}
#'     \item{obs_x}{Covariates from the original dataset.}
#'     \item{obs_y}{Response variable from the original dataset.}
#'     \item{syn_size}{Size of the synthetic dataset.}
#'     \item{syn_data}{Synthetic dataset generated (synthetic covariates and responses).}
#'     \item{syn_x}{Synthetic covariates.}
#'     \item{syn_y}{Synthetic responses.}
#'     \item{syn_x_resample_inform}{Information on resampled synthetic data.}
#'     \item{size}{Size of the original dataset plus the size of the synthetic dataset.}
#'     \item{data}{Combined observation and synthetic data.}
#'     \item{x}{Combined observation and synthetic covariates.}
#'     \item{y}{Combined observation and synthetic responses.}
#'   }
#' @export
#'
#' @examples
#' data("binomial_swim")
#' swim_data <- cbind(binomial_swim$x, binomial_swim$y)
#'
#' # Initialize data for GLM Logistic Regression
#' cat_init <- cat_glm_initialization(
#'   formula = empyr1 ~ 1,
#'   family = binomial,
#'   data = swim_data,
#'   syn_size = 1000,
#'   resample_only = FALSE
#' )
#'
#' cat_init
cat_glm_initialization <- function(formula,
                                   family = gaussian(),
                                   data,
                                   syn_size = NULL,
                                   ...,
                                   x_degree = NULL,
                                   resample_only = FALSE,
                                   na_replace = stats::na.omit) {
  if (!is.formula_glm(formula)) {
    stop(
      paste0(
        "The provided formula <",
        deparse(formula),
        "> is not suitable for a generalized linear model (GLM)."
      )
    )
  }

  if (!is.data.frame(data)) {
    stop(
      paste0(
        "The provided data should be in `data.frame` format, not in <`",
        class(data)[1],
        "`> format."
      )
    )
  }


  family <- get_glm_family(family)

  # Extract the column name of response
  y_col_name <- get_formula_lhs(formula)

  # Extract and standardize observation covariate and response
  data <- get_standardized_data(
    data = data,
    na_replace = na_replace
  )

  obs_x <- as.data.frame(data[, setdiff(names(data), c(y_col_name))])
  obs_y <- as.data.frame(data[, c(y_col_name)])
  colnames(obs_y) <- y_col_name

  obs_data <- cbind(obs_x, obs_y)

  if (is.null(x_degree)) x_degree <- rep(1, ncol(obs_x))
  if (is.null(syn_size)) syn_size <- ncol(obs_x) * 4

  # Generate synthetic covariate by resampling observation covariate
  syn_x_inform <- get_resampled_data(
    data = obs_x,
    resample_size = syn_size,
    data_degree = x_degree,
    resample_only = resample_only
  )

  syn_x <- syn_x_inform$resampled_data


  # Generate synthetic response by using simple model
  simple_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = obs_data,
    ...
  ))

  if (family == "binomial") {
    prob <- stats::predict(simple_model,
      syn_x,
      type = if (family == "binomial") "response" else "link"
    )

    syn_y <- as.data.frame(get_glm_sample_data(
      family = family,
      n = length(prob),
      mean = prob
    ))

    obs_y <- apply(obs_y, 2, function(x) as.integer(as.character(x)))
    syn_y <- apply(syn_y, 2, function(x) as.integer(as.character(x)))
  } else {
    syn_y <- as.data.frame(stats::predict(
      simple_model,
      syn_x
    ))

    obs_y <- apply(obs_y, 2, function(x) as.numeric(as.character(x)))
    syn_y <- apply(syn_y, 2, function(x) as.numeric(as.character(x)))
  }

  colnames(syn_y) <- y_col_name

  cat_init <- list(
    function_name = "cat_glm_initialization",
    formula = formula,
    family = family,
    x_degree = x_degree,
    # Observation data Information
    obs_size = nrow(obs_data),
    obs_data = obs_data,
    obs_x = obs_x,
    obs_y = obs_y,
    # Synthetic data Information
    syn_size = syn_size,
    syn_data = cbind(syn_x, syn_y),
    syn_x = syn_x,
    syn_y = syn_y,
    syn_x_resample_inform = syn_x_inform$resampled_data_record,
    ...
  )

  cat_init <- c(
    cat_init,
    list(
      # Observation + Synthetic data Information
      size = cat_init$obs_size + cat_init$syn_size,
      data = rbind(cat_init$obs_data, cat_init$syn_data),
      x = rbind(cat_init$obs_x, cat_init$syn_x),
      y = rbind(cat_init$obs_y, cat_init$syn_y)
    )
  )

  class(cat_init) <- c(cat_init$class, "cat_initialization")

  return(cat_init)
}
