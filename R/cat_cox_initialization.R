#' Initialize Cox Proportional-Hazards Model with Synthetic Data
#'
#' This function initializes data for fitting Cox Proportional-Hazards models with synthetic data using the catalytic package.
#' Adjust parameters such as the formula, syn_size, and x_degree to customize the initialization process.
#'
#' @param formula A formula specifying the model structure, which is used for generating synthetic response.
#'                See `survival::coxph` for details on formula syntax.
#' @param data The original dataset in `data.frame` format.
#' @param syn_size The size of the synthetic data generated. Default is four times the number of predictors.
#' @param hazard_constant A constant used in the hazard rate calculation. Default is the sum of observed statuses divided by the mean observed time divided by the number of observations.
#' @param entry_points An optional vector specifying the entry points for the data. Default is a vector of zeros.
#' @param x_degree Degrees of covariates in the synthetic data generation process. Defaults to NULL.
#' @param resample_only Logical indicating whether to generate synthetic data only by resampling original data. Default is FALSE.
#' @param na_replace Function used to handle missing values in the original dataset. Defaults to `stats::na.omit`.
#' @param ... Additional arguments to be passed to the `coxph` function.
#'
#' @return An object of class `cat_init` containing:
#'   \describe{
#'     \item{function_name}{Name of the initialization function.}
#'     \item{formula}{Formula used for model initialization.}
#'     \item{hazard_constant}{Hazard rate constant used.}
#'     \item{entry_points}{Entry points for the data.}
#'     \item{x_degree}{Degrees of covariates used in synthetic data generation.}
#'     \item{time_col_name}{Column name of the time variable.}
#'     \item{status_col_name}{Column name of the status variable.}
#'     \item{obs_size}{Size of the original dataset.}
#'     \item{obs_data}{Original dataset (observations).}
#'     \item{obs_x}{Covariates from the original dataset.}
#'     \item{obs_time}{Time variable from the original dataset.}
#'     \item{obs_status}{Status variable from the original dataset.}
#'     \item{syn_size}{Size of the synthetic dataset.}
#'     \item{syn_data}{Synthetic dataset generated (synthetic covariates and responses).}
#'     \item{syn_x}{Synthetic covariates.}
#'     \item{syn_time}{Synthetic times.}
#'     \item{syn_status}{Synthetic statuses.}
#'     \item{syn_x_resample_inform}{Information on resampled synthetic data.}
#'     \item{size}{Size of the original dataset plus the size of the synthetic dataset.}
#'     \item{data}{Combined observation and synthetic data.}
#'     \item{x}{Combined observation and synthetic covariates.}
#'     \item{time}{Combined observation and synthetic times.}
#'     \item{status}{Combined observation and synthetic statuses.}
#'   }
#' @export
#'
#' @examples
#' # Assuming `my_data` is a data frame with columns `time`, `status`, and covariates
#' formula <- survival::Surv(time, status) ~ .
#' cat_init <- cat_cox_initialization(
#'   formula = formula,
#'   data = my_data
#' )
#'
#' cat_init
cat_cox_initialization <- function(formula,
                                   data,
                                   syn_size = NULL,
                                   hazard_constant = NULL,
                                   entry_points = NULL,
                                   x_degree = NULL,
                                   resample_only = FALSE,
                                   na_replace = stats::na.omit) {
  # Validate Input Parameters
  validate_cox_initialization_input(
    formula = formula,
    data = data,
    syn_size = syn_size,
    hazard_constant = hazard_constant,
    entry_points = entry_points,
    x_degree = x_degree
  )

  # Extract Information from Input Parameters
  ## Standardize observation covariate and response
  obs_data <- get_standardized_data(
    data = data,
    na_replace = na_replace
  )

  ## Extract observation covariate and response
  f_lhs <- get_formula_lhs(formula)

  time_col_name <- gsub(
    "`",
    "",
    regmatches(
      f_lhs,
      regexec(
        "\\(([^,]+),",
        f_lhs
      )
    )[[1]][2]
  )

  status_col_name <- gsub(
    "`",
    "",
    regmatches(
      f_lhs,
      regexec(
        ",\\s*([^\\)]+)\\)",
        f_lhs
      )
    )[[1]][2]
  )

  obs_time <- obs_data[, time_col_name]
  obs_status <- obs_data[, status_col_name]
  obs_x <- obs_data[, setdiff(colnames(data), c(time_col_name, status_col_name)), drop = FALSE]
  obs_size <- nrow(obs_x)

  # Generate Synthetic Dataset
  if (is.null(x_degree)) x_degree <- rep(1, ncol(obs_x))
  if (is.null(syn_size)) syn_size <- 4 * ncol(obs_x)
  if (is.null(entry_points)) entry_points <- rep(0, obs_size)

  ## Generate synthetic covariate by resampling observation covariate
  syn_x_inform <- get_resampled_data(
    data = obs_x,
    resample_size = syn_size,
    data_degree = x_degree,
    resample_only = resample_only
  )
  syn_x <- as.data.frame(syn_x_inform$resampled_data)

  # Generate synthetic time by using simple model
  if (get_formula_rhs(formula) == 1) {
    simple_model <- sum(obs_status) / mean(obs_time) / obs_size
    syn_time <- rexp(
      syn_size,
      rate = simple_model
    )
  } else {
    simple_model <- survival::coxph(
      formula = formula,
      data = data
    )
    syn_time <- summary(
      survival::survfit(
        simple_model,
        newdata = as.data.frame(syn_x)
      )
    )$table[, "median"]
  }

  # Delete synthetic data that syn_time is NA
  syn_time_na_idx <- which(is.na(syn_time))

  if (length(syn_time_na_idx) > 0) {
    warning(paste(
      "Found and removed",
      length(syn_time_na_idx),
      "NA when calculate syn_y."
    ))

    syn_time <- syn_time[-syn_time_na_idx]
    syn_x <- syn_x[-syn_time_na_idx, ]
  }
  syn_status <- rep(1, syn_size - length(syn_time_na_idx))

  cat_init <- list(
    function_name = "cat_cox_initialization",
    ## Input/Processed parameters
    formula = formula,
    hazard_constant = if (is.null(hazard_constant)) sum(obs_status) / mean(obs_time) / obs_size else hazard_constant,
    entry_points = entry_points,
    x_degree = x_degree,
    resample_only = resample_only,
    na_replace = na_replace,
    time_col_name = time_col_name,
    status_col_name = status_col_name,
    simple_model = simple_model,
    # Observation data Information
    obs_size = obs_size,
    obs_data = obs_data,
    obs_x = obs_x,
    obs_time = obs_time,
    obs_status = obs_status,
    # Synthetic data Information
    syn_size = syn_size - length(syn_time_na_idx),
    syn_data = cbind(
      syn_x,
      setNames(as.data.frame(syn_time), time_col_name),
      setNames(as.data.frame(syn_status), status_col_name)
    ),
    syn_x = syn_x,
    syn_time = syn_time,
    syn_status = syn_status,
    syn_x_resample_inform = syn_x_inform$resampled_data_record
  )

  cat_init <- c(
    cat_init,
    list(
      ## Whole data information
      size = cat_init$obs_size + cat_init$syn_size,
      data = rbind(cat_init$obs_data, cat_init$syn_data),
      x = rbind(cat_init$obs_x, cat_init$syn_x),
      time = c(cat_init$obs_time, cat_init$syn_time),
      status = c(cat_init$obs_status, cat_init$syn_status)
    )
  )

  class(cat_init) <- c(cat_init$class, "cat_initialization")

  return(cat_init)
}
