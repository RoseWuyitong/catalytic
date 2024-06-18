#' Title
#'
#' @param formula
#' @param family
#' @param data
#' @param ...
#' @param syn_size
#' @param x_degree
#' @param na_replace
#'
#' @return
#' @export
#'
#' @examples
data_initialization_glm <- function(formula,
                                    family = gaussian(),
                                    data,
                                    ...,
                                    syn_size = NULL,
                                    x_degree = NULL,
                                    na_replace = mean) {
  if (!is.glm_formula(formula)) {
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


  family <- validate_glm_family(family)

  # Extract the column name of response
  y_col_name <- get_formula_lhs(formula)

  # Standardize input data
  obs_data <- get_standardized_data(data, na_replace)

  # Extract observation covariate and response
  obs_x <- as.data.frame(obs_data[, setdiff(names(obs_data), c(y_col_name))])
  obs_y <- as.data.frame(obs_data[, c(y_col_name)])

  if (is.null(x_degree)) x_degree <- rep(1, ncol(obs_x))
  if (is.null(syn_size)) syn_size <- ncol(obs_x) * 4

  # Generating synthetic covariate by resampling observation covariate
  syn_x_inform <- get_resampled_data(
    data = obs_x,
    resample_size = syn_size,
    data_degree = x_degree
  )

  # Synthetic covariate
  syn_x <- syn_x_inform$resampled_data

  # Simple model for generate synthetic response
  simple_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = obs_data,
    ...
  ))

  # Synthetic response
  syn_y <- as.data.frame(predict(simple_model, syn_x))

  colnames(syn_y) <- colnames(obs_y) <- y_col_name

  data_ <- list(
    function_name = "data_initilization_glm",
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
    syn_data = cbind(syn_x, syn_y)[, colnames(data)],
    syn_x = syn_x,
    syn_y = syn_y,
    syn_x_resample_inform = syn_x_inform$resampled_data_record
  )

  data_ <- c(
    data_,
    list(
      # Observation + Synthetic data Information
      size = data_$obs_size + data_$syn_size,
      data = rbind(data_$obs_data, data_$syn_data),
      x = rbind(data_$obs_x, data_$syn_x),
      y = rbind(data_$obs_y, data_$syn_y)
    )
  )

  class(data_) <- c(data_$class, "data_initialization")

  return(data_)
}
