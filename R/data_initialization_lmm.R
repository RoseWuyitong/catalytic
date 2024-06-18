#' Title
#'
#' @param formula
#' @param ...
#' @param data
#' @param X_cols
#' @param Y_cols
#' @param Z_cols
#' @param syn_size
#'
#' @return
#' @export
#'
#' @examples
data_initialization_lmm <- function(formula,
                                    data,
                                    X_cols,
                                    Y_cols,
                                    Z_cols,
                                    ...,
                                    syn_size = NULL,
                                    by_group = FALSE) {
  if (!(is.glm_formula(formula) | is.lmm_formula(formula))) {
    stop(
      paste0(
        "The provided formula <",
        deparse(formula),
        "> is not suitable for a linear mixed model (LMM) or generalized linear model (GLM)."
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

  # Extract the column name of group from input formula
  group_col_name <- gsub(
    "`",
    "",
    gsub(
      ".*\\((.*?)\\|\\s*(\\S+)\\)",
      "\\2",
      deparse(formula)
    )
  )


  if (!check_subset_membership(
    colnames(data),
    X_cols,
    Y_cols,
    Z_cols,
    group_col_name
  )) {
    stop("Please enter `X_cols`, `Y_cols`, `Z_cols` and group_col with no overlap,
    and inlcude all column names showed in `data.
    Or, please enter right group column name in the formula")
  }

  if (is.null(syn_size)) syn_size <- length(X_cols) * 4

  # Standardize input data
  obs_data <- standardize_data(data)

  # Extract observation fixed effect, group, response and random effect
  obs_x <- as.data.frame(obs_data[, X_cols])
  obs_group <- as.data.frame(obs_data[, c(group_col_name)])
  obs_y <- as.data.frame(obs_data[, Y_cols])
  obs_z <- as.data.frame(obs_data[, Z_cols])

  # Generating synthetic covariate by resampling observation covariate
  syn_group_inform <- get_resampled_data(
    data = obs_group,
    resample_size = syn_size,
    data_degree = NULL,
    resample_only = TRUE
  )

  # Synthetic group
  syn_group <- syn_group_inform$resampled_data

  colnames(syn_group) <- colnames(obs_group) <- group_col_name

  if (by_group) {
    syn_x <- as.data.frame(matrix(
      nrow = syn_size,
      ncol = dim(obs_x)[2],
      dimnames = list(NULL, c(colnames(obs_x)))
    ))
    syn_z <- as.data.frame(matrix(
      nrow = syn_size,
      ncol = dim(obs_z)[2],
      dimnames = list(NULL, c(colnames(obs_z)))
    ))
    for (g in unique(c(obs_group)[[1]])) {
      idx <- which(syn_group == g)
      syn_x_inform <- get_resampled_data(
        obs_data = obs_x[which(obs_group == g), ],
        resample_size = length(idx),
        data_degree = NULL
      )

      syn_x[idx, ] <- syn_x_inform$resampled_data
      syn_z_inform <- get_resampled_data(
        data = obs_z[which(obs_group == g), ],
        resample_size = length(idx),
        data_degree = NULL
      )
      syn_z[idx, ] <- syn_z_inform$resampled_data
    }
    syn_x_inform <- NULL
    syn_z_inform <- NULL
  } else {
    syn_x_inform <- get_resampled_data(
      data = obs_x,
      resample_size = syn_size,
      data_degree = NULL
    )

    syn_z_inform <- get_resampled_data(
      data = obs_z,
      resample_size = syn_size,
      data_degree = NULL
    )

    syn_x <- syn_x_inform$resampled_data
    syn_z <- syn_z_inform$resampled_data
  }

  simple_model <- if (is.lmm_formula(formula)) {
    do.call(
      lme4::lmer,
      list(
        formula = formula,
        data = obs_data,
        ...
      )
    )
  } else {
    do.call(
      stats::lm, list(
        formula = formula,
        data = data,
        ...
      )
    )
  }

  syn_data <- as.data.frame(cbind(syn_x, syn_group, syn_z))

  syn_y <- as.data.frame(predict(simple_model, newdata = syn_data))

  colnames(syn_y) <- colnames(obs_y) <- Y_cols

  data_ <- list(
    function_name = "data_initilization_lmm",
    formula = formula,
    model = simple_model,
    group_col_name = group_col_name,
    obs_size = nrow(obs_x),
    orig_obs_data = data,
    obs_x = as.data.frame(obs_x),
    obs_group_col = as.data.frame(data[group_col_name]),
    obs_y = as.data.frame(obs_y),
    obs_z = as.data.frame(obs_z),
    syn_size = syn_size,
    syn_x = syn_x,
    syn_group_col = as.data.frame(syn_group),
    syn_y = as.data.frame(syn_y),
    syn_z = as.data.frame(syn_z),
    syn_x_generation_steps = syn_x_inform$recd_resampled_data,
    syn_z_generation_steps = syn_z_inform$recd_resampled_data
  )

  data_ <- c(
    data_,
    list(
      size = data_$obs_size + data_$syn_size,
      obs_data = cbind(data_$obs_x, data_$obs_y, data_$obs_z),
      syn_data = cbind(data_$syn_x, data_$syn_y, data_$syn_z),
      x = rbind(data_$obs_x, data_$syn_x),
      group_col = rbind(data_$obs_group_col, data_$syn_group_col),
      y = rbind(data_$obs_y, data_$syn_y),
      z = rbind(data_$obs_z, data_$syn_z)
    )
  )

  data_$data <- rbind(data_$obs_data, data_$syn_data)
  class(data_) <- c(data_$class, "data_initialization")

  return(data_)
}
