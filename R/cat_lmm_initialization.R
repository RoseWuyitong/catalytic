#' Title
#'
#' @param formula
#' @param ...
#' @param data
#' @param x_cols
#' @param Y_cols
#' @param z_cols
#' @param syn_size
#'
#' @return
#' @export
#'
#' @examples
cat_lmm_initialization <- function(formula,
                                   data,
                                   x_cols,
                                   y_col,
                                   z_cols,
                                   group_col = NULL,
                                   syn_size = NULL,
                                   resample_by_group = FALSE,
                                   resample_only = FALSE,
                                   na_replace = mean) {
  if (is.null(group_cols)) {
    group_col <- gsub(
      "`",
      "",
      gsub(
        ".*\\((.*?)\\|\\s*(\\S+)\\)",
        "\\2",
        deparse(formula)
      )
    )
  }

  # Validate Input Parameters
  validate_lmm_initialization_input(
    formula = formula,
    data = data,
    x_cols = x_cols,
    y_col = y_col,
    z_cols = z_cols,
    group_col = group_col,
    syn_size = syn_size
  )

  # Extract Information from Input Parameters
  ## Standardize observation covariate and response
  obs_data <- get_standardized_data(
    data = data,
    na_replace = na_replace
  )

  ## Extract observation fixed effect, group, response and random effect
  obs_x <- as.data.frame(obs_data[, x_cols], dimnames = list(NULL, x_cols))
  obs_group <- obs_data[, group_cols]
  obs_y <- obs_data[, Y_cols]
  obs_z <- as.data.frame(obs_data[, z_cols], dimnames = list(NULL, z_cols))

  # Generate Synthetic Dataset
  if (is.null(syn_size)) syn_size <- length(x_cols) * 4

  ## Generating synthetic group by resampling observation group
  syn_group_inform <- get_resampled_data(
    data = as.data.frame(obs_group),
    resample_size = syn_size,
    data_degree = NULL,
    resample_only = TRUE
  )
  syn_group <- as.vector(syn_group_inform$resampled_data)

  ## Generate synthetic fixed effect and random effect
  if (resample_by_group) {
    syn_x <- as.data.frame(matrix(
      nrow = syn_size,
      ncol = dim(obs_x)[2],
      dimnames = list(NULL, c(x_cols))
    ))

    syn_z <- as.data.frame(matrix(
      nrow = syn_size,
      ncol = dim(obs_z)[2],
      dimnames = list(NULL, c(z_cols))
    ))

    for (g in unique(c(obs_group)[[1]])) {
      group_idx <- which(syn_group == g)

      syn_x_inform <- get_resampled_data(
        obs_data = obs_x[which(obs_group == g), ],
        resample_size = length(group_idx),
        data_degree = NULL
      )
      syn_x[group_idx, ] <- syn_x_inform$resampled_data

      syn_z_inform <- get_resampled_data(
        data = obs_z[which(obs_group == g), ],
        resample_size = length(group_idx),
        data_degree = NULL
      )
      syn_z[group_idx, ] <- syn_z_inform$resampled_data
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

  ## Generate synthetic response by using simple model
  simple_model <- if (grepl("\\(.*\\|.*\\)", deparse(formula))) {
    do.call(
      lme4::lmer,
      list(
        formula = formula,
        data = obs_data
      )
    )
  } else {
    do.call(
      stats::lm, list(
        formula = formula,
        data = data,
      )
    )
  }

  syn_data <- cbind(
    syn_x,
    syn_z,
    setNames(as.data.frame(syn_group), group_col)
  )

  syn_y <- as.data.frame(predict(simple_model, newdata = syn_data))

  # Finalize Setup and Output
  cat_init <- list(
    function_name = "cat_lmm_initialization",
    ## Input/Processed parameters
    formula = formula,
    data = data,
    x_cols = x_cols,
    y_col = y_col,
    z_cols = z_cols,
    group_col = group_col,
    syn_size = syn_size,
    resample_by_group = resample_by_group,
    resample_only = resample_only,
    na_replace = na_replace,
    model = simple_model,
    ## Observation data information
    obs_size = nrow(obs_x),
    obs_data = data,
    obs_x = obs_x,
    obs_y = obs_y,
    obs_z = obs_z,
    obs_group = obs_group,
    ## Synthetic data information
    syn_size = syn_size,
    syn_data = cbind(syn_data, setNames(as.data.frame(syn_y), y_col)),
    syn_x = syn_x,
    syn_y = syn_y,
    syn_z = syn_z,
    syn_group = syn_group,
    syn_x_resample_inform = syn_x_inform$resampled_data_record,
    syn_x_resample_inform = syn_z_inform$resampled_data_record
  )

  cat_init <- c(
    cat_init,
    list(
      ## Whole data information
      size = cat_init$obs_size + cat_init$syn_size,
      data = rbind(cat_init$obs_data, cat_init$syn_data),
      x = rbind(cat_init$obs_x, cat_init$syn_x),
      y = rbind(cat_init$obs_y, cat_init$syn_y),
      z = rbind(cat_init$obs_z, cat_init$syn_z),
      group = c(cat_init$obs_group, cat_init$syn_group)
    )
  )

  class(cat_init) <- c(cat_init$class, "cat_initialization")

  return(cat_init)
}
