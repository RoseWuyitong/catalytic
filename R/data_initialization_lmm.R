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
                                    group_col = NULL,
                                    ...,
                                    syn_size = NULL,
                                    add_intercept = TRUE) {
  .not_true_stop(
    is.data.frame(data),
    "Please enter `data` in `data.frame` format"
  )

  tryCatch(
    {
      group_col_name <- if (is.null(group_col)) gsub("`",
                             "",
                             gsub(".*\\((.*?)\\|\\s*(\\S+)\\)",
                                  "\\2",
                                  deparse(formula))) else group_col
    }, error = function(e) {
      stop("Please enter `group_col`.", call. = FALSE)
    }
  )

  .not_true_stop(
    check_incul_wo_overlap(colnames(data), X_cols, Y_cols, Z_cols, group_col),
    "Please enter `X_cols`, `Y_cols` and `Z_cols` with no overlap,
    and inlcude all column names showed in `data.
    Or, please enter right group column name in the formula"
  )

  validate_glm_formula(formula)

  syn_size <- as.integer(if (is.null(syn_size)) length(X_cols) * 5 else syn_size)

  obs_data <- standardize_data(data)

  obs_x <- as.data.frame(obs_data[, X_cols])
  obs_group <- as.data.frame(obs_data[, group_col])
  obs_y <- as.data.frame(obs_data[, Y_cols])
  obs_z <- as.data.frame(obs_data[, Z_cols])

  syn_x_inform <- auto_resampling(
    obs_data = obs_x,
    resample_size = syn_size,
    data_degree = NULL
  )

  syn_group_inform <- auto_resampling(
    obs_data = obs_group,
    resample_size = syn_size,
    data_degree = NULL
  )

  syn_z_inform <- auto_resampling(
    obs_data = obs_z,
    resample_size = syn_size,
    data_degree = NULL
  )

  tryCatch({model <- do.call(
    lme4::lmer,
    list(
      formula = formula,
      data = obs_data,
      ...
    )
  )

  simple_model_coef <- lme4::fixef(model)

  X_cols_from_model <- as.formula(paste(
    sub(
      "\\+[^+]*$",
      "",
      get_rhs(formula, TRUE)
    ),
    collapse = " "
  ))},
  error = function(e){
    model <- do.call(stats::glm, list(
      formula = formula,
      family = "gaussian",
      data = data,
      ...
    ))

    simple_model_coef <- stats::coef(model)

    X_cols_from_model <- as.formula(get_rhs(formula, TRUE))
  })


  syn_y <- data.frame(do.call(
    get_sim_y("gaussian"),
    c(
      list(
        X = model.frame(
          X_cols_from_model,
          syn_x_inform$resampled_data
        )
      ),
      if (X_cols_from_model == ~1) {
        list(
          mean = as.double(lme4::fixef(simple_model_coef)),
          sd = 0
        )
      } else {
        list(
          coefficients = lme4::fixef(simple_model_coef),
          model = model
        )
      }
    )
  ))

  syn_x <- syn_x_inform$resampled_data
  colnames(syn_y) <- colnames(obs_y) <- Y_cols
  data_ <- list(
    function_name = "data_initilization_lmm",
    formula = formula,
    model = model,
    group_col_name = group_col_name,
    obs_size = nrow(obs_x),
    orig_obs_data = data,
    obs_x = if (add_intercept) {
      as.matrix(
        cbind("intercept" = 1, obs_x)
      )
    } else {
      as.matrix(obs_x)
    },
    obs_group_col = as.matrix(data[group_col_name]),
    obs_y = as.matrix(obs_y),
    obs_z = as.matrix(obs_z),
    syn_size = syn_size,
    syn_x = if (add_intercept) {
      as.matrix(
        cbind("intercept" = 1, syn_x[,-which(names(syn_x) == group_col_name)])
      )
    } else {
      as.matrix(
        syn_x[,-which(names(syn_x) == group_col_name)]
      )
    },
    syn_group_col = as.matrix(syn_group_inform$resampled_data),
    syn_y = as.matrix(syn_y),
    syn_z = as.matrix(syn_z_inform$resampled_data),
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

  class(data_) <- c(data_$class, "data_")

  return(data_)
}
