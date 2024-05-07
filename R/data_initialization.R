#' Title
#'
#' @param formula
#' @param family
#' @param data
#' @param ...
#' @param syn_size
#' @param x_degree
#'
#' @return
#' @export
#'
#' @examples
data_initialization <- function(formula,
                                family = gaussian(),
                                data,
                                ...,
                                syn_size = NULL,
                                x_degree = NULL) {
  validate_glm_formula(formula)
  family <- validate_family(family)
  fam <- get_family(family)

  .not_true_stop(
    is.data.frame(data),
    "Please use `data` that in `data.frame` class."
  )

  data <- standardize_data(data)

  y_names <- as.character(get_lhs(formula))

  obs_x <- as.data.frame(data[, setdiff(names(data), c(y_names))])
  obs_y <- as.data.frame(data[, c(y_names)])

  x_degree <- if (is.null(x_degree)) rep(1, ncol(obs_x)) else x_degree
  syn_size <- as.integer(if (is.null(syn_size)) ncol(obs_x) * 5 else syn_size)

  syn_x_inform <- auto_resampling(
    resampled_data = obs_x,
    resample_size = syn_size,
    data_degree = data_degree
  )

  model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = data,
    ...
  ))
  syn_y <- data.frame(do.call(
    get_sim_y(fam),
    c(
      list(
        X = model.frame(
          get_rhs(formula, TRUE),
          syn_x_inform$resampled_data
        )
      ),
      if (get_rhs(formula, TRUE) == "~1") {
        list(
          mean = as.numeric(model$coefficients),
          sd = 0
        )
      } else {
        list(
          model = model
        )
      }
    )
  ))

  colnames(syn_y) <- colnames(obs_y) <- y_names

  data_ <- list(
    function_name = "data_initilization",
    formula = formula,
    family = family,
    obs_size = nrow(data),
    obs_data = data,
    obs_x = obs_x,
    obs_y = obs_y,
    syn_size = syn_size,
    syn_x = syn_x_inform$resampled_data,
    syn_y = syn_y,
    syn_x_generation_steps = syn_x_inform$recd_resampled_data,
    x_degree = x_degree
  )

  data_ <- c(
    data_,
    list(
      size = data_$obs_size + data_$syn_size,
      syn_data = cbind(data_$syn_x, data_$syn_y),
      x = rbind(data_$obs_x, data_$syn_x),
      y = rbind(data_$obs_y, data_$syn_y),
    )
  )

  data_$data <- rbind(data_$obs_data, data_$syn_data)

  class(data_) <- c(data_$class, "data_")

  return(data_)
}
