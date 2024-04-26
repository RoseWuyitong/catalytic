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

  data <- standardize_data(data)

  # Check if it is survival function or normal GLM funciton
  if (length(grep("Surv", get_lhs(formula))) == 0) {
    y_names <- get_lhs(formula)
    model <- do.call(stats::glm, list(
      formula = formula,
      family = family,
      data = data,
      ...
    ))
  } else {
    y_names <- unlist(strsplit(
      gsub(
        ".*Surv\\((.*?)\\)", "\\1",
        deparse(get_lhs(formula))
      ),
      ", "
    ))
    model <- do.call(
      survival::coxph,
      list(formula = formula,
        data = data,
        ...
      )
    )
  }


  obs_x <- as.data.frame(data[, setdiff(names(data), y_names)])
  obs_y <- as.data.frame(data[, y_names])

  x_degree <- if (is.null(x_degree)) rep(1, ncol(obs_x)) else x_degree
  syn_size <- as.integer(if (is.null(syn_size)) ncol(obs_x) * 5 else syn_size)
  syn_x_inform <- auto_resampling(
    obs_x = obs_x,
    syn_size = syn_size,
    x_degree = x_degree
  )

  syn_y <- data.frame(do.call(
    get_sim_y(fam),
    c(
      list(
        X = model.frame(
          get_rhs(formula, TRUE),
          syn_x_inform$syn_x
        ),
        model = model
      ))))

  colnames(syn_y) <- c(get_lhs(formula))

  data_ <- list(
    formula = formula,
    family = family,
    obs_size = nrow(data),
    obs_data = data,
    obs_x = obs_x,
    obs_y = obs_y,
    syn_size = syn_size,
    syn_data = cbind(syn_x_inform$syn_x, syn_y),
    syn_x = syn_x_inform$syn_x,
    syn_y = syn_y,
    syn_x_generation_steps = syn_x_inform$recd_syn_x,
    x_degree = x_degree
  )
  data_$data <- rbind(data_$obs_data, data_$syn_data)
  names(data_$obs_y) <- as.character(get_lhs(formula))

  data_ <- c(
    data_,
    list(
      size = data_$obs_size + data_$syn_size,
      x = rbind(data_$obs_x, data_$syn_x),
      y = rbind(data_$obs_y, data_$syn_y)
    )
  )

  class(data_) <- c(data_$class, "data_initialization")

  return(data_)
}
