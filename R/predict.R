#' Predict Outcome for New Data Using a Fitted GLM Model
#'
#' This function generates predictions for new data points based on a fitted categorical Generalized Linear Model (GLM) object.
#' Depending on the type of model, it either uses `stats::predict.glm` or calculates predictions based on the model coefficients.
#'
#' @param cat_model A fitted model object of class `cat_glm`, containing the GLM fit and model details.
#' @param newdata An optional data frame containing new predictor values. If `NULL`,
#' the function uses the observation data from the model's initialization object.
#' @param ... Additional arguments passed to `stats::predict.glm`, if applicable.
#'
#' @return A vector of predicted values for the specified new data.
#' @export
predict.cat_glm <- function(
    cat_model,
    newdata = NULL,
    ...) {
  newdata <- if (is.null(newdata)) cat_model$cat_init$adj_obs_x else newdata[, colnames(cat_model$cat_init$adj_x), drop = FALSE]

  if (cat_model$function_name %in% c("cat_glm", "cat_glm_tune")) {
    return(
      stats::predict.glm(
        cat_model$model,
        newdata = as.data.frame(newdata),
        ...
      )
    )
  }

  return(
    c(get_glm_mean(
      family_string = cat_model$cat_init$family,
      X = newdata,
      coefs = stats::coef(cat_model)
    ))
  )
}

#' Predict Linear Predictor for New Data Using a Fitted Cox Model
#'
#' This function calculates the linear predictor (LP) for new data points based on a fitted Cox proportional hazards model.
#'
#' @param cat_model A fitted model object of class `cat_cox`, containing the COX fit and model details.
#' @param newdata An optional data frame with new predictor values. If `NULL`,
#' the function uses the observation data from the model's initialization object.
#'
#' @return A vector of linear predictor values for the specified new data.
#' @export
predict.cat_cox <- function(
    cat_model,
    newdata = NULL) {
  # Prepare the covariates
  newdata <- if (is.null(newdata)) cat_model$cat_init$adj_obs_x else newdata[, colnames(cat_model$cat_init$adj_x), drop = FALSE]

  return(c(get_linear_predictor(
    X = newdata,
    coefs = stats::coef(cat_model)
  )))
}

#' Predict Linear Predictor for New Data Using a Fitted Linear Mixed Model
#'
#' This function calculates the linear predictor (LP) for new data points based on a fitted linear mixed model (LMM) stored in `cat_model`.
#'
#' @param cat_model A fitted model object of class `cat_lmm`, containing the LMM fit and model details.
#' @param newdata An optional data frame with new predictor values. If `NULL`,
#' the function uses the observation data from the model's initialization object.
#'
#' @return A vector of linear predictor values for the specified new data.
#' @export
predict.cat_lmm <- function(
    cat_model,
    newdata = NULL) {
  # Prepare the covariates
  newdata <- if (is.null(newdata)) cat_model$cat_init$obs_x else newdata[, colnames(cat_model$cat_init$x), drop = FALSE]

  return(c(get_linear_predictor(
    X = newdata,
    coefs = stats::coef(cat_model)
  )))
}
