#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param tau
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cat_glm <-
  function(formula,
           family = gaussian(),
           data_,
           tau = NULL,
           ...) {
    .not_true_stop(
      is.data_(data_),
      paste(
        "Please enter `data_` that generated from",
        "catalytic::data_initialization."
      )
    )

    validate_glm_formula(formula)
    family <- validate_family(family)

    tau <- if (is.null(tau)) ncol(data_$x) / 4 else tau


    .small_number_stop(tau, 0, "Please enter non-negative `tau` value.")

    glm_model <- do.call(stats::glm, list(
      formula = formula,
      family = family,
      data = data_$data,
      weights = c(
        rep(1, data_$obs_size),
        rep(
          tau / data_$syn_size,
          data_$syn_size
        )
      )
    ))

    .large_number_warning(
      summary(glm_model)$coefficients[, grepl(
        "Pr.*",
        colnames(summary(glm_model)$coefficients)
      )],
      0.05,
      paste0(
        "Some(one) of the coefficient estimates in the model are(is) ",
        "not statistically significant at the ",
        "conventional significance level of 0.05.\n",
        "This suggests that there may not be sufficient evidence ",
        "to conclude that the variable(s) has a significant effect ",
        "on the response variables.\n",
        "SUGGESTION: Please consider utilizing `catalytic::cat_glm_tune` ",
        "or `catalytic::cat_glm_bayesian_joint` for better estimate."
      )
    )

    cat_glm_obj <- list(
      function_name = "cat_glm",
      formula = formula,
      family = family,
      data_ = data_,
      model = glm_model,
      tau = tau,
      coefficients = glm_model$coefficients
    )


    class(cat_glm_obj) <- c(cat_glm_obj$class, "none_bayes", "cat_glm")

    return(cat_glm_obj)
  }
