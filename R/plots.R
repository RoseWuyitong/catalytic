#' Title
#'
#' @param cat_glm
#' @param digit
#'
#' @return
#' @export
#'
#' @examples
risk_estimate_plot.cat_glm_tune <- function(cat_glm, digit = 2) {
  plot(cat_glm$taus,
    cat_glm$risk_estimate_collect,
    type = "l",
    col = "black",
    lwd = 2,
    xlab = "taus",
    ylab = paste0("risk estimate (", cat_glm$discrepancy_function, ")"),
    main = paste0("Plot of Risk Estimate from cat_glm_tune (", cat_glm$risk_estimate_method, ")")
  )

  text(cat_glm$tau,
    cat_glm$min_risk_estimate,
    labels = round(cat_glm$min_risk_estimate, digit),
    pos = 3,
    offset = 0.5
  )

  points(cat_glm$tau,
    cat_glm$min_risk_estimate,
    pch = 4,
    col = "red"
  )
}

#' Title
#'
#' @param cat_glm_tune
#' @param digit
#'
#' @return
#' @export
#'
#' @examples
risk_estimate_plot <- function(cat_glm, digit = 2) {
  UseMethod("risk_estimate_plot")
}

#' Title
#'
#' @param cat_glm_bayes
#'
#' @return
#' @export
#'
#' @examples
traceplot.bayes <- function(cat_glm) {
  rstan::traceplot(cat_glm$sample_model)
}


#' Title
#'
#' @param cat_glm
#'
#' @return
#' @export
#'
#' @examples
traceplot.bayes_gibbs <- function(cat_glm, cols_to_plot = NULL) {
  df <- bayes_gibbs_5_model$gibbs_collect
  if (is.null(cols_to_plot)) {
    cols_to_plot <- colnames(df)[1:min(9, ncol(df))]
  }

  .large_number_stop(
    length(cols_to_plot),
    10,
    "Please select less than 10 features to plot."
  )
  par(mfrow = c(3, 3))

  for (col_name in cols_to_plot) {
    plot(df[, col_name], main = col_name, xlab = "iter", ylab = "Value", type = "l")
  }
}


#' Title
#'
#' @param cat_glm_bayes
#'
#' @return
#' @export
#'
#' @examples
traceplot <- function(cat_glm_bayes) {
  UseMethod("traceplot")
}

# contour_plot.bayes <- function(cat_glm,
#                                x,
#                                y,
#                                bin=15) {
#   extract <- rstan::extract(cat_glm$sample_model)
#   df <- data.frame(x=as.double(extract[[x]]),y=as.double(extract[[y]]))
#
#   ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
#     ggplot2::geom_density_2d(bins = bin)
#
# }
