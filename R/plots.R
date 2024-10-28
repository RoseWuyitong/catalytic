#' Title
#'
#' @param digit
#' @param cat_tune
#'
#' @return
#' @export
#'
#' @examples
plot.cat_tune <- function(cat_tune, digit = 2, legend_pos = "topright", text_pos = 3) {
  if (cat_tune$function_name == "cat_glm_tune") {
    y_lab <- paste0("discrepancy error (", cat_tune$discrepancy_method, ")")
    y_plot <- cat_tune$risk_estimate_list
    main <- paste0("Plot of Risk Estimates from cat_glm_tune \n(", cat_tune$risk_estimate_method, ")")
    mark_point <- min(cat_tune$risk_estimate_list)
  } else {
    y_lab <- "likelihood score"
    y_plot <- cat_tune$likelihood_list
    main <- paste0("Plot of Likelihood Scores from cat_cox_tune \n(", cat_tune$likelihood_method, ")")
    mark_point <- max(cat_tune$likelihood_list)
  }


  plot(cat_tune$tau_seq,
    y_plot,
    type = "l",
    col = "black",
    lwd = 2,
    xlab = "tau",
    ylab = y_lab,
    main = main
  )

  # Add text labels to the highlighted point
  text(cat_tune$tau,
    mark_point,
    labels = paste0("(", round(cat_tune$tau, 1), ", ", round(mark_point, 1), ")"),
    pos = text_pos
  )

  # Draw red dotted lines from the axes to the highlighted point
  segments(x0 = par("usr")[1], y0 = mark_point, x1 = cat_tune$tau, y1 = mark_point, col = "red", lty = 2) # Horizontal line
  segments(x0 = cat_tune$tau, y0 = par("usr")[3], x1 = cat_tune$tau, y1 = mark_point, col = "red", lty = 2) # Vertical line

  points(cat_tune$tau,
    mark_point,
    pch = 4,
    col = "red"
  )

  legend(legend_pos,
    legend = "Optimal tau",
    col = "red",
    pch = 4,
    text.col = "black"
  )
}

# TODO: combine some that, clean the unnecessary class
#' Title
#'
#' @param cat_glm
#' @param digit
#'
#' @return
#' @export
#'
#' @examples
plot.cat_lmm_tune <- function(cat_glm, digit = 2) {
  plot(cat_glm$tau_seq,
    cat_glm$risk_estimate_collect,
    type = "l",
    col = "black",
    lwd = 2,
    xlab = "tau",
    ylab = paste0("risk estimate (", cat_glm$discrepancy_function, ")"),
    main = paste0("Plot of Risk Estimate from ", cat_glm$function_name, " (", cat_glm$risk_estimate_method, ")")
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
#' @param bayes_
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
traceplot.bayes <- function(bayes_, ...) {
  rstan::traceplot(bayes_$stan_sample_model, ...)
}


#' Title
#'
#' @param cat_glm
#'
#' @return
#' @export
#'
#' @examples
traceplot.cat_glm_bayes_joint_gibbs <- function(cat_glm_bayes_joint_gibbs_, pars = NULL, ...) {
  df <- cat_glm_bayes_joint_gibbs_$gibbs_collect

  if (is.null(pars)) {
    pars <- colnames(df)[1:min(9, ncol(df))]
  }

  if (length(pars) > 9) {
    stop("Please select less than 10 features to plot.")
  }

  par(mar = c(2, 2, 2, 2))
  par(mfrow = c(3, 3))

  for (col_name in pars) {
    plot(
      x = seq((cat_glm_bayes_joint_gibbs_$gibbs_warmup + 1), cat_glm_bayes_joint_gibbs_$gibbs_iter),
      y = df[, col_name],
      main = col_name,
      xlab = "",
      ylab = "",
      type = "l"
    )
  }
}

#' Title
#'
#' @param cat_lmm
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
traceplot.cat_lmm <- function(cat_lmm, pars = NULL, ...) {
  df <- cat_lmm$iter_collect

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
traceplot <- function(cat_, ...) {
  UseMethod("traceplot")
}

# contour_plot.bayes <- function(cat_glm,
#                                x,
#                                y,
#                                bin=15) {
# extract <- rstan::extract(cat_glm$sample_model)
# df <- data.frame(x=as.double(extract[[x]]),y=as.double(extract[[y]]))
#
# ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
#   ggplot2::geom_density_2d(bins = bin)
#
# }

likelihood_plot <- function(cat_cox, digit = 3) {
  plot(cat_cox$tau_seq,
    cat_cox$likelihood_collect,
    type = "l",
    col = "black",
    lwd = 2,
    xlab = "tau",
    ylab = paste0("risk estimate (", cat_cox$discrepancy_function, ")"),
    main = paste0("Plot of Risk Estimate from ", cat_cox$function_name, " (", cat_cox$likelihood_method, ")")
  )

  text(cat_cox$tau,
    cat_cox$max_likelihood,
    labels = round(cat_cox$max_likelihood, digit),
    pos = 3,
    offset = 0.5
  )

  points(cat_cox$tau,
    cat_cox$max_likelihood,
    pch = 4,
    col = "red"
  )
}
