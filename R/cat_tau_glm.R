#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param tau
#' @param subset
#' @param na.action
#' @param offset
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cat_tau_glm <-
  function(formula,
           family = gaussian(),
           data_,
           tau = NULL,
           subset = NULL,
           na.action = NULL,
           offset = NULL,
           ...) {
    .not_in_stop(
      c(
        "obs_size", "obs_data", "obs_x", "obs_y", "syn_size", "syn_data", "syn_x",
        "syn_y", "syn_model", "data", "size", "x", "y"
      ),
      names(data_),
      "Not all necessary variables are present in data_."
    )

    validate_glm_formula(formula)
    family <- validate_family(family)

    if (is.null(tau)) tau <- ncol(data_$data) / 4

    validate_tau(tau)

    call <- match.call(expand.dots = TRUE)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "family", "subset", "na.action", "offset"),
      table = names(mf), nomatch = 0L
    )
    mt <- mf[c(1L, m)]
    mt$data <- data_$data
    mt$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))

    tau_model <- do.call(stats::glm, as.list(mt[-1]))

    call["data_"] <- NULL
    call["tau"] <- NULL
    call["family"] <- NULL
    cat_tau <- c(list(
      model = tau_model,
      tau = tau,
      est_beta = tau_model$coefficient,
      data_ = data_,
      family = family
    ), as.list(call[-1]))

    class(cat_tau) <- c(cat_tau$class, "cat_tau_glm")
    return(cat_tau)
  }

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_tau_glm <- function(x, ...) {
  cat("\nFormula for tau point estimate:  ")
  print(x$formula)
  cat(
    "\nFamily: ", x$family$family,
    "\n\ntau: ", x$tau,
    "\n\nEstimated beta: ", x$est_beta
  )
  invisible(x)
}
