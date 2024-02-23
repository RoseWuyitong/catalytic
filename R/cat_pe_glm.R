cat_pe_glm <-
  function(formula,
           family = gaussian(),
           data_,
           tau = NULL,
           subset = NULL,
           na.action = NULL,
           offset = NULL,
           ...
           ) {
    validate_glm_formula(formula)
    family <- validate_family(family)

    if (is.null(tau)) tau <- ncol(data_$data)/ 4

    validate_tau(tau)

    call <- match.call(expand.dots = TRUE)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "family", "subset", "na.action", "offset"),
               table = names(mf), nomatch = 0L)
    mt <- mf[c(1L, m)]
    mt$data <- data_$data
    mt$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))

    pe_model <- do.call(stats::glm, as.list(mt[-1]))

    cat_pe <- list(
      model = pe_model,
      est_beta = pe_model$coefficient,
      data_ = data_,
      tau = tau,
      formula = formula,
      family = family)

    class(cat_pe) <- c(cat_pe$class, "cat_pe_glm")
    return(cat_pe)
  }

print.cat_pe_glm <- function(x, ...) {
  cat("Formula for tau point estimate:  ")
  print(x$formula)
  cat(
    "\nFamily: ", x$family,
    "\n\ntau: ", x$tau,
    "\n\nEstimated beta: ", x$est_beta
  )
  invisible(x)
}
