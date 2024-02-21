cat_pe_glm <-
  function(formula,
           family = gaussian(),
           data,
           syn_size,
           tau = NULL,
           subset = NULL,
           na.action = NULL,
           offset = NULL,
           ...
           ) {
    validate_glm_formula(formula)
    family <- validate_family(family)
    syn_size <- as.integer(syn_size)

    if (is.null(tau)) tau <- ncol(data)/ 4

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "family", "data", "subset", "na.action", "offset"),
               table = names(mf), nomatch = 0L)
    mf <- mf[c(1L, m)]
    mf$weight <- c(rep(1, nrow(data)-syn_size), rep(tau / syn_size, syn_size))

    model <- do.call(stats::glm, as.list(mf[-1]))

    cat_pe <- list(
      formula = formula,
      family = family,
      data = data,
      syn_size = syn_size,
      tau = tau,
      model = model,
      est_beta = model$coefficient
    )

    class(cat_pe) <- c(cat_pe$class, "cat_pe_glm")
    return(cat_pe)
  }

print.cat_pe_glm <- function(x, ...) {
  cat("Formula for tau point estimate:   ")
  print(x$formula)
  cat(
    "\nFamily: ", x$family$family,
    "\n\nSynthetic data size: ", x$syn_size,
    "\n\ntau: ", x$tau,
    "\n"
  )
  invisible(x)
}
