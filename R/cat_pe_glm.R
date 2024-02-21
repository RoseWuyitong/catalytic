cat_pe_glm <-
  function(formula,
           family = gaussian(),
           data,
           syn_size,
           tau = NULL,
           subset = NULL,
           na.action = NULL,
           offset = NULL,
           ...,
           discrepancy_function = "square_error"
           ) {
    validate_glm_formula(formula)
    family <- validate_family(family)
    syn_size <- as.integer(syn_size)

    if (is.null(tau)) tau <- ncol(data)/ 4

    call <- match.call(expand.dots = TRUE)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "subset", "na.action", "offset"),
               table = names(mf), nomatch = 0L)
    mf <- mf[c(1L, m)]
    mf$weight <- c(rep(1, nrow(data)-syn_size), rep(tau / syn_size, syn_size))
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mf <- check_constant_vars(mf)
    mt <- attr(mf, "terms")
    return(mt)

    modl <- do.call(stats::glm, as.list(match.call()[-1]))

    return(list(

    ))
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
