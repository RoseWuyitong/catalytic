standardize_data <- function(data, if_missing = NULL) {
  if (missing(data) || is.null(data)) {
    warn_data_arg_missing()
    return(if_missing)
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  drop_redundant_dims(data)

  for (i in 1:ncol(data)) {
    if (!(NA %in% suppressWarnings(as.numeric(data[, i])))) {
      data[, i] <- as.numeric(data[, i])
    } else if (!(NA %in% suppressWarnings(as.logical(data[, i])))) {
      data[, i] <- as.logical(data[, i])
    } else {
      data[, i] <- as.factor(data[, i])
    }
  }

  return(data)
}

warn_data_arg_missing <- function() {
  warning(
    "Omitting the 'data' argument is not recommended ",
    "and may not be allowed in future versions of rstanarm. ",
    "Some post-estimation functions (in particular 'update', 'loo', 'kfold') ",
    "are not guaranteed to work properly unless 'data' is specified as a data frame.",
    call. = FALSE
  )
}


drop_redundant_dims <- function(data) {
  drop_dim <- sapply(data, function(v) is.matrix(v) && NCOL(v) == 1)
  data[, drop_dim] <- lapply(data[, drop_dim, drop = FALSE], drop)
  return(data)
}

validate_glm_formula <- function(f) {
  if (!any(grepl("~", f))) {
    stop("please enter ~ for formula",
      call. = FALSE
    )
  }
}




is.binomial <- function(x) x$family == "binomial"
is.gaussian <- function(x) x$family == "gaussian"

validate_family <- function(f) {
  if (is.character(f)) {
    f <- get(f, mode = "function", envir = parent.frame(2))
  }
  if (is.function(f)) {
    f <- f()
  }
  if (!is(f, "family")) {
    stop("'family' must be a family.", call. = FALSE)
  }

  return(f)
}

.validate_glm_formula <- function(f) {
  if (any(grepl("\\|", f))) {
    stop("Using '|' in model formula not allowed. ",
      "Maybe you meant to use 'stan_(g)lmer'?",
      call. = FALSE
    )
  }
}

rhs <- function(formula) {
  # Split the formula string by "~" and extract the first element
  return(rlang::f_rhs(formula))
}

lhs <- function(formula) {
  return(rlang::f_lhs(formula))
}

.validate_regression <- function(family) {
  if (.is.binomial(family)) {
    return("logistic")
  }
  return("linear")
}

.large_number_warning <- function(num, threshold, suggestion) {
  if (num > threshold) warning(suggestion)
}

.small_number_warning <- function(num, threshold, suggestion) {
  if (num < threshold) warning(paste(suggestion))
}

.not_same_stop <- function(num1, num2, suggestion) {
  if (num1 != num2) stop(suggestion)
}

.not_in_stop <- function(obj, expected_objs, suggestion) {
  if (!all(obj %in% expected_objs)) stop(suggestion)
}

.not_single_str_stop <- function(obj, suggestion) {
  if (!is.character(obj) | length(obj) != 1) stop(suggestion)
}

.large_number_ask <- function(num1, threshold, question, stop_ans) {
  if (num1 > threshold) {
    user_input <- readline(prompt = question)
    if (identical(user_input, stop_ans)) stop("program stop")
  }
}


is.continuous <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
}
