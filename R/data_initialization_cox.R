data_initialization_cox <- function(formula,
                                    data,
                                    ...,
                                    syn_size = NULL,
                                    x_degree = NULL) {
  .not_true_stop(
    is.data.frame(data),
    "Please use `data` that in `data.frame` class."
  )

  .not_true_stop(
    is.surv_formula(formula),
    "Please use `Surv` for formula."
  )

  # Extract time, event and other columns(x) from `data`
  y_lhs <- as.character(get_lhs(formula))
  time_col_name <- gsub("`", "", y_lhs[2])
  event_col_name <- gsub("`", "", y_lhs[3])

  obs_time <- matrix(data[, time_col_name], dimnames = list(NULL, time_col_name))
  obs_event <- matrix(data[, event_col_name], dimnames = list(NULL, event_col_name))
  obs_x <- as.matrix(data[, setdiff(
    colnames(data),
    c(time_col_name, event_col_name)
  )])
  obs_size <- nrow(data)

  x_degree <- if (is.null(x_degree)) rep(1, ncol(obs_x)) else x_degree
  syn_size <- as.integer(if (is.null(syn_size)) ncol(obs_x) * 5 else syn_size)


  # Generate synthetic data
  syn_x_inform <- auto_resampling(
    data = obs_x,
    resample_size = syn_size,
    data_degree = x_degree
  )
  syn_x <- as.matrix(syn_x_inform$resampled_data)

  # Generate synthetic time

  if (get_rhs(formula) == 1) {
    syn_time <- rexp(syn_size,
      rate = sum(obs_event) / mean(obs_time) / obs_size
    )
  } else {
    simple_model <- do.call(
      survival::coxph,
      list(formula = formula, data = data)
    )

    syn_time <- summary(survival::survfit(simple_model,
      newdata = as.data.frame(syn_x)
    ))$table[, "median"]
  }

  syn_NA_idx <- as.integer(which(is.na(syn_time)))
  #
  #   # Calculate the baseline hazard function for the Cox model
  #   baseline_survival <- basehaz(simple_model, centered = FALSE)
  #   baseline_time <- baseline_survival$time
  #
  #   # Get Baseline Hazard and Cumulative Baseline Hazard
  #   baseline_hazard <- baseline_survival$hazard
  #   cumulative_baseline_hazard <- cumsum(baseline_hazard)
  #   browser()
  #   syn_time <- matrix(0, nrow = syn_size, dimnames = list(NULL, time_col_name))
  #   # for each predicted linear predictor
  #   for (i in seq_along(linear_preds)) {
  #     risk_score <- exp(linear_preds[i])
  #
  #     # Find the time where the cumulative hazard matches the expected value
  #     expected_hazard <- -log(runif(1)) / risk_score
  #     predicted_time <- baseline_time[which.min(abs(cumulative_baseline_hazard - expected_hazard))]
  #     syn_time[i] <- predicted_time
  #   }
  if (length(syn_NA_idx) > 0) {
    warning(paste("Found", length(syn_NA_idx), "NA when calculate syn_y, removed."))
    syn_event <- matrix(rep(1, syn_size)[-syn_NA_idx], dimnames = list(NULL, event_col_name))
    syn_x <- syn_x[-syn_NA_idx, ]
    syn_time <- matrix(syn_time[-syn_NA_idx],
      ncol = 1, dimnames = list(NULL, time_col_name)
    )
  } else {
    syn_event <- matrix(rep(1, syn_size), dimnames = list(NULL, event_col_name))
    syn_time <- matrix(syn_time,
      ncol = 1, dimnames = list(NULL, time_col_name)
    )
  }

  syn_data <- as.matrix(cbind(syn_x, syn_time, syn_event))[, colnames(data)]
  data_ <- list(
    function_name = "data_initilization_cox",
    formula = formula,
    time_col_name = time_col_name,
    event_col_name = event_col_name,
    obs_size = obs_size,
    obs_data = data,
    obs_x = obs_x,
    obs_time = obs_time,
    obs_event = obs_event,
    syn_size = syn_size - length(syn_NA_idx),
    syn_data = syn_data,
    syn_x = syn_x,
    syn_time = syn_time,
    syn_event = syn_event,
    syn_x_generation_steps = syn_x_inform$resampled_data_record,
    x_degree = x_degree
  )

  data_ <- c(
    data_,
    list(
      size = data_$obs_size + data_$syn_size,
      x = as.matrix(rbind(data_$obs_x, data_$syn_x)),
      time = as.matrix(rbind(data_$obs_time, data_$syn_time)),
      event = as.matrix(rbind(data_$obs_event, data_$syn_event))
    )
  )

  data_$data <- rbind(data_$obs_data, data_$syn_data)

  class(data_) <- c(data_$class, "data_")

  return(data_)
}
