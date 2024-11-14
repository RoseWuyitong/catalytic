cox_data <- mock_cox_data()
cat_init <- mock_cat_cox_initialization(cox_data)

chains <- 1
iter <- 1000
warmup <- 500

suppressWarnings(
  cat_model <- cat_cox_bayes(
    formula = ~.,
    cat_init = cat_init,
    chains = chains,
    iter = iter,
    warmup = warmup
  )
)
test_that("cat_cox_bayes runs without errors for valid input", {
  expected_components <- c(
    "function_name", "cat_init", "tau", "hazard_beta",
    "chains", "iter", "warmup", "stan_data", "stan_model",
    "stan_sample_model", "coefficients", "increment_cumulative_baseline_hazard"
  )

  cat_model$cat_init$adj_obs_x <- cat_model$cat_init$adj_syn_x <- cat_model$cat_init$adj_x <- NULL

  expect_type(cat_model, "list")
  expect_equal(cat_model$function_name, "cat_cox_bayes")
  expect_equal(cat_model$cat_init, cat_init)
  expect_true(all(expected_components %in% names(cat_model)))

  expect_equal(cat_model$chains, chains)
  expect_equal(cat_model$iter, iter)
  expect_equal(cat_model$warmup, warmup)
  expect_true(inherits(cat_model$stan_model, "stanmodel"))
  expect_true(inherits(cat_model$stan_sample_model, "stanfit"))
})

test_that("cat_cox_bayes returns different results with different tau values", {
  cat_model_low_tau <- cat_cox_bayes(
    formula = ~.,
    cat_init = cat_init,
    tau = 0.1,
    chains = chains,
    iter = iter,
    warmup = warmup
  )
  suppressWarnings(
    cat_model_high_tau <- cat_cox_bayes(
      formula = ~.,
      cat_init = cat_init,
      tau = 100,
      chains = chains,
      iter = iter,
      warmup = warmup
    )
  )

  expect_true(all(coef(cat_model_low_tau) != coef(cat_model_high_tau)))
})

test_that("cat_cox_bayes can provide stable result with high iteration", {
  result_coefs <- matrix(NA, nrow = 10, ncol = 3)
  for (i in 1:10) {
    suppressWarnings(
      cat_model <- cat_cox_bayes(
        formula = ~.,
        cat_init = cat_init,
        chains = 2,
        iter = 5000,
        warmup = 3000
      )
    )

    result_coefs[i, ] <- coef(cat_model)
  }

  expect_true(all(apply(result_coefs, 2, sd) < 0.1))
})

test_that("cat_cox_bayes can provide stable result with same observation data but different synthetic data", {
  result_coefs <- matrix(NA, nrow = 10, ncol = 3)
  for (i in 1:10) {
    cat_init <- mock_cat_cox_initialization(cox_data)
    suppressWarnings(
      cat_model <- cat_cox_bayes(
        formula = ~.,
        cat_init = cat_init,
        chains = 2,
        iter = 5000,
        warmup = 3000
      )
    )
    result_coefs[i, ] <- coef(cat_model)
  }

  expect_true(all(apply(result_coefs, 2, sd) < 0.5))
})
