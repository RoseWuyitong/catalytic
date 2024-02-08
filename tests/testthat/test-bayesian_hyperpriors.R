test_that("bayesian_hyperpriors: linear none-joint - normal run", {
  data_ <- initialization(
    matrix(rnorm(100, 1, 1), nrow = 25),
    rnorm(25, 0, 1),
  )

  data_ <- bayesian_hyperpriors(data_,
    iter = 3000
  )

  expect_equal(data_$est_beta[[1]]$ess_pass, "PASS")
  expect_equal(data_$est_beta[[1]]$rhat_pass, "PASS")
})

test_that("bayesian_hyperpriors: logistic none-joint - normal run", {
  data_ <- initialization(
    matrix(rnorm(100, 1, 1), nrow = 25),
    sample(c(0, 1), 25, replace = TRUE),
    synthetic_response_generation_family = "binomial",
    regression_type = "logistic"
  )

  data_ <- bayesian_hyperpriors(data_,
    joint = TRUE
  )

  expect_equal(data_$est_beta[[1]]$ess_pass, "PASS")
  expect_equal(data_$est_beta[[1]]$rhat_pass, "PASS")
})

test_that("bayesian_hyperpriors: linear joint - normal run", {
  obs_x <- matrix(rnorm(100, 1, 1), nrow = 25)
  obs_y <- obs_x %*% c(1,-1,1,-1) + 1

  data_ <- initialization(
    obs_x,
           obs_y)

  data_ <- bayesian_hyperpriors(data_, joint = TRUE)

  expect_equal(data_$est_beta[[1]]$ess_pass, "PASS")
  expect_equal(data_$est_beta[[1]]$rhat_pass, "PASS")
})

test_that("bayesian_hyperpriors: logistic joint - weak run", {
  data_ <- initialization(
    matrix(rnorm(100, 1, 1), nrow = 25),
    sample(c(0, 1), 25, replace = TRUE),
    synthetic_response_generation_family = "binomial",
    regression_type = "logistic"
  )

  data_ <- bayesian_hyperpriors(data_,
    joint = TRUE,
    iter = 30,
    warmup = 10
  )

  expect_equal(data_$est_beta[[1]]$ess_pass, "FAIL")
  expect_equal(data_$est_beta[[1]]$rhat_pass, "FAIL")
})
