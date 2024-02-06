test_that("frequentist_prior: normal run - parametric_boostrap", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3),
  )

  data_ <- frequentist_prior(data_, taus = 1:5)

  expect_true(grepl("risk estimate method = parametric_bootstrap",
    data_$est_beta[[1]]$method,
    fixed = TRUE
  ))
})


test_that("frequentist_prior: normal run - steinian_estimate", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    as.integer(c(1, 0, 1)),
    regression_type = "logistic",
    synthetic_response_generation_family = "binomial"
  )

  data_ <- frequentist_prior(data_,
    taus = 1:5,
    risk_estimate_method = "steinian_estimate"
  )

  expect_true(grepl("risk estimate method = steinian_estimate",
    data_$est_beta[[1]]$method,
    fixed = TRUE
  ))
})

test_that("frequentist_prior: normal run - cross_validation", {
  data_ <- initialization(
    matrix(rnorm(100, 1, 1), nrow = 25),
    rnorm(25, 0, 1),
  )

  data_ <- frequentist_prior(data_,
    taus = 1:5,
    risk_estimate_method = "cross_validation"
  )

  expect_true(grepl("risk estimate method = cross_validation",
    data_$est_beta[[1]]$method,
    fixed = TRUE
  ))
})

test_that("frequentist_prior: normal run - mallows_estimate", {
  data_ <- initialization(
    matrix(rnorm(100, 1, 1), nrow = 25),
    rnorm(25, 0, 1),
  )

  data_ <- frequentist_prior(data_, taus = 1:5, risk_estimate_method = "mallows_estimate")

  expect_true(grepl("risk estimate method = mallows_estimate",
    data_$est_beta[[1]]$method,
    fixed = TRUE
  ))
})


test_that("frequentist_prior: store all beta", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    as.integer(c(1, 0, 1)),
    regression_type = "logistic",
    synthetic_response_generation_family = "binomial"
  )

  data_ <- frequentist_prior(data_, taus = 1:5)
  data_ <- frequentist_prior(data_,
    taus = 1:5,
    risk_estimate_method = "steinian_estimate"
  )

  expect_true(grepl("risk estimate method = parametric_bootstrap",
    data_$est_beta[[1]]$method,
    fixed = TRUE
  ))
  expect_true(grepl("risk estimate method = steinian_estimate",
    data_$est_beta[[2]]$method,
    fixed = TRUE
  ))
})

test_that("frequentist_prior: not data_ stop", {
  data_ <- rlang::env()

  expect_error(frequentist_prior(data_))
})

test_that("frequentist_prior: single tau point estimate", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3)
  )

  data_ <- frequentist_prior(data_, taus = 1)

  expect_equal(data_$est_beta[[1]]$method, "single point estimate: tau =  1")
})

test_that("frequentist_prior: null tau", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3)
  )

  data_ <- frequentist_prior(data_)

  expect_equal(data_$est_beta[[1]]$method, "single point estimate: tau =  0.75")
})
