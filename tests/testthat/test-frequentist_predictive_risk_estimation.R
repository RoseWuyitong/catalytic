test_that("frequentist_predictive_risk_estimation: normal run", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3)
  )

  data_ <- frequentist_predictive_risk_estimation(data_)

  expect_equal(data_$est_beta[[1]]$method, "frequentist: omega estimate method = parametric_bootstrap ; discrepancy function = square_error ; bootstrap iteration times = 100 ; tau 0 =  1 ; optimal tau =  25")
})

test_that("frequentist_predictive_risk_estimation: store all beta", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3)
  )

  data_ <- frequentist_predictive_risk_estimation(data_)
  data_ <- frequentist_predictive_risk_estimation(data_, omega_estimate_method = "steinian_estimate")

  expect_equal(data_$est_beta[[1]]$method, "frequentist: omega estimate method = parametric_bootstrap ; discrepancy function = square_error ; bootstrap iteration times = 100 ; tau 0 =  1 ; optimal tau =  25")

  expect_equal(data_$est_beta[[2]]$method, "frequentist: omega estimate method = steinian_estimate ; discrepancy function = square_error ; bootstrap iteration times = 100 ; tau 0 =  1 ; optimal tau =  25")
})

test_that("frequentist_predictive_risk_estimation: not data_ stop", {
  data_ <- rlang::env()

  expect_error(frequentist_predictive_risk_estimation(data_))
})

test_that("frequentist_predictive_risk_estimation: wrong taus format stop", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3)
  )

  expect_error(frequentist_predictive_risk_estimation(data_, taus = c("a", "b")))
})
