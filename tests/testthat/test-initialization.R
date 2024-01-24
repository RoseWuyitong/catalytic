test_that("initialization: normal run - linear regression", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3)
  )
  expect_equal(as.double(data_$syn_y_gen_model$coefficients), as.double(2))
})

test_that("initialization: normal run - logistic regression", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 0,1),
    synthetic_data_size = 30,
    synthetic_response_generation_family = "binomial",
    regression_type = "logistic",
    synthetic_response_generation_model = c(1,2)
  )
  expect_equal(as.double(data_$syn_y_gen_model$coefficients),
  as.double(c(149.2517687, 94.2642749, -109.9749873)))
})


test_that("initialization: nrow not same stop", {
  expect_error(initialization(
    matrix(1:9, nrow = 3),
    c(1)
  ))
})

test_that("initialization: length of degree not same stop", {
  expect_error(initialization(matrix(1:9, nrow = 3),
    c(1, 0, 1),
    observed_covariate_degree = c(1)
  ))
})
