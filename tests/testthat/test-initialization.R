test_that("initialization: normal run", {
  data_ <- initialization(
    matrix(c(1, 4, 8, 2, 5, 8, 3, 6, 9), nrow = 3),
    c(1, 2, 3)
  )
  expect_equal(data_$syn_y, rep(2, 100))
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
