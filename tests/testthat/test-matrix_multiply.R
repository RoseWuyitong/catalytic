test_that("matrix_multiply: with interval", {
  expect_equal(
    matrix_multiply(matrix(1:9, nrow = 3), c(1, 2, 3, 4)),
    1 + matrix(1:9, nrow = 3) %*% c(2, 3, 4)
  )
})

test_that("matrix_multiply: without interval", {
  expect_equal(
    matrix_multiply(matrix(1:9, nrow = 3), c(2, 3, 4), intl = FALSE),
    matrix(1:9, nrow = 3) %*% c(2, 3, 4)
  )
})
