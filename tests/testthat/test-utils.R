test_that("multiply_matrix", {
  matrix1 <- matrix(data = c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  matrix2 <- matrix(data = c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4)
  df2 <- setNames(as.data.frame(matrix2), colnames(matrix2))
  matrix3 <- matrix(c(1))

  expected_m1_m2 <- matrix1 %*% matrix2

  expect_equal(multiply_matrix(matrix1, matrix2), expected_m1_m2)
  expect_equal(multiply_matrix(matrix2, matrix1), expected_m1_m2)
  expect_equal(multiply_matrix(matrix1, df2), expected_m1_m2)
  expect_error(multiply_matrix(matrix1, matrix3))
})
