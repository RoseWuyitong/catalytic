test_that("auto_resampling", {
  test_data_matrix <- matrix(
    data = c(
      10000, 0.1, rep(0, 5),
      2, 3, rep(0, 5),
      sample(c("red", "green", "blue"), 7, replace = TRUE),
      "a", rep("b", 6)
    ),
    ncol = 4
  )
  test_data_df <- as.data.frame(test_data_matrix)
  auto_resampling(
    data = test_data_matrix,
    resample_size = 20,
    data_degree = NULL,
    resample_only = FALSE
  )
})


test_that("resampling_methods", {
  test_x <- c(10000, rep(0, 9))
  test_x_skew <- abs(e1071::skewness(test_x))

  # Test symmetrizing can help with high skew vector
  expect_true(sum(replicate(5, test_x_skew >
    abs(e1071::skewness(resampling_methods(
      x = test_x,
      sample_x = test_x,
      method = "symmetrizing"
    )))), na.rm = TRUE) > 3)

  # Test smoothing can help with high skew vector
  expect_true(sum(replicate(5, test_x_skew >
    abs(e1071::skewness(resampling_methods(
      x = test_x,
      sample_x = test_x,
      method = "smoothing"
    )))), na.rm = TRUE) > 3)

  # Test flattening can help with high skew vector
  expect_true(sum(replicate(5, test_x_skew >
    abs(e1071::skewness(resampling_methods(
      x = test_x,
      sample_x = test_x,
      method = "flattening"
    )))), na.rm = TRUE) > 3)

  # Test IQR can help with high skew vector
  expect_true(sum(replicate(5, test_x_skew >
    abs(e1071::skewness(resampling_methods(
      x = test_x,
      sample_x = test_x,
      method = "IQR"
    )))), na.rm = TRUE) > 3)
})
