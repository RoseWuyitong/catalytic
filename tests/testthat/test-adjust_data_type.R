test_that("adjust_data_type: data type transformation", {
  test_matrix <- matrix(c(
    "1", "0", "1", "1", "0", "0",
    "TRUE", "FALSE", "TRUE", "FALSE", "TRUE", "TRUE",
    "red", "blue", "green", "red", "green", "green",
    "1.1", "2.0", "3.1", "4.2", "5.6", "6.0"
  ), ncol = 4)
  target_df <- data.frame(
    V1 = c(1, 0, 1, 1, 0, 0),
    V2 = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
    V3 = as.factor(c("red", "blue", "green", "red", "green", "green")),
    V4 = c(1.1, 2.0, 3.1, 4.2, 5.6, 6.0)
  )
  expect_equal(adjust_data_type(test_matrix), target_df)
})
