test_that("initial.resampling: symmetrizing method reduce skewness", {
  test_x_lst <- list(
    c(1000, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1000)
  )

  for (test_x in test_x_lst) {
    expect_true(
      abs(e1071::skewness(test_x)) > abs(e1071::skewness(initial.resampling(test_x, test_x, "symmetrizing")))
    )
  }
})

test_that("initial.resampling: smoothing method increase diversity", {
  test_x_lst <- list(
    c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )

  for (test_x in test_x_lst) {
    expect_true(
      length(unique(test_x)) < length(unique(initial.resampling(test_x, test_x, "smoothing")))
    )
  }
})

test_that("initial.resampling: flattening method increase balance", {
  test_x_lst <- list(
    as.factor(c("red", "green", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue")),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )

  for (test_x in test_x_lst) {
    expect_true(
      sd(table(test_x)) > sd(table(initial.resampling(test_x, test_x, "flattening")))
    )
  }
})

test_that("initial.resampling: IQR method reduce skewness", {
  test_x_lst <- list(
    c(-1000, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1000)
  )

  for (test_x in test_x_lst) {
    expect_true(
      abs(e1071::skewness(test_x)) > abs(e1071::skewness(initial.resampling(test_x, test_x, "IQR")))
    )
  }
})
