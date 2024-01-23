test_that("initial.auto_resampling: data type", {
  test_x <- data.frame(
    "1" = c(1, 0, 1, 0, 1, 0),
    "2" = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
    "3" = as.factor(c("red", "green", "blue", "red", "green", "blue")),
    "4" = c(1, 2, 3, 4, 5, 1),
    "5" = c(1.1, 2.2, 1.3, 2.4, 1.5, 2.6)
  )

  data_ <- rlang::env(
    syn_x = test_x,
    obs_x = test_x,
    syn_n = 10,
    x_deg = c(1, 1, 1, 1, 1),
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    )
  )

  data_ <- initial.auto_resampling(data_)
  expect_equal(data_$recd_x[, "Type"], c(
    "Category",
    "Category",
    "Category",
    "Continuous",
    "Continuous"
  ))
})

test_that("initial.auto_resampling: no trigger", {
  test_x <- data.frame(
    "1" = c(1, 0, 1, 0, 1, 0),
    "2" = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
    "3" = as.factor(c("red", "green", "blue", "red", "green", "blue")),
    "4" = c(1, 2, 3, 4, 5, 1),
    "5" = c(1.1, 2.2, 1.3, 2.4, 1.5, 2.6)
  )
  data_ <- rlang::env(
    syn_x = test_x,
    obs_x = test_x,
    syn_n = 10,
    x_deg = c(0, 1, 2, 3, 1),
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    )
  )

  data_ <- initial.auto_resampling(data_)
  expect_equal(data_$recd_x[, "Process"], c(
    "Coordinate",
    "Coordinate",
    "Coordinate",
    "Coordinate",
    "Coordinate"
  ))
})


test_that("initial.auto_resampling: smoothing trigger", {
  test_x <- data.frame(
    "1" = c(1, 2, 1, 2, 1, 2),
    "2" = c(1, 2, 1, 2, 1, 2),
    "3" = c(1, 2, 1, 2, 1, 2),
    "4" = c(1, 2, 1, 2, 1, 2),
    "5" = c(1, 2, 3, 1, 2, 3)
  )
  data_ <- rlang::env(
    syn_x = test_x,
    obs_x = test_x,
    syn_n = 10,
    x_deg = c(1, 1, 1, 1, 2),
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    )
  )

  data_ <- initial.auto_resampling(data_)
  expect_equal(data_$recd_x[, "Process"], c(
    "Coordinate",
    "Coordinate",
    "Coordinate",
    "Coordinate",
    "Coordinate->Smoothing"
  ))
})


test_that("initial.auto_resampling: IQR trigger", {
  test_x <- data.frame(
    "1" = c(1000, 1000, 1000, -1000, -1000, -1000),
    "2" = c(-1000, -1000, -1000, 1000, 1000, 1000),
    "3" = c(1, 2, 1, 2, 1, 2),
    "4" = c(0, 0, 0, 0, 1, 100000),
    "5" = c(100000, 1, 0, 0, 0, 0)
  )
  data_ <- rlang::env(
    syn_x = test_x,
    obs_x = test_x,
    syn_n = 10,
    x_deg = c(1, 1, 1, 1, 1),
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    )
  )

  data_ <- initial.auto_resampling(data_)
  expect_equal(data_$recd_x[, "Process"], c(
    "Coordinate",
    "Coordinate",
    "Coordinate",
    "Coordinate->IQR",
    "Coordinate->IQR"
  ))
})

test_that("initial.auto_resampling: flattening trigger", {
  test_x <- data.frame(
    "1" = c(1, 2, 1, 2, 1, 2),
    "2" = c(1, 2, 1, 2, 1, 2),
    "3" = c(1, 2, 1, 2, 1, 2),
    "4" = c(1, 2, 1, 2, 1, 2),
    "5" = c(1, 0, 0, 0, 0, 0)
  )

  data_ <- rlang::env(
    syn_x = test_x,
    obs_x = test_x,
    obs_n = 5,
    syn_n = 10,
    x_deg = c(1, 1, 1, 1, 1),
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    )
  )

  data_ <- initial.auto_resampling(data_)
  expect_equal(data_$recd_x[, "Process"], c(
    "Coordinate",
    "Coordinate",
    "Coordinate",
    "Coordinate",
    "Coordinate->Flattening"
  ))
})

test_that("initial.auto_resampling: all trigger", {
  test_x <- data.frame(
    "1" = c(1, 2, 1, 2, 1, 2),
    "2" = c(1, 2, 3, 1, 2, 3),
    "3" = c(0, 0, 0, 0, 1, 100000),
    "4" = c(10000, 1, 0, 0, 0, 0),
    "5" = c(1, 0, 0, 0, 0, 0)
  )

  data_ <- rlang::env(
    syn_x = test_x,
    obs_x = test_x,
    syn_n = 10,
    obs_n = 5,
    x_deg = c(1, 3, 1, 3, 1),
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    )
  )

  data_ <- initial.auto_resampling(data_)
  expect_equal(data_$recd_x[, "Process"], c(
    "Coordinate",
    "Coordinate->Smoothing",
    "Coordinate->IQR",
    "Coordinate->Smoothing->IQR",
    "Coordinate->Flattening"
  ))
})
