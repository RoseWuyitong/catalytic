test_that("memerised successfully", {
  expect_equal(
    resampling(c(1, 2, 3, 4, 5)),
    resampling(c(1, 2, 3, 4, 5))
  )
})

test_that("return NA for wrong resample method", {
  expect_equal(
    resampling(c(1, 2, 3, 4, 5), resample_method = "Not Exist"),
    NA
  )
})
