test_that("get_iqr: normal test", {
  expect_equal(get_iqr(c(100, 0, 0, 0)), 25)
})


test_that("get_iqr: limit test", {
  expect_equal(get_iqr(c(100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), 20)
})
