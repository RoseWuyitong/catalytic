test_that("is_continuous: continuous", {
  expect_equal(is_continuous(c(1.1, 1.2, 1.1, 1.3, 0, 1)), TRUE)
})

test_that("is_continuous: category", {
  test_lst <- list(
    c(1.1, 0, 1.1, 0, 1.1, 1.1),
    as.factor(c("red", "blue", "green", "red", "blue")),
    c(TRUE, FALSE, TRUE, FALSE, FALSE),
    c(-1, 1, -1, 1, 1)
  )
  for (lst in test_lst) expect_equal(is_continuous(lst), FALSE)
})
