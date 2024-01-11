test_that("memerised successfully", {
  expect_equal(
    generate_data_from_dist(from_dist = "gaussian"),
    generate_data_from_dist(from_dist = "gaussian")
  )
})
