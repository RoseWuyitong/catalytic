lmm_data <- mock_lmm_data()

cat_init <- cat_lmm_initialization(
  formula = Y ~ X1 + (1 | group),
  data = lmm_data,
  x_cols = c("X1", "X2", "X3"),
  y_col = "Y",
  z_cols = c("Z1", "Z2", "Z3"),
  group_col = NULL,
  syn_size = 100
)

# Basic functionality test
test_that("cat_lmm_initialization returns a list with expected structure and contents", {
  expected_components <- c("formula", "data", "x_cols", "y_col", "z_cols", "group_col",
                           "syn_size", "resample_by_group", "resample_only", "na_replace",
                           "simple_model", "obs_size", "obs_data", "obs_x", "obs_y",
                           "obs_z", "obs_group", "syn_size", "syn_data", "syn_x",
                           "syn_y", "syn_z", "syn_group")

  expect_type(cat_init, "list")
  expect_equal(cat_init$function_name, "cat_lmm_initialization")
  expect_true(all(expected_components %in% names(cat_init)))
})

# Test observation data still the same
test_that("Observation data has csame size and structuer", {
  expect_equal(cat_init$obs_size, 10)
  expect_equal(nrow(cat_init$obs_data), cat_init$obs_size)
  expect_equal(nrow(cat_init$obs_x), cat_init$obs_size)
  expect_equal(nrow(cat_init$obs_z), cat_init$obs_size)
  expect_equal(length(cat_init$obs_y), cat_init$obs_size)
  expect_equal(length(cat_init$obs_group), cat_init$obs_size)
  expect_true(all(colnames(cat_init$obs_x) %in% c("X1", "X2", "X3")))
  expect_true(all(colnames(cat_init$obs_z) %in% c("Z1", "Z2", "Z3")))
  expect_true(all(cat_init$obs_x == lmm_data[, c("X1", "X2", "X3")]))
  expect_true(all(cat_init$obs_z == lmm_data[, c("Z1", "Z2", "Z3")]))
  expect_true(all(cat_init$obs_y == lmm_data$Y))
  expect_true(all(cat_init$obs_group == lmm_data$group))
})

# Test synthetic data generation
test_that("Synthetic data has correct size and structure", {
  expect_equal(cat_init$syn_size, 100)
  expect_equal(nrow(cat_init$syn_data), cat_init$syn_size)
  expect_equal(nrow(cat_init$syn_x), cat_init$syn_size)
  expect_equal(nrow(cat_init$syn_z), cat_init$syn_size)
  expect_equal(length(cat_init$syn_y), cat_init$syn_size)
  expect_equal(length(cat_init$syn_group), cat_init$syn_size)
  expect_true(all(colnames(cat_init$syn_x) %in% c("X1", "X2", "X3")))
  expect_true(all(colnames(cat_init$syn_z) %in% c("Z1", "Z2", "Z3")))
})

# Test whole data combination
test_that("Synthetic data has correct size and structure", {
  expect_equal(cat_init$size, 110)
  expect_equal(cat_init$size, cat_init$obs_size + cat_init$syn_size)
  expect_equal(nrow(cat_init$data), cat_init$size)
  expect_equal(nrow(cat_init$x), cat_init$syn_size)
  expect_equal(nrow(cat_init$z), cat_init$syn_size)
  expect_equal(length(cat_init$y), cat_init$size)
  expect_equal(length(cat_init$group), cat_init$size)
  expect_true(all(colnames(cat_init$x) %in% c("X1", "X2", "X3")))
  expect_true(all(colnames(cat_init$z) %in% c("Z1", "Z2", "Z3")))
  expect_true(cat_init$y_col == "Y")
  expect_true(cat_init$group_col == "group")
})
