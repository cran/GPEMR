library(testthat)

csv_path <- system.file("extdata", "sample_data.csv", package = "GPEMR")
txt_path <- system.file("extdata", "sample_data2.txt", package = "GPEMR")


# Define parameters for testing
test_params <- list(r = 0.5, K = 100, c = 0.5, b = 0.2)

# Define a test for the real_data function
test_that("real_data works correctly for Logistic model with CSV input", {
  result <- real_data(
    data_path = csv_path,
    window_size = 3,
    model = "Logistic",
    parameter = test_params,
    cov = TRUE,
    Plot_est = FALSE,
    p_value_plot = FALSE,
    tolerance = 0.05
  )
  expect_type(result, "list")
  expect_true(all(c("est_local", "cov_local") %in% names(result)))
})

test_that("real_data works correctly for Logistic model with TXT input", {
  result <- real_data(
    data_path = txt_path,
    window_size = 3,
    model = "Logistic",
    parameter = test_params,
    cov = FALSE,
    Plot_est = FALSE,
    p_value_plot = FALSE,
    tolerance = 0.05
  )
  expect_type(result, "list")
  expect_true(all(c("est_local", "logistic_p_value") %in% names(result)))
})

test_that("real_data works correctly for Von-bertalanffy model with CSV input", {
  result <- real_data(
    data_path = csv_path,
    window_size = 3,
    model = "Von-bertalanffy",
    parameter = test_params,
    cov = TRUE,
    Plot_est = FALSE,
    p_value_plot = FALSE,
    tolerance = 0.05
  )
  expect_type(result, "list")
  expect_true(all(c("est_local", "cov_local", "VB_p_value") %in% names(result)))
})

test_that("real_data works correctly for Gompertz model with CSV input", {
  result <- real_data(
    data_path = csv_path,
    window_size = 3,
    model = "Gompertz",
    parameter = test_params,
    cov = FALSE,
    Plot_est = FALSE,
    p_value_plot = FALSE,
    tolerance = 0.05
  )
  expect_type(result, "list")
  expect_true(all(c("est_local", "Gom_p_value") %in% names(result)))
})
