library(testthat)

csv_path <- system.file("extdata", "sample_data.csv", package = "GPEMR")
txt_path <- system.file("extdata", "sample_data2.txt", package = "GPEMR")


# Test the comparison_plot function with the temporary CSV file
test_that("comparison_plot works with CSV data", {
  expect_silent(comparison_plot(data_path = csv_path,
                                window_size = 3,
                                parameter = list(log_r = 0.7, von_r = 0.2, b = 0.3, c = 0.1),
                                p_val_method = "Parametric",
                                repetition = 10))
})

# Test the comparison_plot function with the temporary TXT file
test_that("comparison_plot works with TXT data", {
  expect_silent(comparison_plot(data_path = txt_path,
                                window_size = 3,
                                parameter = list(log_r = 0.7, von_r = 0.2, b = 0.3, c = 0.1),
                                p_val_method = "Non-Parametric",
                                repetition = 10))
})
