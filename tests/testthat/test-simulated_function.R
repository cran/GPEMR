library(testthat)

test_that("simulated_function works for logistic model", {
  time_points <- 1:10
  n <- 10
  window_size <- 3
  model <- "Logistic"
  parameter <- list(r = 0.1, K = 50)
  sigma2 <- 2
  rho <- 0.5
  x_0 <- 10
  cov <- FALSE
  Plot_est <- FALSE

  result <- simulated_function(time_points, n, window_size, model, parameter, sigma2, rho, x_0, cov, Plot_est)

  expect_true(is.list(result))
  expect_named(result, c("data", "est_global", "cov_global", "est_local"))
  expect_true(all(sapply(result$data, is.numeric)))
})

test_that("simulated_function works for exponential model", {
  time_points <- 1:10
  n <- 10
  window_size <- 3
  model <- "Exponential"
  parameter <- list(r = 0.1)
  sigma2 <- 2
  rho <- 0.5
  x_0 <- 10
  cov <- FALSE
  Plot_est <- FALSE

  result <- simulated_function(time_points, n, window_size, model, parameter, sigma2, rho, x_0, cov, Plot_est)

  expect_true(is.list(result))
  expect_named(result, c("data", "est_global", "cov_global", "est_local"))
  expect_true(all(sapply(result$data, is.numeric)))
})
