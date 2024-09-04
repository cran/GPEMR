#' Simulate and Estimate Parameters of Growth Models
#'
#' This function simulates data for independent trajectories based on specified growth models and estimates the model parameters.
#' The simulation can be performed for several models, including Logistic, Exponential, Theta-logistic, Von Bertalanffy, and Gompertz.
#' It also allows for calculation of global and local parameter estimates using the negative log-likelihood function.
#'
#' @param time_points A numeric vector representing the time points for which the data should be simulated.
#' @param n An integer specifying the number of independent trajectories to simulate.
#' @param window_size An integer specifying the window size for local estimation.
#' @param model A character string specifying the growth model to use. Options include 'Logistic', 'Exponential', 'Theta-logistic', 'Von-bertalanffy', and 'Gompertz'.
#' @param parameter A list of model-specific parameters required for the mean function.
#' @param sigma2 A numeric value for the variance of the process.
#' @param rho A numeric value for the correlation coefficient.
#' @param x_0 A numeric value for the initial state.
#' @param cov A logical value indicating whether to print the covariance matrix. Default is FALSE.
#' @param Plot_est A logical value indicating whether to plot the parameter estimates. Default is FALSE.
#' @return A list containing the simulated data, global parameter estimates, global covariance matrix, local parameter estimates, and optionally local covariance matrices.
#' @details
#' The function first checks if the parameters are provided as a list. It then calculates the mean function based on the specified model and forms the covariance matrix. Multivariate normal data for the specified number of trajectories is generated using the mvtnorm::rmvnorm function. The negative log-likelihood function is defined and minimized using the optim function to estimate global parameters. Local parameter estimation is performed using a sliding window approach.
#'
#' The available models are:
#' - Logistic: Requires parameters r (growth rate) and K (carrying capacity).
#' - Exponential: Requires parameter r (growth rate).
#' - Theta-logistic: Requires parameters r (growth rate), theta, and K (carrying capacity).
#' - Von-bertalanffy: Requires parameters r (growth rate) and K (asymptotic size).
#' - Gompertz: Requires parameters b and c.
#'
#' @examples
#'res <- simulated_function(
#'  model = 'Logistic',
#'  parameter = list(r = 0.2, K = 100),
#'  cov = TRUE,
#'  Plot_est = TRUE)
#'

#' @importFrom stats optim sd pnorm
#' @importFrom graphics par abline legend lines points segments
#' @importFrom utils read.table read.csv
#' @export


simulated_function <- function(
    time_points = 1:10,
    n = 10,
    window_size = 3,
    model,
    parameter,
    sigma2 = 2,
    rho = 0.5,
    x_0 = 10,
    cov = FALSE,
    Plot_est = FALSE
) {

  if (!is.list(parameter)) stop("parameter must be a list.")

  t <- time_points
  q <- length(t)
  P <- window_size
  h <- numeric(length = length(t) - 1)
  for (i in 1:(length(t) - 1)) {
    h[i] <- t[i + 1] - t[i]
  }

  # Calculate mean function based on the model
  if (model == 'Logistic') {
    r <- parameter$r
    K <- parameter$K
    mu <- K / (1 + (K / x_0 - 1) * exp(-r * t))   # Logistic mean function
  } else if (model == 'Exponential') {
    r <- parameter$r
    mu <- x_0 * exp(r * t) # Exponential mean function
  } else if (model == 'Theta-logistic') {
    r <- parameter$r
    theta <- parameter$theta
    K <- parameter$K
    mu <- K / (1 + ((K / x_0)^theta - 1) * exp(-r * theta * t))
  } else if (model == 'Von-bertalanffy') {
    r <- parameter$r
    K <- parameter$K
    mu <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
  } else if (model == 'Gompertz') {
    b <- parameter$b
    c <- parameter$c
    mu <- x_0 * exp((b / c) * (1 - exp(-c * t)))  # Gompertz mean function
  } else {
    stop("Invalid model specified.")
  }

  # Formation of covariance matrix
  cov_mat <- matrix(data = NA, nrow = q, ncol = q)
  for (i in 1:q) {
    for (j in 1:q) {
      cov_mat[i, j] <- sigma2 * rho^(abs(i - j))
    }
  }

  # Generation of multivariate normal data for n trajectories
  data_l <- mvtnorm::rmvnorm(n, mean = mu, sigma = cov_mat)

  # Calculation of negative log-likelihood
  fun_likelihood <- function(param) {
    if (model == 'Logistic') {
      r <- param[1]
      K <- param[2]
      sigma2 <- param[3]
      rho <- param[4]
      mu <- K / (1 + (K / x_0 - 1) * exp(-r * t))
    } else if (model == 'Exponential') {
      r <- param[1]
      sigma2 <- param[2]
      rho <- param[3]
      mu <- x_0 * exp(r * t)
    } else if (model == 'Theta-logistic') {
      r <- param[1]
      K <- param[2]
      theta <- param[3]
      sigma2 <- param[4]
      rho <- param[5]
      mu <- K / (1 + ((K / x_0)^theta - 1) * exp(-r * theta * t))
    } else if (model == 'Von-bertalanffy') {
      r <- param[1]
      K <- param[2]
      sigma2 <- param[3]
      rho <- param[4]
      mu <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
    } else if (model == 'Gompertz') {
      c <- param[1]
      b <- param[2]
      sigma2 <- param[3]
      rho <- param[4]
      mu <- x_0 * exp((b / c) * (1 - exp(-c * t)))
    }

    cov_mat <- matrix(data = NA, nrow = q, ncol = q)
    for (i in 1:q) {
      for (j in 1:q) {
        cov_mat[i, j] <- sigma2 * rho^(abs(i - j))
      }
    }

    exponent <- 0
    for (i in 1:nrow(data_l)) {
      exponent <- exponent + t(data_l[i, ] - mu) %*% solve(cov_mat) %*% (data_l[i, ] - mu)
    }

    log_lik <- - (nrow(data_l) * q / 2) * log(2 * pi) - (nrow(data_l) / 2) * log(det(cov_mat)) - (1 / 2) * exponent
    return(-log_lik)
  }

  if (model == 'Logistic') {
    param_0 <- c(0.3, parameter$K, sigma2, rho)
    lower <- c(0.001, parameter$K - (parameter$K * 0.2), 0.0001, -0.001)
    upper <- c(3, parameter$K + (parameter$K * 0.2), 10, 0.999)
  } else if (model == 'Exponential') {
    param_0 <- c(parameter$r, sigma2, rho)
    lower <- c(0.001, 0.0001, -0.001)
    upper <- c(3, 10, 0.999)
  } else if (model == 'Theta-logistic') {
    param_0 <- c(0.3, parameter$K, parameter$theta, sigma2, rho)
    lower <- c(0.001, parameter$K - (parameter$K * 0.2), 0.1, 0.0001, -0.001)
    upper <- c(3, parameter$K + (parameter$K * 0.2), 2, 10, 0.999)
  } else if (model == 'Von-bertalanffy') {
    param_0 <- c(0.8, parameter$K, sigma2, rho)
    lower <- c(0.001, parameter$K - (parameter$K * 0.2), 0.0001, -0.001)
    upper <- c(3, parameter$K + (parameter$K * 0.2), 10, 0.999)
  } else if (model == 'Gompertz') {
    param_0 <- c(parameter$c, parameter$b, sigma2, rho)
    lower <- c(0.05, 0.05, 0.005, -0.001)
    upper <- c(2, 2, 2, 0.999)
  }

  out <- optim(param_0, fun_likelihood, method = "L-BFGS-B", hessian = TRUE, lower = lower, upper = upper)
  global_estimate <- out$par
  cov_mat_global <- solve(out$hessian)

  if (model == 'Logistic' || model == 'Von-bertalanffy' || model == 'Gompertz') {
    colnames(cov_mat_global) <- c("r", "K", "sigma2", "rho")
  } else if (model == 'Exponential') {
    colnames(cov_mat_global) <- c("r", "sigma2", "rho")
  } else if (model == 'Theta-logistic') {
    colnames(cov_mat_global) <- c("r", "K", "theta", "sigma2", "rho")
  }else if(model == 'Gompertz'){
    colnames(cov_mat_global) <- c("b", "v", "sigma2", "rho")
  }

  fun_local_likelihood <- function(loc_parameter) {
    # Extract model-specific parameters
    if (model == 'Logistic') {
      r <- loc_parameter[1]
      K <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu_loc <- K / (1 + (K / x_0 - 1) * exp(-r * t))
    } else if (model == 'Exponential') {
      r <- loc_parameter[1]
      sigma2 <- loc_parameter[2]
      rho <- loc_parameter[3]
      mu_loc <- x_0 * exp(r * t)
    } else if (model == 'Theta-logistic') {
      r <- loc_parameter[1]
      K <- loc_parameter[2]
      theta <- loc_parameter[3]
      sigma2 <- loc_parameter[4]
      rho <- loc_parameter[5]
      mu_loc <- K / (1 + ((K / x_0)^theta - 1) * exp(-r * theta * t))
    } else if (model == 'Von-bertalanffy') {
      r <- loc_parameter[1]
      K <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu_loc <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
    } else if (model == 'Gompertz') {
      b <- loc_parameter[1]
      c <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu_loc <- x_0 * exp((b / c) * (1 - exp(-c * t)))
    } else {
      stop("Invalid model specified.")
    }

    # Ensure indices are within bounds
    if (ind + P - 1 > length(t)) {
      stop("Window size is too large for the remaining data points.")
    }
    if (ind + P - 1 > ncol(data_l)) {
      stop("Window size is too large for the data dimensions.")
    }

    # Create covariance matrix
    cov_mat_ = matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:q) {
      for (j in 1:q) {
        cov_mat_[i,j] = sigma2*rho^(abs(i-j))
      }
    }

    # localized estimates
    mu_local = mu_loc[ind:(ind + P-1)]
    data_local = data_l[,ind:(ind + P - 1)]
    cov_mat_local = cov_mat_[ind:(ind+P-1), ind:(ind+P-1)] # Ensure proper dimension for covariance matrix

    # Check if covariance matrix dimensions match
    if (nrow(cov_mat_local) != P || ncol(cov_mat_local) != P) {
      stop("Covariance matrix dimensions do not match the expected size.")
    }

    exponent <- 0
    for (i in 1:n) {
      exponent <- exponent + t(data_local[i, ] - mu_local) %*% solve(cov_mat_local) %*% (data_local[i, ] - mu_local)
    }

    log_lik <- - (n * length(mu_local) / 2) * log(2 * pi) - (n / 2) * log(det(cov_mat_local)) - (1 / 2) * exponent
    return(-log_lik)
  }

  # Determine the number of columns and column names for est_local based on the model
  if (model == 'Logistic' || model == 'Von-bertalanffy') {
    est_local <- matrix(data = NA, nrow = q - P + 1, ncol = 4)
    cov_local=list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("r", "K", "sigma2", "rho")
  } else if (model == 'Exponential') {
    est_local <- matrix(data = NA, nrow = q - P + 1, ncol = 3)
    cov_local=list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("r", "sigma2", "rho")
  } else if (model == 'Theta-logistic') {
    est_local <- matrix(data = NA, nrow = q - P + 1, ncol = 5)
    cov_local=list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("r", "K", "theta", "sigma2", "rho")
  }else if(model == 'Gompertz'){
    est_local <- matrix(data = NA, nrow = q - P + 1, ncol = 4)
    cov_local=list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("b", "c", "sigma2", "rho")
  }

  # ISRP estimates by maximizing the localized likelihood function
  for (j in 1:(q - P + 1)) {
    ind <- j
    if (model == 'Logistic') {
      param_0 <- c(0.3, parameter$K, sigma2, rho)
      lower <- c(0.001, parameter$K - (parameter$K * 0.2), 0.0001, -0.001)
      upper <- c(3, parameter$K + (parameter$K * 0.2), 10, 0.999)
    } else if (model == 'Exponential') {
      param_0 <- c(parameter$r, sigma2, rho)
      lower <- c(0.001, 0.0001, -0.001)
      upper <- c(3, 10, 0.999)
    } else if (model == 'Theta-logistic') {
      param_0 <- c(0.3, parameter$K, parameter$theta, sigma2, rho)
      lower <- c(0.001, parameter$K - (parameter$K * 0.2), 0.1, 0.0001, -0.001)
      upper <- c(3, parameter$K + (parameter$K * 0.2), 2, 10, 0.999)
    } else if (model == 'Von-bertalanffy') {
      param_0 <- c(0.8, parameter$K, sigma2, rho)
      lower <- c(0.001, parameter$K - (parameter$K * 0.2), 0.0001, -0.001)
      upper <- c(3, parameter$K + (parameter$K * 0.2), 10, 0.999)
    } else if (model == 'Gompertz') {
      param_0 <- c(parameter$b, parameter$c, sigma2, rho)
      lower <- c(0.05, 0.05, 0.005, -0.001)
      upper <- c(2, 2, 2, 0.999)
    }

    opt_results <- optim(
      par = param_0,
      fn = fun_local_likelihood,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      hessian = TRUE
    )

    cov_est_local <- solve(opt_results$hessian)
    cov_local[[j]] <- cov_est_local
    est_local[ind, ] <- opt_results$par

  }

  if (cov==TRUE) {
    res <- list(data = data_l,est_global = global_estimate,cov_global = cov_mat_global,est_local = est_local,cov_local=cov_local)
  }else{
    res <- list(data = data_l,est_global = global_estimate,cov_global = cov_mat_global,est_local = est_local)
  }

  est_plot <- if(Plot_est == TRUE){
    if(model == 'Logistic'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(1:(q-P+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(r(Delta~t))[MLE]),
           ylab = expression(hat(r)))
      plot(1:(q-P+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(K(Delta~t))[MLE]),
           ylab = expression(hat(K)))
      plot(1:(q-P+1), est_local[,3], type = "b", col = "red", lwd = 2,
           main = expression(widehat(sigma^2)),
           ylab = expression(hat(sigma^2)))
      plot(1:(q-P+1), est_local[,4], type = "b", col = "red", lwd = 2,
           main = expression(widehat(rho)),
           ylab = expression(hat(rho)))
    }else if(model == 'Exponential'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(1:(q-P+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(r(Delta~t))[MLE]),
           ylab = expression(hat(r)))
      plot(1:(q-P+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(sigma^2)),
           ylab = expression(hat(sigma^2)))
      plot(1:(q-P+1), est_local[,3], type = "b", col = "red", lwd = 2,
           main = expression(widehat(rho)),
           ylab = expression(hat(rho)))
    }else if(model == 'Theta-logistic'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(3,2))
      plot(1:(q-P+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(r(Delta~t))[MLE]),
           ylab = expression(hat(r)))
      plot(1:(q-P+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(K(Delta~t))[MLE]),
           ylab = expression(hat(K)))
      plot(1:(q-P+1), est_local[,3], type = "b", col = "red", lwd = 2,
           main = expression(widehat(theta(Delta~t))[MLE]),
           ylab = expression(hat(theta)))
      plot(1:(q-P+1), est_local[,4], type = "b", col = "red", lwd = 2,
           main = expression(widehat(sigma^2)),
           ylab = expression(hat(sigma^2)))
      plot(1:(q-P+1), est_local[,5], type = "b", col = "red", lwd = 2,
           main = expression(widehat(rho)),
           ylab = expression(hat(rho)))
    }else if(model == 'Von-bertalanffy'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(1:(q-P+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(r(Delta~t))[MLE]),
           ylab = expression(hat(r)))
      plot(1:(q-P+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(K(Delta~t))[MLE]),
           ylab = expression(hat(K)))
      plot(1:(q-P+1), est_local[,3], type = "b", col = "red", lwd = 2,
           main = expression(widehat(sigma^2)),
           ylab = expression(hat(sigma^2)))
      plot(1:(q-P+1), est_local[,4], type = "b", col = "red", lwd = 2,
           main = expression(widehat(rho)),
           ylab = expression(hat(rho)))
    }else if(model == 'Gompertz'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(1:(q-P+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(c(Delta~t))[MLE]),
           ylab = expression(hat(c)))
      plot(1:(q-P+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(b(Delta~t))[MLE]),
           ylab = expression(hat(b)))
      plot(1:(q-P+1), est_local[,3], type = "b", col = "red", lwd = 2,
           main = expression(widehat(sigma^2)),
           ylab = expression(hat(sigma^2)))
      plot(1:(q-P+1), est_local[,4], type = "b", col = "red", lwd = 2,
           main = expression(widehat(rho)),
           ylab = expression(hat(rho)))

    }
  }
  return(c(res,est_plot))

}



#' Real Data Analysis with Model Fitting and Visualization
#'
#' This function performs parameter estimation for specified models (Logistic, Von-Bertalanffy, or Gompertz) using real data from a TXT or PDF file. It calculates global and local estimates of the model parameters, including their covariance matrices, and optionally generates plots of the estimates and p-values.
#'
#' @param data_path A character string specifying the path to the data file. Supported file types are TXT and CSV.
#' @param window_size An integer specifying the size of the moving window for localized estimation. Default is 3.
#' @param model A character string specifying the model to fit. Options are "Logistic", "Von-Bertalanffy", and "Gompertz".
#' @param parameter A list containing initial parameter estimates. For "Logistic" and "Von-Bertalanffy" models, the list should include `r` and `K`. For the "Gompertz" model, it should include `c` and `b`.
#' @param cov A logical value indicating whether to return the covariance matrices. Default is FALSE.
#' @param Plot_est A logical value indicating whether to generate plots of the estimated parameters. Default is FALSE.
#' @param p_value_plot A logical value indicating whether to generate plots of p-values for local estimates. Default is FALSE.
#' @param tolerance A numeric value specifying the alpha level for p-value calculation. Default is 0.05.
#'
#' @return A list containing:
#' \item{est_global}{A vector of global parameter estimates.}
#' \item{cov_global}{A matrix of global covariance estimates (if `cov` is TRUE).}
#' \item{est_local}{A matrix of local parameter estimates.}
#' \item{cov_local}{A list of local covariance matrices (if `cov` is TRUE).}
#' \item{est_plot}{A plot of the estimated parameters (if `Plot_est` is TRUE).}
#' \item{p_value_plot}{A plot of p-values for local estimates (if `p_value_plot` is TRUE).}
#'
#' @examples
#' data_csv <- system.file("extdata", "sample_data.csv", package = "GPEMR")
#' results_logistic <- real_data(data_csv, window_size = 5, model = "Logistic",
#'                               parameter = list(r= 0.7), cov = TRUE, Plot_est = TRUE)
#'
#' @importFrom stats optim sd pnorm
#' @importFrom graphics par abline legend lines points segments
#' @importFrom utils read.table read.csv
#'
#' @export


real_data <- function(
    data_path,
    window_size,
    model,
    parameter,
    cov = FALSE,
    Plot_est = FALSE,
    p_value_plot = FALSE,
    tolerance = 0.05) {

  if (!is.list(parameter)) stop("parameter must be a list.")

  ext <- tools::file_ext(data_path)

  # Read data based on file type
  if (ext == "txt") {
    data_r <- as.matrix(read.table(data_path, header = FALSE))
  } else if (ext == "csv") {
    data_r <- as.matrix(read.csv(data_path))
  } else {
    stop("Unsupported file type. Only TXT and CSV files are supported.")
  }

  P <- window_size

  row_number=nrow(data_r)
  col_number=ncol(data_r)
  q=ncol(data_r)
  t=1:q


  initial=mean(data_r[,1])
  final=mean(data_r[,q])

  K1<-final-final*0.05
  K2<-final+final*0.05
  x_0<- initial
  K<- final

  sigma2= 1 # We will fix it for al  the model
  rho= 0.5 # We will fix it for all the

  n<- row_number

  # Calculation of negative log-likelihood
  fun_likelihood <- function(param) {
    if (model == 'Logistic') {
      r <- param[1]
      K <- param[2]
      sigma2 <- param[3]
      rho <- param[4]
      mu <- K / (1 + (K / x_0 - 1) * exp(-r * t))
    } else if (model == 'Von-bertalanffy') {
      r <- param[1]
      K <- param[2]
      sigma2 <- param[3]
      rho <- param[4]
      mu <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
    } else if (model == 'Gompertz') {
      c <- param[1]
      b <- param[2]
      sigma2 <- param[3]
      rho <- param[4]
      mu <- x_0 * exp((b / c) * (1 - exp(-c * t)))
    }

    cov_mat <- matrix(data = NA, nrow = q, ncol = q)
    for (i in 1:q) {
      for (j in 1:q) {
        cov_mat[i, j] <- sigma2 * rho^(abs(i - j))
      }
    }

    exponent <- 0
    for (i in 1:nrow(data_r)) {
      exponent <- exponent + t(data_r[i, ] - mu) %*% solve(cov_mat) %*% (data_r[i, ] - mu)
    }

    log_lik <- - (nrow(data_r) * q / 2) * log(2 * pi) - (nrow(data_r) / 2) * log(det(cov_mat)) - (1 / 2) * exponent
    return(-log_lik)
  }

  if (model == 'Logistic') {
    param_0 <- c(parameter$r, K, 3,0.3)
    lower <- c(0.001, K1, 0.0001, -0.001)
    upper <- c(3, K2, 3, 0.999)
  } else if (model == 'Von-bertalanffy') {
    param_0 <- c(parameter$r, K, 2,0.3)
    lower <- c(0.001, K1, 0.0001, -0.001)
    upper <- c(3, K2, 3, 0.999)
  } else if (model == 'Gompertz') {
    param_0 <- c(parameter$c, parameter$b, 2,0.3)
    lower <- c(0.01, 0.01, 0.001, -0.001)
    upper <- c(5, 10, 3, 0.999)
  }

  out_global <- optim(param_0, fun_likelihood, method = "L-BFGS-B", hessian = TRUE, lower = lower, upper = upper)
  global_estimate <- out_global$par
  cov_mat_global <- solve(out_global$hessian)

  if (model == 'Logistic' || model == 'Von-bertalanffy') {
    colnames(cov_mat_global) <- c("r", "K", "sigma2", "rho")
  }else if(model == 'Gompertz'){
    colnames(cov_mat_global) <- c("b", "v", "sigma2", "rho")
  }

  fun_local_likelihood <- function(loc_parameter) {
    # Extract model-specific parameters
    if (model == 'Logistic') {
      r <- loc_parameter[1]
      K <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu <- K / (1 + (K / x_0 - 1) * exp(-r * t))
    }  else if (model == 'Von-bertalanffy') {
      r <- loc_parameter[1]
      K <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
    } else if (model == 'Gompertz') {
      b <- loc_parameter[1]
      c <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu <- x_0 * exp((b / c) * (1 - exp(-c * t)))
    } else {
      stop("Invalid model specified.")
    }

    # Ensure indices are within bounds
    if (ind + window_size - 1 > length(t)) {
      stop("Window size is too large for the remaining data points.")
    }
    if (ind + window_size - 1 > ncol(data_r)) {
      stop("Window size is too large for the data dimensions.")
    }

    # Create covariance matrix
    cov_mat_ = matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:q) {
      for (j in 1:q) {
        cov_mat_[i,j] = sigma2*rho^(abs(i-j))
      }
    }

    # localized estimates
    mu_local = mu[ind:(ind + window_size-1)]
    data_local = data_r[,ind:(ind + window_size - 1)]
    cov_mat_local = cov_mat_[ind:(ind+window_size-1), ind:(ind+window_size-1)] # Ensure proper dimension for covariance matrix

    # Check if covariance matrix dimensions match
    if (nrow(cov_mat_local) != window_size || ncol(cov_mat_local) != window_size) {
      stop("Covariance matrix dimensions do not match the expected size.")
    }

    exponent <- 0
    for (i in 1:n) {
      exponent <- exponent + t(data_local[i, ] - mu_local) %*% solve(cov_mat_local) %*% (data_local[i, ] - mu_local)
    }

    log_lik <- - (n * length(mu_local) / 2) * log(2 * pi) - (n / 2) * log(det(cov_mat_local)) - (1 / 2) * exponent
    return(-log_lik)
  }

  # Determine the number of columns and column names for est_local based on the model
  if (model == 'Logistic' || model == 'Von-bertalanffy' ) {
    est_local <- matrix(data = NA, nrow = q - window_size + 1, ncol = 4)
    cov_local<-list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("r", "K", "sigma2", "rho")
  }else if(model == 'Gompertz'){
    est_local <- matrix(data = NA, nrow = q - window_size + 1, ncol = 4)
    cov_local<-list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("b", "v", "sigma2", "rho")
  }

  # ISRP estimates by maximizing the localized likelihood function
  for (j in 1:(q - window_size + 1)) {
    ind <- j
    if (model == 'Logistic') {
      param_0 <- c(parameter$r,K,3,0.3)
      lower <- c(0.001, K1, 0.0001, -0.001)
      upper <- c(3, K2, 3, 0.999)
    } else if (model == 'Von-bertalanffy') {
      param_0 <- c(parameter$r,K,3,0.3)
      lower <- c(0.001, K1, 0.0001, -0.001)
      upper <- c(3, K2, 3, 0.999)
    } else if (model == 'Gompertz') {
      param_0 <- c(parameter$b, parameter$c, 2, 0.3)
      lower <- c(0.01, 0.01, 0.001, -0.001)
      upper <- c(5, 10, 3, 0.999)
    }

    opt_results <- optim(
      par = param_0,
      fn = fun_local_likelihood,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      hessian = TRUE
    )

    cov_est_local <- solve(opt_results$hessian)
    cov_local[[j]] <- cov_est_local
    est_local[j, ] <- opt_results$par

  }

  # window_size VALUE

  if (model == 'Logistic'){
    RGR_hat_log <- numeric(length = q-window_size+1)
    for(i in 1:(q-window_size+1)){
      RGR_hat_log[i] <- est_local[i,1]*(1- mean(data_r[,i])/est_local[i,2])
    }

    var_RGR_hat_log <- numeric(length = q-window_size+1)
    for(i in 1:(q-P+1)){
      A = cov_local[[i]][1,1]*(1-mean(data_r[,i])/est_local[i,2])^2
      B = cov_local[[i]][2,2]*(est_local[i,1]*mean(data_r[,i])/est_local[i,2]^2)^2
      C = 2*(1-mean(data_r[,i])/est_local[i,2])*(est_local[i,1]*mean(data_r[,i])/est_local[i,2]^2)*cov_local[[i]][1,2]
      var_RGR_hat_log[i] = A+B+C
    }

    V_log <- numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      V_log[i] = log(1/RGR_hat_log[i] - 1/est_local[i,1])
    }

    var_V_log <- numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      A = (-est_local[i,1]/(RGR_hat_log[i]*(est_local[i,1]-RGR_hat_log[i])))^2*var_RGR_hat_log[i]
      B = (RGR_hat_log[i]/(est_local[i,1]*(est_local[i,1]-RGR_hat_log[i])))^2*cov_local[[i]][1,1]
      C = 2*(-1/(est_local[i,1]-RGR_hat_log[i])^2)*(cov_local[[i]][1,1]*(1-mean(data_r[,i])/est_local[i,2]) + (est_local[i,1]*mean(data_r[,i])/(est_local[i,2])^2)*cov_local[[i]][1,2])
      var_V_log[i] = A + B + C
    }

    test_stat_log <- numeric(length = q-window_size-1)
    var_test_stat_log <- numeric(length = q-window_size-1)
    for(i in 1:(q-window_size-1)){
      var_test_stat_log[i] = var_V_log[i+2] + 4*var_V_log[i+1] + var_V_log[i]
      test_stat_log[i] = (V_log[i+2]-2*V_log[i+1] + V_log[i])
    }

    p_value_local_log_l <- numeric(q-window_size-1)
    for(i in 1:(q-window_size-1)){
      p_value_local_log_l[i] = 1- pnorm(abs(test_stat_log[i]/sqrt(var_test_stat_log[i])))
    }


  }else if(model == 'Von-bertalanffy'){
    RGR_hat_VB = numeric(length = q-window_size+1)
    for(i in 1:(q-window_size+1)){
      RGR_hat_VB[i] = est_local[i,1]*((est_local[i,2]/mean(data_r[,i]))^(1/3))*(1- (mean(data_r[,i])/est_local[i,2])^(1/3))
    }

    var_RGR_hat_VB <- numeric(length = q-window_size+1)
    for(i in 1:(q-P+1)){
      A = cov_local[[i]][1,1]*((est_local[i,2]/mean(data_r[,i]))^(1/3)-1)^2
      B = cov_local[[i]][2,2]*(est_local[i,1]/(3*(mean(data_r[,i])*est_local[i,2]^2)^(1/3)))^2
      C = 2*((est_local[i,2]/mean(data_r[,i]))^(1/3)-1)*(est_local[i,1]/(3*(mean(data_r[,i])*est_local[i,2]^2)^(1/3)))*cov_local[[i]][1,2]
      var_RGR_hat_VB[i] = A+B+C
    }
    V_VB <- numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      V_VB[i] = log(1/RGR_hat_VB[i] + 1/est_local[i,1])
    }

    var_V_VB <- numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      A = (-est_local[i,1]/(RGR_hat_VB[i]*(est_local[i,1]+RGR_hat_VB[i])))^2*var_RGR_hat_VB[i]
      B = (-RGR_hat_VB[i]/(est_local[i,1]*(est_local[i,1]+RGR_hat_VB[i])))^2*cov_local[[i]][1,1]
      C = 2*(1/(est_local[i,1]+RGR_hat_VB[i])^2)*(cov_local[[i]][1,1]*(1-mean(data_r[,i])/est_local[i,2]) + (est_local[i,1]*mean(data_r[,i])/(est_local[i,2])^2)*cov_local[[i]][1,2])
      var_V_VB[i] = A + B + C
    }

    test_stat_VB <- numeric(length = q-window_size-1)
    var_test_stat_VB = numeric(length = q-window_size-1)
    for(i in 1:(q-window_size-1)){
      var_test_stat_VB[i] = var_V_VB[i+2] + 4*var_V_VB[i+1] + var_V_VB[i]
      test_stat_VB[i] = (V_VB[i+2]-2*V_VB[i+1] + V_VB[i])
    }

    # Computation of the p-value
    p_value_local_VB <- numeric(q-window_size-1)
    for(i in 1:(q-window_size-1)){
      p_value_local_VB[i] = 1- pnorm(test_stat_VB[i]/sqrt(var_test_stat_VB[i]))
    }

  }else if (model == 'Gompertz'){
    RGR_hat_gom <- numeric(length = q-window_size+1)
    for(i in 1:(q-window_size+1)){
      RGR_hat_gom[i] = est_local[i,2]*exp(-est_local[i,1]*t[i])
    }

    var_RGR_hat_gom <- numeric(length = q-window_size+1)
    for(i in 1:(q-P+1)){
      A = cov_local[[i]][1,1]*(-est_local[i,2]*t[i]*exp(-est_local[i,1]*t[i]))^2
      B = cov_local[[i]][2,2]*(exp(-est_local[i,1]*t[i]))^2
      C = 2*(-est_local[i,2]*t[i]*exp(-est_local[i,1]*t[i]))*(exp(-est_local[i,1]*t[i]))*cov_local[[i]][1,2]
      var_RGR_hat_gom[i] = A+B+C
    }

    # Computation of V
    V_gom = numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      V_gom[i] = log(1/RGR_hat_gom[i])
    }

    var_V_gom <- numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      A = cov_local[[i]][1,1]*(-t[i])^2
      B = cov_local[[i]][2,2]*(1/est_local[i,2])^2
      C = 2*(-t[i]/est_local[i,2])*cov_local[[i]][1,2]
      var_V_gom[i] = A + B + C
    }


    # Computation of the test statistic
    test_stat_gom <- numeric(length = q-window_size-1)
    var_test_stat_gom <- numeric(length = q-window_size-1)
    for(i in 1:(q-window_size-1)){
      var_test_stat_gom[i] = var_V_gom[i+2] + 4*var_V_gom[i+1] + var_V_gom[i]
      test_stat_gom[i] = (V_gom[i+2]-2*V_gom[i+1] + V_gom[i])
    }

    p_value_local_gom <- numeric(q-window_size-1)
    for(i in 1:(q-window_size-1)){
      p_value_local_gom[i] = 1- pnorm(abs(test_stat_gom[i]/sqrt(var_test_stat_gom[i])))
    }
  }

  if(model == 'Logistic'){
    p_val <- list(logistic_p_value = p_value_local_log_l)


  }else if(model == 'Von-bertalanffy'){
    p_val <- list(VB_p_value = p_value_local_VB)

  }else if(model == 'Gompertz'){
    p_val <- list(Gom_p_value = p_value_local_gom)

  }

  if (cov==TRUE) {
    res <- list(est_global = global_estimate,cov_global = cov_mat_global,est_local = est_local,cov_local=cov_local)
  }else{
    res <- list(est_global = global_estimate,cov_global = cov_mat_global,est_local = est_local)
  }

  est_plot <- if(Plot_est == TRUE){
    if(model == 'Logistic'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(1:(q-window_size+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(r(Delta~t))[MLE]),
           ylab = expression(hat(r)))
      plot(1:(q-window_size+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(K(Delta~t))[MLE]),
           ylab = expression(hat(K)))}
    else if(model == 'Von-bertalanffy'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(1:(q-window_size+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(r(Delta~t))[MLE]),
           ylab = expression(hat(r)))
      plot(1:(q-window_size+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(K(Delta~t))[MLE]),
           ylab = expression(hat(K)))
    }else if(model == 'Gompertz'){
      oldpar <- par(no.readonly = TRUE)  # Save current par settings
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(1:(q-window_size+1), est_local[,1], type = "b", col = "red", lwd = 2,
           main = expression(widehat(c(Delta~t))[MLE]),
           ylab = expression(hat(c)))
      plot(1:(q-window_size+1), est_local[,2], type = "b", col = "red", lwd = 2,
           main = expression(widehat(b(Delta~t))[MLE]),
           ylab = expression(hat(b)))}

  }

  P_value_plot <- if(p_value_plot == TRUE){
    if(model == 'Logistic'){
      plot(1:(q-window_size-1), p_value_local_log_l, type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'window_size Value Plot',xlab = 'Time Points')
      points(1:(q-window_size-1), p_value_local_log_l, col = "red", pch = 19)
      segments(1:(q-window_size-1), 0, 1:(q-window_size-1), p_value_local_log_l, col = "black", lwd = 2)
      abline(h = tolerance, lwd = 2, col = "black", lty = 2)
    } else if(model == 'Von-bertalanffy'){
      plot(1:(q-window_size-1), p_value_local_VB, type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'window_size Value Plot')
      points(1:(q-window_size-1), p_value_local_VB, col = "red", pch = 19)
      segments(1:(q-window_size-1), 0, 1:(q-window_size-1), p_value_local_VB, col = "black", lwd = 2)
      abline(h = tolerance, lwd = 2, col = "black", lty = 2)
    } else if(model == 'Gompertz'){
      plot(1:(q-window_size-1), p_value_local_gom, type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'window_size Value Plot',xlab = 'Time Points')
      points(1:(q-window_size-1), p_value_local_gom, col = "red", pch = 19)
      segments(1:(q-window_size-1), 0, 1:(q-window_size-1), p_value_local_gom, col = "black", lwd = 2)
      abline(h = tolerance, lwd = 2, col = "black", lty = 2)
    }
  }
  return(c(res,est_plot,P_value_plot,p_val))
}



#' Comparison Plot for p-values Across Growth Models
#'
#' This function generates a comparison plot of p-values obtained from different growth models (Logistic, Von-Bertalanffy, and Gompertz) using real data. It visualizes how p-values vary across time points for the specified models and provides a horizontal reference line for significance.
#'
#' @param data_path A character string specifying the path to the data file. The file should contain the necessary data for p-value calculations.
#' @param window_size An integer specifying the size of the moving window for calculations. Default is 3.
#' @param parameter A list of parameters used in p-value calculation.
#' @param p_val_method A character string specifying the method for p-value calculation. The options are Parametric and Non - Parametric methods.
#' @param repetition An integer specifying the number of repetitions for calculations.
#'
#' @return A plot comparing p-values across the Logistic, Von-Bertalanffy, and Gompertz models over time points. The plot includes lines for each model's p-values, a legend, and a reference line of tolerance = 0.05.
#'
#' @examples
#'\donttest{
#' data_path <- system.file("extdata", "sample_data.csv", package = "GPEMR")
#' params <- list(log_r = 0.7,von_r = 0.3,c=0.2,b=0.1)
#' comparison_plot(data_path = data_path, window_size = 5, parameter = params,
#'                 p_val_method = "Parametric", repetition = 2)}
#'
#' @importFrom stats optim sd pnorm
#' @importFrom graphics par abline legend lines points segments
#' @importFrom utils read.table read.csv
#' @export



comparison_plot <- function(data_path,
                            window_size = 3,
                            parameter,
                            p_val_method,
                            repetition){

  ext <- tools::file_ext(data_path)

  # Read data based on file type
  if (ext == "txt") {
    data_r <- as.matrix(read.table(data_path, header = FALSE))
  } else if (ext == "csv") {
    data_r <- as.matrix(read.csv(data_path))
  } else {
    stop("Unsupported file type. Only TXT and CSV files are supported.")
  }

  P<- window_size
  q <- ncol(data_r)


  log_p_val <- p_value_calc(data_path,window_size,parameter,p_val_method,repetition,model='Logistic')
  von_p_val <- p_value_calc(data_path,window_size,parameter,p_val_method,repetition,model='Von-bertalanffy')
  gom_p_val <- p_value_calc(data_path,window_size,parameter,p_val_method,repetition,model='Gompertz')


  comp_plot <-
    plot(1:(q-P-1),log_p_val, col = "red", type = "b",
         ylim = c(0,1), ylab = "p-value", lwd = 2, pch = 19,xlab = 'Time Points')
  lines(1:(q - P - 1), von_p_val, col = "green", type = "b")
  lines(1:(q - P - 1), gom_p_val, col = "magenta", type = "b")

  # Add a single legend for all lines
  legend("topright", legend = c("Logistic P value Plot", "VB P value Plot", "Gompertz P value Plot"),
         col = c("red", "green", "magenta"), lty = 1, lwd = 2, pch = 19)

  # Add a horizontal line at the specified tolerance level
  abline(h = 0.05, lwd = 2, col = "black", lty = 2)


}




#' P value calculation function for models used in comparison_plot() function
#' @param data_path The path to the data file.
#' @param window_size The size of the window for analysis.
#' @param model The model to be used for p-value calculation.
#' @param parameter The parameters for the model.
#' @param p_val_method The method to be used for p-value calculation ("Parametric" or "Non-Parametric").
#' @param repetition The number of repetitions for the analysis.
#' @return A list containing the calculated p-values and other related information.
#' @export
#' @importFrom stats optim sd pnorm
#' @importFrom graphics par abline legend lines points segments
#' @importFrom utils read.table read.csv

p_value_calc <- function(data_path,
                    window_size = 3,
                    model,
                    parameter,
                    p_val_method,
                    repetition){

  if (!is.list(parameter)) stop("parameter must be a list.")

  ext <- tools::file_ext(data_path)

  # Read data based on file type
  if (ext == "txt") {
    data_r <- as.matrix(read.table(data_path, header = FALSE))
  } else if (ext == "csv") {
    data_r <- as.matrix(read.csv(data_path))
  } else {
    stop("Unsupported file type. Only TXT and CSV files are supported.")
  }

  P<- window_size
  q <- ncol(data_r)
  n <- nrow(data_r)
  t<- 1:q
  x_0 <- mean(data_r[,1])
  x01 <- x_0-x_0*0.1
  x02 <- x_0+x_0*0.1
  K <- mean(data_r[,q])
  K1 <- K-K*0.3
  K2 <- K+K*0.3
  B_ <- repetition

  fun_likelihood = function(param){
    if (model == 'Logistic'){
      r = param[1]
      K = param[2]
      x_0 = param[3]
      sigma2  = param[4]
      rho = param[5]
      mu_lik = K/(1+(K/x_0 - 1)*exp(-r*t))
    }
    else if(model == 'Von-bertalanffy'){
      r = param[1]
      K = param[2]
      x_0 = param[3]
      sigma2  = param[4]
      rho = param[5]
      mu_lik = K*((1 + ((x_0/K)^(1/3) - 1)*exp(-r*(t/3)))^3)
    } else if(model == 'Gompertz'){
      b = param[1]
      c = param[2]
      x_0 = param[3]
      sigma2  = param[4]
      rho = param[5]
      mu_lik = x_0*exp((b/c)*(1-exp(-c*t)))
    }

    cov_mat = matrix(data = NA, nrow = q, ncol = q)
    for (i in 1:q) {
      for (j in 1:q) {
        cov_mat[i,j] = sigma2*rho^(abs(i-j))
      }
    }
    exponent = 0
    for(i in 1:n){
      exponent = exponent + t(data_r[i,]-mu_lik)%*%(solve(cov_mat))%*%(data_r[i,]-mu_lik)
    }
    log_lik = - (q/2)*log(2*pi) - (1/2)*log(det(cov_mat)) - (1/2)*exponent
    return(-log_lik)
  }


  if (model == 'Logistic') {
    param_0 <- c(parameter$log_r, K, x_0, 3, 0.3)
    lower <- c(0.001, K1, x01, 0.0001, -0.001)
    upper <- c(3, K2, x02, 10, 0.999)
  } else if (model == 'Von-bertalanffy') {
    param_0 <- c(parameter$von_r, K, x_0, 3, 0.3)
    lower <- c(0.001, K1, x01, 0.0001, -0.001)
    upper <- c(3, K2, x02, 10, 0.999)
  } else if (model == 'Gompertz') {
    param_0 <- c(parameter$b, parameter$c,x_0,3,0.3)
    lower <-  c(0.01, 0.01, x01, 0.001, -0.001)
    upper <- c(10, 10, x02, 10, 0.999)
  }

  out_global <- optim(param_0, fun_likelihood, method = "L-BFGS-B", hessian = TRUE, lower = lower, upper = upper)
  global_estimate <- out_global$par
  cov_mat_global <- solve(out_global$hessian)


  fun_local_likelihood <- function(loc_parameter) {
    # Extract model-specific parameters
    if (model == 'Logistic') {
      r <- loc_parameter[1]
      K <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu <- K / (1 + (K / x_0 - 1) * exp(-r * t))
    }  else if (model == 'Von-bertalanffy') {
      r <- loc_parameter[1]
      K <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
    } else if (model == 'Gompertz') {
      b <- loc_parameter[1]
      c <- loc_parameter[2]
      sigma2 <- loc_parameter[3]
      rho <- loc_parameter[4]
      mu <- x_0 * exp((b / c) * (1 - exp(-c * t)))
    } else {
      stop("Invalid model specified.")
    }

    # Ensure indices are within bounds
    if (ind + window_size - 1 > length(t)) {
      stop("Window size is too large for the remaining data points.")
    }
    if (ind + window_size - 1 > ncol(data_r)) {
      stop("Window size is too large for the data dimensions.")
    }

    # Create covariance matrix
    cov_mat_ = matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:q) {
      for (j in 1:q) {
        cov_mat_[i,j] = sigma2*rho^(abs(i-j))
      }
    }

    # localized estimates
    mu_local = mu[ind:(ind + window_size-1)]
    data_local = data_r[,ind:(ind + window_size - 1)]
    cov_mat_local = cov_mat_[ind:(ind+window_size-1), ind:(ind+window_size-1)] # Ensure proper dimension for covariance matrix

    # Check if covariance matrix dimensions match
    if (nrow(cov_mat_local) != window_size || ncol(cov_mat_local) != window_size) {
      stop("Covariance matrix dimensions do not match the expected size.")
    }

    exponent <- 0
    for (i in 1:n) {
      exponent <- exponent + t(data_local[i, ] - mu_local) %*% solve(cov_mat_local) %*% (data_local[i, ] - mu_local)
    }

    log_lik <- - (n * length(mu_local) / 2) * log(2 * pi) - (n / 2) * log(det(cov_mat_local)) - (1 / 2) * exponent
    return(-log_lik)
  }

  # Determine the number of columns and column names for est_local based on the model
  if (model == 'Logistic' || model == 'Von-bertalanffy' ) {
    est_local <- matrix(data = NA, nrow = q - window_size + 1, ncol = 4)
    cov_local<-list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("r", "K", "sigma2", "rho")
  }else if(model == 'Gompertz'){
    est_local <- matrix(data = NA, nrow = q - window_size + 1, ncol = 4)
    cov_local<-list() # to store the local variance-covariance matrix.}
    colnames(est_local) <- c("b", "c", "sigma2", "rho")
  }

  # ISRP estimates by maximizing the localized likelihood function
  for (j in 1:(q - window_size + 1)) {
    ind <- j
    if (model == 'Logistic') {
      param_0 <- c(parameter$log_r,K,3,0.3)
      lower <- c(0.001, K1, 0.0001, -0.001)
      upper <- c(3, K2, 3, 0.999)
    } else if (model == 'Von-bertalanffy') {
      param_0 <- c(parameter$von_r,K,3,0.3)
      lower <- c(0.001, K1, 0.0001, -0.001)
      upper <- c(3, K2, 3, 0.999)
    } else if (model == 'Gompertz') {
      param_0 <- c(parameter$b, parameter$c, 2, 0.3)
      lower <- c(0.01, 0.01, 0.001, -0.001)
      upper <- c(5, 10, 3, 0.999)
    }

    opt_results <- optim(
      par = param_0,
      fn = fun_local_likelihood,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      hessian = TRUE
    )

    cov_est_local <- solve(opt_results$hessian)
    cov_local[[j]] <- cov_est_local
    est_local[ind, ] <- opt_results$par

  }

  if (model == 'Logistic'){
    # Computation of Model Specific RGR for each time points (1,2,3,...,q-P+1)
    RGR_hat_log = numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      RGR_hat_log[i] = est_local[i,1]*(1- mean(data_r[,i])/est_local[i,2])
    }

    # Computation of V
    V_log = numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      V_log[i] = log(1/RGR_hat_log[i] - 1/est_local[i,1])
    }


    # Computation of the test statistic
    test_stat = numeric(length = q-P-1)
    for(i in 1:(q-P-1)){
      test_stat[i] = (V_log[i+2]-2*V_log[i+1] + V_log[i])
    }

  }else if (model == 'Von-bertalanffy'){
    # Computation of Model Specific RGR for each time points (1,2,3,...,q-P+1)
    RGR_hat_VB = numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      RGR_hat_VB[i] = est_local[i,1]*((est_local[i,2]/mean(data_r[,i]))^(1/3))*(1- (mean(data_r[,i])/est_local[i,2])^(1/3))
    }

    # Computation of V
    V_VB = numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      V_VB[i] = log(1/RGR_hat_VB[i] + 1/est_local[i,1])
    }

    # Computation of the test statistic
    test_stat = numeric(length = q-P-1)
    for(i in 1:(q-P-1)){
      test_stat[i] = (V_VB[i+2]-2*V_VB[i+1] + V_VB[i])
    }
  }else if (model == 'Gompertz'){
    # Computation of Model Specific RGR for each time points (1,2,3,...,q-P+1)
    RGR_hat_gom = numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      RGR_hat_gom[i] = est_local[i,2]*exp(-est_local[i,1]*t[i])
    }

    # Computation of V
    V_gom = numeric(length = q-P+1)
    for(i in 1:(q-P+1)){
      V_gom[i] = log(1/RGR_hat_gom[i])
    }

    # Computation of the test statistic
    test_stat = numeric(length = q-P-1)

    for(i in 1:(q-P-1)){
      test_stat[i] = (V_gom[i+2]-2*V_gom[i+1] + V_gom[i])
    }
  }
  #######################

  # standard deviation computation of test statistics by bootstrapping

  ############ Parametric bootstrapping
  # Parametric bootstrapping to derive the standard deviation (sd) of the test statistics

  if (model == 'Logistic') {
    r = out_global$par[1]
    K = out_global$par[2]
    x_0 = out_global$par[3]
    sigma2 = out_global$par[4]
    rho = out_global$par[5]
    mu_r <- K / (1 + (K / x_0 - 1) * exp(-r * t))
  } else if (model == 'Von-bertalanffy') {
    r = out_global$par[1]
    K = out_global$par[2]
    x_0 = out_global$par[3]
    sigma2 = out_global$par[4]
    rho = out_global$par[5]
    mu_r <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
  } else if (model == 'Gompertz') {
    b = out_global$par[1]
    c = out_global$par[2]
    x_0 = out_global$par[3]
    sigma2 = out_global$par[4]
    rho = out_global$par[5]
    mu_r <- x_0 * exp((b / c) * (1 - exp(-c * t)))
  }

  cov_mat = matrix(data = NA, nrow = q, ncol = q)
  for (i in 1:q) {
    for (j in 1:q) {
      cov_mat[i,j] = sigma2*rho^(abs(i-j))
    }
  }


  #parametric
  if (p_val_method == 'Parametric'){
    ############ Parametric bootstrapping
    B <- repetition # Parametric bootstrapping
    est_local_log_list = list()
    x0_mean_vals_log = numeric(length = B)
    data_boot_list = list()
    for(l in 1:B){
      #
      data_boot = mvtnorm::rmvnorm(n=1, mean = mu_r, sigma = cov_mat)
      data_boot_list[[l]] = data_boot

      est_local_log_boot = matrix(data = NA, nrow = q-P+1, ncol = length(param_0))
      x0_mean_vals_log = mean(data_boot[,1])

      fun_local_likelihood <- function(loc_parameter) {
        # Extract model-specific parameters
        if (model == 'Logistic') {
          r <- loc_parameter[1]
          K <- loc_parameter[2]
          sigma2 <- loc_parameter[3]
          rho <- loc_parameter[4]
          mu <- K / (1 + (K / x_0 - 1) * exp(-r * t))
        }  else if (model == 'Von-bertalanffy') {
          r <- loc_parameter[1]
          K <- loc_parameter[2]
          sigma2 <- loc_parameter[3]
          rho <- loc_parameter[4]
          mu <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
        } else if (model == 'Gompertz') {
          b <- loc_parameter[1]
          c <- loc_parameter[2]
          sigma2 <- loc_parameter[3]
          rho <- loc_parameter[4]
          mu <- x_0 * exp((b / c) * (1 - exp(-c * t)))
        } else {
          stop("Invalid model specified.")
        }

        # Ensure indices are within bounds
        if (ind + P - 1 > length(t)) {
          stop("Window size is too large for the remaining data points.")
        }
        if (ind + P - 1 > ncol(data_r)) {
          stop("Window size is too large for the data dimensions.")
        }

        # Create covariance matrix
        cov_mat_ = matrix(data = NA, nrow = length(t), ncol = length(t))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat_[i,j] = sigma2*rho^(abs(i-j))
          }
        }

        # localized estimates
        mu_local_boot = mu[ind:(ind + P-1)]
        data_local = data_boot[,ind:(ind + P - 1)]
        cov_mat_local = cov_mat_[ind:(ind+P-1), ind:(ind+P-1)] # Ensure proper dimension for covariance matrix

        # Check if covariance matrix dimensions match
        if (nrow(cov_mat_local) != P || ncol(cov_mat_local) != P) {
          stop("Covariance matrix dimensions do not match the expected size.")
        }

        exponent = 0

        exponent = exponent + t(data_local - mu_local_boot)%*%(solve(cov_mat_local))%*%(data_local - mu_local_boot)


        log_lik <- - (n * length(mu_local_boot) / 2) * log(2 * pi) - (n / 2) * log(det(cov_mat_local)) - (1 / 2) * exponent
        return(-log_lik)
      }

      # Determine the number of columns and column names for est_local based on the model
      if (model == 'Logistic' || model == 'Von-bertalanffy') {
        colnames(est_local_log_boot) <- c("r", "K", "sigma2", "rho")
      }else if(model == 'Gompertz'){
        colnames(est_local_log_boot) <- c("b", "c", "sigma2", "rho")
      }

      # ISRP estimates by maximizing the localized likelihood function
      for (j in 1:(q - P + 1)) {
        ind <- j
        if (model == 'Logistic') {
          param_0 <- c(parameter$log_r, K, 3, 0.3)
          lower = c(0.001, K1, 0.5, -0.001)
          upper = c(3, K2, 10, 0.9)
        } else if (model == 'Von-bertalanffy') {
          param_0 <- c(parameter$von_r, K, 3, 0.3)
          lower = c(0.001, K1, 0.5, -0.001)
          upper = c(3, K2, 10, 0.9)
        } else if (model == 'Gompertz') {
          param_0 <- c(parameter$b, parameter$c, 3, 0.3)
          lower = c(0.01, 0.01, 0.5, -0.001)
          upper = c(10, 10, 10, 0.9)
        }
        local_out_log_boot <- optim(
          par = param_0,
          fn = fun_local_likelihood,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          hessian = TRUE
        )
        est_local_log_boot[j, ] =  local_out_log_boot$par

        }
      est_local_log_list[[l]] = est_local_log_boot
    }


    if (model == "Logistic") {
      # Estimates
      mat_est_r_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      mat_est_K_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (j in 1:B) {
        mat_est_r_MLE[,j] = est_local_log_list[[j]][,1]
        mat_est_K_MLE[,j] = est_local_log_list[[j]][,2]
      }

      # Estimation of modified RGR and its distribution
      modified_RGR = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        numerator = mat_est_r_MLE[i,]*exp(-mat_est_r_MLE[i,]*t[i])*(mat_est_K_MLE[i,] - x0_mean_vals_log)
        denominator = x0_mean_vals_log + exp(-mat_est_r_MLE[i,]*t[i])*(mat_est_K_MLE[i,] - x0_mean_vals_log)
        modified_RGR[i,] = numerator/denominator
      }
      # Computation of V
      V_boot = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        V_boot[i, ] = log(1/modified_RGR[i,] - 1/ mat_est_r_MLE[i,])
      }

    } else if (model == "Von-bertalanffy") {
      # Estimates
      mat_est_r_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      mat_est_K_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (j in 1:B) {
        mat_est_r_MLE[,j] = est_local_log_list[[j]][,1]
        mat_est_K_MLE[,j] = est_local_log_list[[j]][,2]
      }
      # Estimation of modified RGR and its distribution
      modified_RGR = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        numerator = mat_est_r_MLE[i,]*exp(-mat_est_r_MLE[i,]*(t[i]/3))*
          ((mat_est_K_MLE[i,])^(1/3) - (x0_mean_vals_log)^(1/3))
        denominator = (mat_est_K_MLE[i,])^(1/3) + exp(-mat_est_r_MLE[i,]*(t[i]/3))*
          ((x0_mean_vals_log)^(1/3) - (mat_est_K_MLE[i,])^(1/3))
        modified_RGR[i,] = numerator/denominator
      }
      # Computation of V
      V_boot = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        V_boot[i, ] = log(1/modified_RGR[i,] + 1/ mat_est_r_MLE[i,])
      }
    } else if (model == "Gompertz") {
      # Estimates
      mat_est_b_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      mat_est_c_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (j in 1:B) {
        mat_est_b_MLE[,j] = est_local_log_list[[j]][,1]
        mat_est_c_MLE[,j] = est_local_log_list[[j]][,2]
      }
      # Estimation of modified RGR and its distribution
      modified_RGR = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        modified_RGR[i,] = mat_est_b_MLE[i,]*exp(-mat_est_c_MLE[i,]*t[i])
      }
      # Computation of V
      V_boot = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        V_boot[i, ] = log(1/modified_RGR[i,])
      }


    } else {
      stop("Invalid model type. Please specify 'Logistic', 'Von-Bertalanffy', or 'Gompertz'.")
    }

    # Computation of test stat
    test_stat_boot = matrix(data = NA, nrow = q-P-1, ncol = B)
    for (i in 1:(q-P-1)) {
      test_stat_boot[i,] = V_boot[i+2,] - 2*V_boot[i+1,] + V_boot[i,]
    }

    sd_test_stat = numeric(q-P-1)
    for (i in 1:(q-P-1)) {
      sd_test_stat[i] = sd(test_stat_boot[i,])
    }
    p_value_param = numeric(q-P-1)
    for (i in 1:(q-P- 1)) {
      p_value_param[i] = 1- pnorm(abs(test_stat[i]/sd_test_stat[i]))
    }

    #p_val_plot <- # Plotting of the p-values
     # plot(1:(q-P-1), p_value_param, col = "red", type = "b", ylim = c(0,1), ylab = "p-value",
      #     lwd = 2, pch = 19)
    #abline(h = 0.05, lwd = 2, col = "black", lty = 2)

    return(p_value_param)

  }
  # NON PARAMETRIC
  else if (p_val_method == 'Non - Parametric'){

    ############ Non - Parametric bootstrapping
    B <- repetition
    est_local_log_list_np = list()
    x0_mean_vals_log = numeric(length = B)
    data_boot_list = list()
    for(l in 1:B){
      #
      boot = sample(1:nrow(data_r), size = nrow(data_r), replace = TRUE)
      data_np_boot  = data_r[boot, ]
      data_boot_list[[l]] = data_np_boot

      est_local_log_boot = matrix(data = NA, nrow = q-P+1, ncol = length(param_0))
      x0_mean_vals_log = mean(data_np_boot[,1])

      fun_local_likelihood <- function(loc_parameter) {
        # Extract model-specific parameters
        if (model == 'Logistic') {
          r <- loc_parameter[1]
          K <- loc_parameter[2]
          sigma2 <- loc_parameter[3]
          rho <- loc_parameter[4]
          mu <- K / (1 + (K / x_0 - 1) * exp(-r * t))
        }  else if (model == 'Von-bertalanffy') {
          r <- loc_parameter[1]
          K <- loc_parameter[2]
          sigma2 <- loc_parameter[3]
          rho <- loc_parameter[4]
          mu <- K * ((1 + ((x_0 / K)^(1 / 3) - 1) * exp(-r * (t / 3)))^3)
        } else if (model == 'Gompertz') {
          b <- loc_parameter[1]
          c <- loc_parameter[2]
          sigma2 <- loc_parameter[3]
          rho <- loc_parameter[4]
          mu <- x_0 * exp((b / c) * (1 - exp(-c * t)))
        } else {
          stop("Invalid model specified.")
        }

        # Ensure indices are within bounds
        if (ind + P - 1 > length(t)) {
          stop("Window size is too large for the remaining data points.")
        }
        if (ind + P - 1 > ncol(data_r)) {
          stop("Window size is too large for the data dimensions.")
        }

        # Create covariance matrix
        cov_mat_ = matrix(data = NA, nrow = length(t), ncol = length(t))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat_[i,j] = sigma2*rho^(abs(i-j))
          }
        }

        # localized estimates
        mu_local = mu[ind:(ind + P-1)]
        data_local = data_np_boot[,ind:(ind + P - 1)]
        cov_mat_local = cov_mat_[ind:(ind+P-1), ind:(ind+P-1)] # Ensure proper dimension for covariance matrix

        # Check if covariance matrix dimensions match
        if (nrow(cov_mat_local) != P || ncol(cov_mat_local) != P) {
          stop("Covariance matrix dimensions do not match the expected size.")
        }

        exponent <- 0
        for (i in 1:n) {
          exponent <- exponent + t(data_local[i, ] - mu_local) %*% solve(cov_mat_local) %*% (data_local[i, ] - mu_local)
        }

        log_lik <- - (n * length(mu_local) / 2) * log(2 * pi) - (n / 2) * log(det(cov_mat_local)) - (1 / 2) * exponent
        return(-log_lik)
      }

      # Determine the number of columns and column names for est_local based on the model
      if (model == 'Logistic' || model == 'Von-bertalanffy') {
        colnames(est_local_log_boot) <- c("r", "K", "sigma2", "rho")
      }else if(model == 'Gompertz'){
        colnames(est_local_log_boot) <- c("b", "c", "sigma2", "rho")
      }

      # ISRP estimates by maximizing the localized likelihood function
      for (j in 1:(q - P + 1)) {
        ind <- j
        if (model == 'Logistic') {
          param_0 <- c(r, K, 3, 0.3)
          lower = c(0.001, K1, 0.5, -0.001)
          upper = c(3, K2, 10, 0.9)
        } else if (model == 'Von-bertalanffy') {
          param_0 <- c(r, K, 3, 0.3)
          lower = c(0.001, K1, 0.5, -0.001)
          upper = c(3, K2, 10, 0.9)}
        else if(model == 'Gompertz') {
          param_0 <- c(b, c, 3, 0.3)
          lower = c(0.01, 0.01, 0.5, -0.001)
          upper = c(10, 10, 10, 0.9)
        }
        local_out_log_boot <- optim(
          par = param_0,
          fn = fun_local_likelihood,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          hessian = TRUE
        )
        est_local_log_boot[j, ] =  local_out_log_boot$par

      }
      est_local_log_list_np[[l]] = est_local_log_boot
      }


    if (model == "Logistic") {
      # Estimates
      mat_est_r_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      mat_est_K_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (j in 1:B) {
        mat_est_r_MLE[,j] = est_local_log_list_np[[j]][,1]
        mat_est_K_MLE[,j] = est_local_log_list_np[[j]][,2]
      }
      # Estimation of modified RGR and its distribution
      modified_RGR = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        numerator = mat_est_r_MLE[i,]*exp(-mat_est_r_MLE[i,]*t[i])*(mat_est_K_MLE[i,] - x0_mean_vals_log)
        denominator = x0_mean_vals_log + exp(-mat_est_r_MLE[i,]*t[i])*(mat_est_K_MLE[i,] - x0_mean_vals_log)
        modified_RGR[i,] = numerator/denominator
      }
      # Computation of V
      V_boot = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        V_boot[i, ] = log(1/modified_RGR[i,] - 1/ mat_est_r_MLE[i,])
      }
    } else if (model == "Von-bertalanffy") {
      # Estimates
      mat_est_r_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      mat_est_K_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (j in 1:B) {
        mat_est_r_MLE[,j] = est_local_log_list_np[[j]][,1]
        mat_est_K_MLE[,j] = est_local_log_list_np[[j]][,2]
      }
      # Estimation of modified RGR and its distribution
      modified_RGR = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        numerator = mat_est_r_MLE[i,]*exp(-mat_est_r_MLE[i,]*(t[i]/3))*
          ((mat_est_K_MLE[i,])^(1/3) - (x0_mean_vals_log)^(1/3))
        denominator = (mat_est_K_MLE[i,])^(1/3) + exp(-mat_est_r_MLE[i,]*(t[i]/3))*
          ((x0_mean_vals_log)^(1/3) - (mat_est_K_MLE[i,])^(1/3))
        modified_RGR[i,] = numerator/denominator
      }
      # Computation of V
      V_boot = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        V_boot[i, ] = log(1/modified_RGR[i,] + 1/ mat_est_r_MLE[i,])
      }
    } else if (model == "Gompertz") {
      # Estimates
      mat_est_b_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      mat_est_c_MLE = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (j in 1:B) {
        mat_est_b_MLE[,j] = est_local_log_list_np[[j]][,1]
        mat_est_c_MLE[,j] = est_local_log_list_np[[j]][,2]
      }
      # Estimation of modified RGR and its distribution
      modified_RGR = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        modified_RGR[i,] = mat_est_b_MLE[i,]*exp(-mat_est_c_MLE[i,]*t[i])
      }
      # Computation of V
      V_boot = matrix(data = NA, nrow = q-P+1, ncol = B)
      for (i in 1:(q-P+1)) {
        V_boot[i, ] = log(1/modified_RGR[i,])
      }
    } else {
      stop("Invalid model type. Please specify 'Logistic', 'Von-bertalanffy', or 'Gompertz'.")
    }

    # Computation of test stat
    test_stat_boot = matrix(data = NA, nrow = q-P-1, ncol = B)
    for (i in 1:(q-P-1)) {
      test_stat_boot[i,] = V_boot[i+2,] - 2*V_boot[i+1,] + V_boot[i,]
    }
    sd_test_stat = numeric(q-P-1)
    for (i in 1:(q-P-1)) {
      sd_test_stat[i] = sd(test_stat_boot[i,])
    }
    p_value_non_param = numeric(q-P-1)
    for (i in 1:(q-P- 1)) {
      p_value_non_param[i] = 1- pnorm(abs(test_stat[i]/sd_test_stat[i]))
    }

   # p_val_plot <- # Plotting of the p-values
    #  plot(1:(q-P-1), p_value_non_param, col = "red", type = "b", ylim = c(0,1), ylab = "p-value",
     #      lwd = 2, pch = 19)
    #abline(h = 0.05, lwd = 2, col = "black", lty = 2)

    return(p_value_non_param)

  }

}


