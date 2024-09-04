
server <- function(input, output, session) {




  observeEvent(input$button1,{
    show_toast('Important',
               text = 'Dont click this button twice, Computation May take little time',
               timer = 11000,
               position = "center",
               type = 'success')
  })




  # <-------------------------------LOGISTIC (STAART) -------------------------------------->
  # plot function




  observeEvent(input$t_range,{
    if ((length(seq(from = input$t_range[1], to = input$t_range[2], by = 1))) >=80){
      (show_toast(
        title = "CAUTION ⚠️",
        text = "Please Select Time points range less than 80",
        timer = 5000,
        position = 'center',
        type = 'info'
      ))
    }})

  observeEvent(input$button1,{
    output$model_name <- renderText({
      selected_model <- input$model

      if (selected_model == "Logistic Model") {
        return("Estimation using Logistic Model")
      } else if (selected_model == "Exponential Model") {
        return("Estimation using Exponential Model")
      } else if (selected_model == "Von-Bertallanfy Model") {
        return("Estimation using Von-Bertallanfy Model")
      } else if (selected_model == "Theta - Logistic Model") {
        return("Estimation using Theta - Logistic Model")
      } else if (selected_model == "Gompertz Model") {
        return("Estimation using Gompertz Model")
      } else {
        return("Please select a model.")
      }
    })
  })

  traj_plot <- function(K, x0, r, sigma2, rho, n, t_range) {

    if (K == 0 || x0 == 0 || r == 0 || sigma2 == 0 || rho == 0 || n == 0) {
      stop(show_toast(
        title = "Error ❌",
        text = "All / Some input values should not be zero.",
        timer = 10000,
        position = 'center',
        type = 'warning'
      ))
    }

    t <- seq(from = t_range[1], to = t_range[2], by = 1)

    h <- numeric(length = length(t) - 1)
    for (i in 1:(length(t) - 1)) {
      h[i] <- t[i + 1] - t[i]
    }

    mu <- K / (1 + (K / x0 - 1) * exp(-r * t))

    cov_mat <- matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:length(t)) {
      for (j in 1:length(t)) {
        cov_mat[i, j] <- sigma2 * rho^(abs(i - j))
      }
    }



    data <- rmvnorm(n, mean = mu, sigma = cov_mat)

    par(mar = c(4,4,2,1))
    plot <- matplot(t, t(data), type = "l", xlab = "Time", ylab = "Population", main = "Population Growth Simulation",cex.lab = 1.2)
    return(list(plot = plot, data = data))
  }



  traj_plot_expo <- function(x0, r,sigma2,rho, n, t_range){

    if (x0 == 0 || r == 0 || sigma2 == 0 || rho == 0 || n == 0) {
      stop(show_toast(
        title = "Error ❌",
        text = "All / Some input values should not be zero.",
        timer = 10000,
        position = 'center',
        type = 'warning'
      ))
    }

    t <- seq(from = t_range[1], to = t_range[2], by = 1)

    h <- numeric(length = length(t) - 1)
    for (i in 1:(length(t) - 1)) {
      h[i] <- t[i + 1] - t[i]
    }

    mu <-  x0*exp(r*t)

    cov_mat <- matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:length(t)) {
      for (j in 1:length(t)) {
        cov_mat[i, j] <- sigma2 * rho^(abs(i - j))
      }
    }


    data <- rmvnorm(n, mean = mu, sigma = cov_mat)




    par(mar = c(4,4,2,1))
    plot <- matplot(t, t(data), type = "l", xlab = "Time", ylab = "Population", main = "Population Growth Simulation",cex.lab = 1.3)
    return(list(plot = plot, data = data))

  }


  traj_plot_theta <- function (K, x0, theta,r, sigma2, rho, n, t_range){
    if (K == 0 || x0 == 0 || r == 0 || sigma2 == 0 || rho == 0 || theta == 0 || n == 0) {
      stop(show_toast(
        title = "Error ❌",
        text = "All / Some input values should not be zero.",
        timer = 10000,
        position = 'center',
        type = 'warning'
      ))
    }

    t <- seq(from = t_range[1], to = t_range[2], by = 1)

    h = numeric(length = length(t)-1)
    for(i in 1:(length(t)-1)){
      h[i] = t[i+1] - t[i]
    }

    # theta-logistic mean function
    mu_theta = K/(1 + ((K/x0)^theta - 1)*exp(-r*theta*t))


    cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:length(t)) {
      for (j in 1:length(t)) {
        cov_mat[i,j] = sigma2*rho^(abs(i-j))
      }
    }


    data <- rmvnorm(n, mean = mu_theta, sigma = cov_mat)
    par(mar = c(4,4,2,1))
    plot <- matplot(t, t(data), type = "l", xlab = "Time", ylab = "Population", main = "Population Growth Simulation",cex.lab = 1.3)
    return(list(plot = plot, data = data))
  }


  traj_plot_von <- function (K, x0, r, sigma2, rho, n, t_range){


    if (K == 0 || x0 == 0 || r == 0 || sigma2 == 0 || rho == 0 || n == 0) {
      stop(show_toast(
        title = "Error ❌",
        text = "All / Some input values should not be zero.",
        timer = 10000,
        position = 'center',
        type = 'warning'
      ))
    }

    t <- seq(from = t_range[1], to = t_range[2], by = 1)

    h <- numeric(length = length(t) - 1)
    for (i in 1:(length(t) - 1)) {
      h[i] <- t[i + 1] - t[i]
    }

    mu <- K*((1 + ((x0/K)^(1/3) - 1)*exp(-r*(t/3)))^3)

    cov_mat <- matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:length(t)) {
      for (j in 1:length(t)) {
        cov_mat[i, j] <- sigma2 * rho^(abs(i - j))
      }
    }


    data <- rmvnorm(n, mean = mu, sigma = cov_mat)
    par(mar = c(4,4,2,1))
    plot <- matplot(t, t(data), type = "l", xlab = "Time", ylab = "Population", main = "Population Growth Simulation",cex.lab = 1.3)
    return(list(plot = plot, data = data))

  }

  traj_plot_gom <- function (b, x0, c, sigma2, rho, n, t_range){

    if (b == 0 || x0 == 0 || c == 0 || sigma2 == 0 || rho == 0 || n == 0) {
      stop(show_toast(
        title = "Error ❌",
        text = "All / Some input values should not be zero.",
        timer = 10000,
        position = 'center',
        type = 'warning'
      ))
    }

    t <- seq(from = t_range[1], to = t_range[2], by = 1)

    h <- numeric(length = length(t) - 1)
    for (i in 1:(length(t) - 1)) {
      h[i] <- t[i + 1] - t[i]
    }

    mu <- x0*exp((b/c)*(1-exp(-c*t)))

    cov_mat <- matrix(data = NA, nrow = length(t), ncol = length(t))
    for (i in 1:length(t)) {
      for (j in 1:length(t)) {
        cov_mat[i, j] <- sigma2 * rho^(abs(i - j))
      }
    }



    data <- rmvnorm(n, mean = mu, sigma = cov_mat)
    par(mar = c(4,4,2,1))
    plot <- matplot(t, t(data), type = "l", xlab = "Time", ylab = "Population", main = "Population Growth Simulation")
    return(list(plot = plot, data = data))


  }

  # traj_plot using logistic


  generatePlot <- eventReactive(input$button1, {
    if (input$model == 'Logistic Model'){

      traj_plot(input$logi_k, input$logi_in, input$logi_r, input$logi_s, input$logi_rho, input$trajectory, input$t_range)$plot}
    else if (input$model == 'Exponential Model'){
      traj_plot_expo(input$expo_in, input$expo_r, input$expo_s, input$expo_rho, input$trajectory, input$t_range)$plot
    } else if (input$model == 'Theta - Logistic Model'){
      traj_plot_theta(input$theta_k, input$theta_in,  input$theta_th,input$theta_r, input$theta_s, input$theta_rho,input$trajectory, input$t_range)$plot
    } else if (input$model == 'Von-Bertallanfy Model'){
      traj_plot_von(input$von_k, input$von_in, input$von_r, input$von_s, input$von_rho, input$trajectory, input$t_range)$plot
    } else if (input$model == 'Gompertz Model'){
      traj_plot_gom(input$gom_b, input$gom_in, input$gom_c, input$gom_s, input$gom_rho, input$trajectory, input$t_range)$plot
    }
  })

  output$plot1 <- renderPlot({
    if (input$model == 'Logistic Model'){
      generatePlot()}
    else if (input$model == 'Exponential Model'){
      generatePlot()
    }
    else if (input$model == 'Theta - Logistic Model'){
      generatePlot()
    } else if (input$model == 'Von-Bertallanfy Model'){
      generatePlot()
    }
    else if (input$model == 'Gompertz Model'){
      generatePlot()
    }
  })
  # zooming in the plot



  observeEvent(input$zoom, {
    if (input$model == 'Logistic Model'){

      showModal(modalDialog(

        renderPlot({
          traj_plot(input$logi_k, input$logi_in, input$logi_r, input$logi_s, input$logi_rho, input$trajectory, input$t_range)$plot
        }, height = 600),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))}
    else if (input$model == 'Exponential Model'){
      showModal(modalDialog(

        renderPlot({
          traj_plot_expo(input$expo_in, input$expo_r, input$expo_s, input$expo_rho, input$trajectory, input$t_range)$plot
        }, height = 600),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    }
    else if (input$model == 'Theta - Logistic Model'){
      showModal(modalDialog(

        renderPlot({
          traj_plot_theta(input$theta_k, input$theta_in, input$theta_th,input$theta_r, input$theta_s, input$theta_rho, input$trajectory, input$t_range)$plot
        }, height = 600),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    }
    else if (input$model == 'Von-Bertallanfy Model'){
      showModal(modalDialog(

        renderPlot({
          traj_plot_von(input$von_k, input$von_in, input$von_r, input$von_s, input$von_rho, input$trajectory, input$t_range)$plot
        }, height = 600),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    }
    else if (input$model == 'Gompertz Model'){
      showModal(modalDialog(

        renderPlot({

          traj_plot_gom(input$gom_b, input$gom_in, input$gom_c, input$gom_s, input$gom_rho, input$trajectory, input$t_range)$plot

        }, height = 600),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    }

  })

  # download handler for plot
  output$plot_down <- downloadHandler(
    filename = function() {
      if (input$model == 'Logistic Model') {
        if (input$format == "png") {
          paste("data_plot_logi-", Sys.Date(), ".png", sep = "")
        } else if (input$format == "pdf") {
          paste("data_plot_logi-", Sys.Date(), ".pdf", sep = "")
        } else if (input$format == "csv") {
          paste("data_logi-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Exponential Model') {
        if (input$format == "png") {
          paste("data_plot_expo-", Sys.Date(), ".png", sep = "")
        } else if (input$format == "pdf") {
          paste("data_plot_expo-", Sys.Date(), ".pdf", sep = "")
        } else if (input$format == "csv") {
          paste("data_expo-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format == "png") {
          paste("data_plot_theta-", Sys.Date(), ".png", sep = "")
        } else if (input$format == "pdf") {
          paste("data_plot_theta-", Sys.Date(), ".pdf", sep = "")
        } else if (input$format == "csv") {
          paste("data_theta-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format == "png") {
          paste("data_plot_von-", Sys.Date(), ".png", sep = "")
        } else if (input$format == "pdf") {
          paste("data_plot_von-", Sys.Date(), ".pdf", sep = "")
        } else if (input$format == "csv") {
          paste("data_von-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Gompertz Model'){
        if (input$format == "png") {
          paste("data_plot_gom-", Sys.Date(), ".png", sep = "")
        } else if (input$format == "pdf") {
          paste("data_plot_gom-", Sys.Date(), ".pdf", sep = "")
        } else if (input$format == "csv") {
          paste("data_gom-", Sys.Date(), ".csv", sep = "")
        }
      }
    },
    content = function(file) {
      if (input$model == 'Logistic Model') {
        if (input$format == "png") {
          png(file)
          traj_plot(input$logi_k, input$logi_in, input$logi_r, input$logi_s, input$logi_rho, input$trajectory, input$t_range)$plot
          dev.off()
        } else if (input$format == "pdf") {
          pdf(file)
          traj_plot(input$logi_k, input$logi_in, input$logi_r, input$logi_s, input$logi_rho, input$trajectory, input$t_range)$plot
          dev.off()
        } else if (input$format == "csv") {
          write.csv(traj_plot(input$logi_k, input$logi_in, input$logi_r, input$logi_s, input$logi_rho, input$trajectory, input$t_range)$data,file)

        }
      } else if (input$model == 'Exponential Model') {
        if (input$format == "png") {
          png(file)
          traj_plot_expo(input$expo_in, input$expo_r, input$expo_s, input$expo_rho, input$trajectory, input$t_range)$plot
          dev.off()
        } else if (input$format == "pdf") {
          pdf(file)
          traj_plot_expo(input$expo_in, input$expo_r, input$expo_s, input$expo_rho, input$trajectory, input$t_range)$plot
          dev.off()
        }  else if (input$format == "csv") {

          write.csv(traj_plot_expo(input$expo_in, input$expo_r, input$expo_s, input$expo_rho, input$trajectory, input$t_range)$data,file)

        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format == "png") {
          png(file)
          traj_plot_theta(input$theta_k, input$theta_in, input$theta_th,input$theta_r, input$theta_s, input$theta_rho, input$trajectory, input$t_range)$plot
          dev.off()
        } else if (input$format == "pdf") {
          pdf(file)
          traj_plot_theta(input$theta_k, input$theta_in, input$theta_th,input$theta_r, input$theta_s, input$theta_rho, input$trajectory, input$t_range)$plot
          dev.off()
        } else if (input$format == "csv") {

          write.csv(traj_plot_theta(input$theta_k, input$theta_in, input$theta_th,input$theta_r, input$theta_s, input$theta_rho, input$trajectory, input$t_range)$data,file)

        }
      } else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format == "png") {
          png(file)
          traj_plot_von(input$von_k, input$von_in, input$von_r, input$von_s, input$von_rho, input$trajectory, input$t_range)$plot
          dev.off()
        } else if (input$format == "pdf") {
          pdf(file)
          traj_plot_von(input$von_k, input$von_in, input$von_r, input$von_s, input$von_rho, input$trajectory, input$t_range)$plot
          dev.off()
        } else if (input$format == "csv") {

          write.csv(traj_plot_von(input$von_k, input$von_in, input$von_r, input$von_s, input$von_rho, input$trajectory, input$t_range)$data,file)

        }
      } else if (input$model == 'Gompertz Model'){
        if (input$format == "png") {
          png(file)
          traj_plot_gom(input$gom_b, input$gom_in, input$gom_c, input$gom_s, input$gom_rho, input$trajectory, input$t_range)$plot

          dev.off()
        } else if (input$format == "pdf") {
          pdf(file)
          traj_plot_gom(input$gom_b, input$gom_in, input$gom_c, input$gom_s, input$gom_rho, input$trajectory, input$t_range)$plot

          dev.off()
        } else if (input$format == "csv") {

          write.csv(traj_plot_gom(input$gom_b, input$gom_in, input$gom_c, input$gom_s, input$gom_rho, input$trajectory, input$t_range)$data,file)

        }
      }
    }
  )





  # Likelihood function of Logistic
  # GLOBAL ESTIMATES

  data_global <- reactiveVal(NULL)
  data_global_expo <- reactiveVal(NULL)
  data_global_theta <- reactiveVal(NULL)
  data_global_von <- reactiveVal(NULL)
  data_global_gom <- reactiveVal(NULL)

  observeEvent(input$button1, {
    if (input$model == 'Logistic Model'){

      param <- c(0.3, input$logi_k, 2,0.5)
      data0 <- traj_plot(input$logi_k, input$logi_in, input$logi_r, input$logi_s, input$logi_rho, input$trajectory, input$t_range)$data
      t0 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n0 <- input$trajectory
      x00 <- input$logi_in

      # Global Likelihood

      fun_likelihood = function(param){
        r = param[1]
        K = param[2]
        sigma2  = param[3]
        rho = param[4]
        q = length(t0)

        mu = K/(1 + (K/x00 - 1)*exp(-r*t0))

        cov_mat = matrix(data = NA, nrow = length(t0), ncol = length(t0))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }



        exponent = 0
        for(i in 1:n0){
          exponent = exponent + t((data0[i,])-(mu))%*%(solve(cov_mat))%*%(data0[i,]-mu)
        }


        log_lik = - (n0*q/2)*log(2*pi) - (n0/2)*log(det(cov_mat)) - (1/2)*exponent
        return(-log_lik)
      }

      param0 <- (param)

      out0 <- optim(param0, fun_likelihood, method = "L-BFGS-B",
                    hessian = TRUE, lower = c(0.001, input$logi_k-(input$logi_k*0.2), 0.0001, -0.001),
                    upper = c(3, input$logi_k+(input$logi_k*0.2), 10, 0.999))
      cov <- data.frame(solve(out0$hessian))
      colnames(cov) <- c('r','K','Sigma','Rho')

      # Create a data frame with parameter names
      result_df <- data.frame(r = out0$par[1], K = out0$par[2], sigma = out0$par[3], rho = out0$par[4])


      result_list <- list(estimates = result_df, cov_matrix = cov)
      data_global(result_list)


      output$table1 <- renderTable({
        result_df_formatted <- result_df  # Create a copy of the data frame to format
        # Format numeric columns to display 5 decimal places
        numeric_columns <- sapply(result_df_formatted, is.numeric)
        result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
        result_df_formatted
      }, width = "80%", align = 'c',height = 1000)


      output$covM <- renderTable({
        # Create a data frame from the covariance matrix
        cov_df <- as.data.frame(cov)

        # Rename the columns
        colnames(cov_df) <- c('r', 'K', 'Sigma', 'rho')

        # Format the numeric columns to display 5 decimal places
        cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

        cov_df
      }, width = "80%", align = 'c',height= 1000)}

    else if (input$model == 'Exponential Model'){

      param_exp <- c(input$expo_r, 2, 0.5)
      data0 <- traj_plot_expo(input$expo_in, input$expo_r, input$expo_s, input$expo_rho, input$trajectory, input$t_range)$data
      t0 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n0 <- input$trajectory
      x00 <- input$expo_in
      q <- length(t0)

      fun_likelihood = function(param){
        r = param[1]
        sigma2  = param[2]
        rho = param[3]
        mu = x00*exp(r*t0) # Exponential mean function

        cov_mat = matrix(data = NA, nrow = length(t0), ncol = length(t0))
        for (i in 1:length(t0)) {
          for (j in 1:length(t0)) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }


        exponent = 0
        for(i in 1:n0){
          exponent = exponent + t(data0[i,]-mu)%*%(solve(cov_mat))%*%(data0[i,]-mu)
        }





        log_lik = - (n0*length(t0)/2)*log(2*pi) - (n0/2)*log(det(cov_mat)) - (1/2)*exponent
        return(-log_lik)
      }
      param_0 <- param_exp
      out0 <- optim(param_0, fun_likelihood, method = "L-BFGS-B",
                    hessian = TRUE, lower = c(0.001, 0.0001, -0.001),
                    upper = c(3, 10, 0.999))
      cov <- data.frame(solve(out0$hessian))
      colnames(cov) <- c('r','Sigma','Rho')

      # Create a data frame with parameter names
      result_df <- data.frame(r = out0$par[1], sigma = out0$par[2], rho = out0$par[3])


      result_list <- list(estimates_expo = result_df, cov_matrix_expo = cov)
      data_global_expo(result_list)


      output$table1 <- renderTable({
        result_df_formatted <- result_df  # Create a copy of the data frame to format
        # Format numeric columns to display 5 decimal places
        numeric_columns <- sapply(result_df_formatted, is.numeric)
        result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
        result_df_formatted
      }, width = "80%", align = 'c',height = 1000)


      output$covM <- renderTable({
        # Create a data frame from the covariance matrix
        cov_df <- as.data.frame(cov)

        # Rename the columns
        colnames(cov_df) <- c('r', 'Sigma', 'rho')

        # Format the numeric columns to display 5 decimal places
        cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

        cov_df
      }, width = "100%", align = 'c',height= 1000)}

    else if (input$model == 'Theta - Logistic Model'){


      param <- c(0.3, input$theta_k, input$theta_th,2, 0.5)
      data0 <- traj_plot_theta(input$theta_k, input$theta_in, input$theta_th,input$theta_r, input$theta_s, input$theta_rho,input$trajectory, input$t_range)$data
      t0 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n0 <- input$trajectory
      x00 <- input$theta_in
      q  <- length(t0)

      # Global Likelihood

      fun_likelihood = function(param){
        r = param[1]
        K = param[2]
        theta=param[3]
        sigma2  = param[4]
        rho = param[5]
        mu_theta = K/(1 + ((K/x00)^theta - 1)*exp(-r*theta*t0))

        cov_mat = matrix(data = NA, nrow = length(t0), ncol = length(t0))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }
        exponent = 0
        for(i in 1:n0){
          exponent = exponent + t(data0[i,]-mu_theta)%*%(solve(cov_mat))%*%(data0[i,]-mu_theta)
        }


        log_lik = - (n0*q/2)*log(2*pi) - (n0/2)*log(det(cov_mat)) - (1/2)*exponent
        return(-log_lik)
      }

      param_0 <- param
      out0 <- optim(param_0, fun_likelihood, method = "L-BFGS-B",
                    hessian = TRUE, lower = c(0.001, input$theta_k-(input$theta_k*0.2), 0.1, 0.0001, -0.001),
                    upper = c(3, input$theta_k+(input$theta_k*0.2), 2, 10, 0.999))

      cov <- data.frame(solve(out0$hessian))
      colnames(cov) <- c('r','K','Theta','Sigma','Rho')

      # Create a data frame with parameter names
      result_df <- data.frame(r = out0$par[1], K = out0$par[2], Theta = out0$par[3],sigma = out0$par[4], rho = out0$par[5])


      result_list <- list(estimates_theta = result_df, cov_matrix_theta = cov)
      data_global_theta(result_list)

      output$table1 <- renderTable({
        result_df_formatted <- data_global_theta()$estimates_theta  # Create a copy of the data frame to format
        # Format numeric columns to display 5 decimal places
        numeric_columns <- sapply(result_df_formatted, is.numeric)
        result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))

        result_df_formatted
      }, width = "80%", align = 'c',height = 1000)


      output$covM <- renderTable({
        # Create a data frame from the covariance matrix
        cov_df <- as.data.frame(data_global_theta()$cov_matrix_theta)

        # Rename the columns
        colnames(cov_df) <- c('r', 'K','Theta', 'Sigma', 'rho')

        # Format the numeric columns to display 5 decimal places
        cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

        cov_df[,1:4]
      }, width = "80%", align = 'c')

    }
    else if (input$model == 'Von-Bertallanfy Model'){
      param <- c(0.8, input$von_k, 2, 0.5)
      data0 <- traj_plot_von(input$von_k, input$von_in, input$von_r, input$von_s, input$von_rho, input$trajectory, input$t_range)$data
      t0 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n0 <- input$trajectory
      x00 <- input$von_in
      q <- length(t0)

      # Global Likelihood

      fun_likelihood = function(param){
        r = param[1]
        K = param[2]
        sigma2  = param[3]
        rho = param[4]
        mu = K*((1 + ((x00/K)^(1/3) - 1)*exp(-r*(t0/3)))^3)

        cov_mat = matrix(data = NA, nrow = length(t0), ncol = length(t0))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }
        exponent = 0
        for(i in 1:n0){
          exponent = exponent + t(data0[i,]-mu)%*%(solve(cov_mat))%*%(data0[i,]-mu)
        }


        log_lik = - (n0*q/2)*log(2*pi) - (n0/2)*log(det(cov_mat)) - (1/2)*exponent
        return(-log_lik)
      }

      param0 <- c(param)

      out0 <- optim(param0, fun_likelihood, method = "L-BFGS-B",
                    hessian = TRUE, lower = c(0.001, input$von_k-(input$von_k*0.2), 0.0001, -0.001),
                    upper = c(3, input$von_k+(input$von_k*0.2), 10, 0.999))
      cov <- data.frame(solve(out0$hessian))
      colnames(cov) <- c('r','K','Sigma','Rho')

      # Create a data frame with parameter names
      result_df <- data.frame(r = out0$par[1], K = out0$par[2], sigma = out0$par[3], rho = out0$par[4])


      result_list <- list(estimates_von = result_df, cov_matrix_von = cov)
      data_global_von(result_list)


      output$table1 <- renderTable({
        result_df_formatted <- data_global_von()$estimates_von  # Create a copy of the data frame to format
        # Format numeric columns to display 5 decimal places
        numeric_columns <- sapply(result_df_formatted, is.numeric)
        result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
        result_df_formatted
      }, width = "80%", align = 'c',height = 1000)


      output$covM <- renderTable({
        # Create a data frame from the covariance matrix
        cov_df <- as.data.frame(data_global_von()$cov_matrix_von)

        # Rename the columns
        colnames(cov_df) <- c('r', 'K', 'Sigma', 'rho')

        # Format the numeric columns to display 5 decimal places
        cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

        cov_df
      }, width = "80%", align = 'c',height= 1000)
    }
    else if (input$model == 'Gompertz Model'){
      param <- c(input$gom_c, input$gom_b, 2,0.5)
      data0 <- traj_plot_gom(input$gom_b, input$gom_in, input$gom_c, input$gom_s, input$gom_rho, input$trajectory, input$t_range)$data
      t0 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n0 <- input$trajectory
      x00 <- input$gom_in
      q <- length(t0)

      # Global Likelihood

      fun_likelihood = function(param){
        c = param[1]
        b = param[2]
        sigma2  = param[3]
        rho = param[4]
        mu = x00*exp((b/c)*(1-exp(-c*t0)))

        cov_mat = matrix(data = NA, nrow = length(t0), ncol = length(t0))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }
        exponent = 0
        for(i in 1:n0){
          exponent = exponent + t(data0[i,]-mu)%*%(solve(cov_mat))%*%(data0[i,]-mu)
        }


        log_lik = - (n0*q/2)*log(2*pi) - (n0/2)*log(det(cov_mat)) - (1/2)*exponent
        return(-log_lik)
      }
      param_0 <- param
      out0 <- optim(param_0, fun_likelihood, method = "L-BFGS-B",
                    hessian = TRUE, lower = c(0.05, 0.05, 0.005, -0.001),
                    upper = c(5, 10, 10, 0.999))

      cov <- data.frame(solve(out0$hessian))
      colnames(cov) <- c("c", "b", "Sigma", "Rho")

      # Create a data frame with parameter names
      result_gom <- data.frame(c = out0$par[1], b = out0$par[2], Sigma = out0$par[3], Rho = out0$par[4])


      result_list <- list(estimates_gom = result_gom, cov_matrix_gom = cov)
      data_global_gom(result_list)


      output$table1 <- renderTable({
        result_df_formatted <- data_global_gom()$estimates_gom  # Create a copy of the data frame to format
        # Format numeric columns to display 5 decimal places
        numeric_columns <- sapply(result_df_formatted, is.numeric)
        result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
        result_df_formatted
      }, width = "80%", align = 'c',height = 1000)


      output$covM <- renderTable({
        # Create a data frame from the covariance matrix
        cov_df <- as.data.frame(data_global_gom()$cov_matrix_gom)

        # Rename the columns
        colnames(cov_df) <- c("c", "b", "Sigma", "Rho")

        # Format the numeric columns to display 5 decimal places
        cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

        cov_df
      }, width = "80%", align = 'c',height= 1000)
    }


    show_toast('Progress',
               text = 'Your Global Estimates Are Done...',
               timer = 5000,
               position = "center",
               type = 'success')

  }



  )




  # Zooming globall estimates

  observeEvent(input$zoom2, {
    if (input$model == 'Logistic Model'){
      showModal(modalDialog(
        paste('Global Estimation Table'),
        br(),
        renderTable({
          result_df_formatted <- data_global()$estimates  # Create a copy of the data frame to format
          # Format numeric columns to display 5 decimal places
          numeric_columns <- sapply(result_df_formatted, is.numeric)
          result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
          result_df_formatted
        },height = 1000,width = "100%",align = 'c'),
        br(),
        br(),
        paste('Covariance Matrix'),
        br(),
        renderTable({
          # Create a data frame from the covariance matrix
          cov_df <- as.data.frame(data_global()$cov_matrix)

          # Rename the columns
          colnames(cov_df) <- c('r', 'K', 'Sigma', 'rho')

          # Format the numeric columns to display 5 decimal places
          cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

          cov_df
        }, width = "100%", align = 'c',height = 1000),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    }else if(input$model == 'Exponential Model'){

      showModal(modalDialog(
        paste('Global Estimation Table'),
        br(),
        renderTable({
          result_df_formatted <- data_global_expo()$estimates_expo  # Create a copy of the data frame to format
          # Format numeric columns to display 5 decimal places
          numeric_columns <- sapply(result_df_formatted, is.numeric)
          result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
          result_df_formatted
        },height = 1000,width = "100%",align = 'c'),
        br(),
        br(),
        paste('Covariance Matrix'),
        br(),
        renderTable({
          # Create a data frame from the covariance matrix
          cov_df <- as.data.frame(data_global_expo()$cov_matrix_expo)

          # Rename the columns
          colnames(cov_df) <- c('r', 'Sigma', 'rho')

          # Format the numeric columns to display 5 decimal places
          cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

          cov_df
        }, width = "100%", align = 'c',height = 1000),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))

    }
    else if (input$model == 'Theta - Logistic Model'){


      showModal(modalDialog(
        paste('Global Estimation Table'),
        br(),
        renderTable({
          result_df_formatted <- data_global_theta()$estimates_theta  # Create a copy of the data frame to format
          # Format numeric columns to display 5 decimal places
          numeric_columns <- sapply(result_df_formatted, is.numeric)
          result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
          result_df_formatted
        },height = 1000,width = "100%",align = 'c'),
        br(),
        br(),
        paste('Covariance Matrix'),
        br(),
        renderTable({
          # Create a data frame from the covariance matrix
          cov_df <- as.data.frame(data_global_theta()$cov_matrix_theta)

          # Rename the columns
          colnames(cov_df) <- c('r','K',"Theta", 'Sigma', 'rho')

          # Format the numeric columns to display 5 decimal places
          cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

          cov_df
        }, width = "100%", align = 'c',height = 1000),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))


    } else if (input$model == 'Von-Bertallanfy Model'){

      showModal(modalDialog(
        paste('Global Estimation Table'),
        br(),
        renderTable({
          result_df_formatted <- data_global_von()$estimates_von  # Create a copy of the data frame to format
          # Format numeric columns to display 5 decimal places
          numeric_columns <- sapply(result_df_formatted, is.numeric)
          result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
          result_df_formatted
        },height = 1000,width = "100%",align = 'c'),
        br(),
        br(),
        paste('Covariance Matrix'),
        br(),
        renderTable({
          # Create a data frame from the covariance matrix
          cov_df <- as.data.frame(data_global_von()$cov_matrix_von)

          # Rename the columns
          colnames(cov_df) <- c('r','K', 'Sigma', 'rho')

          # Format the numeric columns to display 5 decimal places
          cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

          cov_df
        }, width = "100%", align = 'c',height = 1000),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))


    } else if (input$model == 'Gompertz Model'){

      showModal(modalDialog(
        paste('Global Estimation Table'),
        br(),
        renderTable({
          result_df_formatted <- data_global_gom()$estimates_gom  # Create a copy of the data frame to format
          # Format numeric columns to display 5 decimal places
          numeric_columns <- sapply(result_df_formatted, is.numeric)
          result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
          result_df_formatted
        },height = 1000,width = "100%",align = 'c'),
        br(),
        br(),
        paste('Covariance Matrix'),
        br(),
        renderTable({
          # Create a data frame from the covariance matrix
          cov_df <- as.data.frame(data_global_gom()$cov_matrix_gom)

          # Rename the columns
          colnames(cov_df) <- c("c", "b", "Sigma", "Rho")

          # Format the numeric columns to display 5 decimal places
          cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

          cov_df
        }, width = "100%", align = 'c',height = 1000),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))


    }
  })

  # download handler for Global estimation


  output$c_mat <- downloadHandler(
    filename = function() {
      if (input$model == 'Logistic Model' && input$format1_1 == 'csv') {
        paste("Cov_Matrix_logi-", Sys.Date(), ".csv", sep = "")
      } else if (input$model == 'Exponential Model' && input$format1_1 == 'csv') {
        paste("Cov_Matrix_expo-", Sys.Date(), ".csv", sep = "")
      } else if (input$model == 'Theta - Logistic Model' && input$format1_1 == 'csv'){
        paste("Cov_Matrix_theta-", Sys.Date(), ".csv", sep = "")
      } else if (input$model == 'Von-Bertallanfy Model' && input$format1_1 == 'csv'){
        paste("Cov_Matrix_von-", Sys.Date(), ".csv", sep = "")
      }  else if (input$model == 'Gompertz Model' && input$format1_1 == 'csv'){
        paste("Cov_Matrix_gom-", Sys.Date(), ".csv", sep = "")
      }
    },
    content = function(file) {
      if (input$model == 'Logistic Model' && input$format1_1 == 'csv') {
        write.csv(data_global()$cov_matrix, file)
      } else if (input$model == 'Exponential Model' && input$format1_1 == 'csv') {
        write.csv(data_global_expo()$cov_matrix_expo, file)
      } else if (input$model == 'Theta - Logistic Model' && input$format1_1 == 'csv'){
        write.csv(data_global_theta()$cov_matrix_theta, file)
      } else if (input$model == 'Von-Bertallanfy Model' && input$format1_1 == 'csv'){
        write.csv(data_global_von()$cov_matrix_von, file)
      }  else if (input$model == 'Gompertz Model' && input$format1_1 == 'csv'){
        write.csv(data_global_gom()$cov_matrix_gom, file)
      }
    }
  )



  output$down_global <- downloadHandler(
    filename = function() {
      if (input$model == 'Logistic Model') {
        if (input$format1_2 == "csv") {
          paste("Global_Estimates_logi-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Exponential Model') {
        if (input$format1_2 == "csv") {
          paste("Global_Estimates_expo-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format1_2 == "csv") {
          paste("Global_Estimates_theta-", Sys.Date(), ".csv", sep = "")
        }
      }  else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format1_2 == "csv") {
          paste("Global_Estimates_von-", Sys.Date(), ".csv", sep = "")
        }
      }   else if (input$model == 'Gompertz Model'){
        if (input$format1_2 == "csv") {
          paste("Global_Estimates_gom-", Sys.Date(), ".csv", sep = "")
        }
      }
    },
    content = function(file) {
      if (input$model == 'Logistic Model') {
        if (input$format1_2 == "csv") {
          write.csv(data_global()$estimates, file)
        }
      } else if (input$model == 'Exponential Model') {
        if (input$format1_2 == "csv") {
          write.csv(data_global_expo()$estimates_expo, file)
        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format1_2 == "csv") {
          write.csv(data_global_theta()$estimates_theta, file)
        }
      } else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format1_2 == "csv") {
          write.csv(data_global_von()$estimates_von, file)
        }
      } else if (input$model == 'Gompertz Model') {
        if (input$format1_2 == "csv") {
          write.csv(data_global_gom()$estimates_gom, file)
        }
      }

    }
  )


  # Local Estimation


  data_local <- reactiveVal(NULL)
  data_local_expo <- reactiveVal(NULL)
  data_local_theta <- reactiveVal(NULL)
  data_local_von <- reactiveVal(NULL)
  data_local_gom <- reactiveVal(NULL)

  observeEvent(input$button1,{


    if (input$model == 'Logistic Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n2 <- input$trajectory
      data2 <- traj_plot(input$logi_k, input$logi_in, input$logi_r, input$logi_s, input$logi_rho, input$trajectory, input$t_range)$data
      x02 <- input$logi_in
      q <- length(t2)

      est_local = matrix(data = NA, nrow = length(t2)-P+1, ncol = 4)
      cov_local = list()

      fun_local_likelihood = function(param){
        r = param[1]
        K = param[2]
        sigma2  = param[3]
        rho = param[4]
        mu = K/(1 + (K/x02 - 1)*exp(-r*t2))

        cov_mat = matrix(data = NA, nrow = length(t2), ncol = length(t2))
        for (i in 1:length(t2)) {
          for (j in 1:length(t2)) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }



        # localized estimates
        mu_local = mu[ind:(ind + P-1)]
        data_local = data2[,ind:(ind + P - 1)]
        cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

        exponent = 0
        for(i in 1:n2){
          exponent = exponent + t((data_local[i,]) - (mu_local))%*%(solve(cov_mat_local))%*%((data_local[i,]) - (mu_local))
        }

        log_lik = - (n2*length(ind:(ind+P-1))/2)*log(2*pi) - (n2/2)*log(det(cov_mat_local)) - (1/2)*exponent
        return(-log_lik)
      }




      for(j in 1:(length(t2)-P+1)){
        ind = j
        param_0 = c(0.3, input$logi_k, 2,0.5)
        param = c(param_0)
        out2 = optim(param_0, fun_local_likelihood,method = "L-BFGS-B",
                     hessian = TRUE, lower = c(0.001, input$logi_k-(input$logi_k*0.2), 0.0001, -0.001),
                     upper = c(3, input$logi_k+(input$logi_k*0.2), 10, 0.999))
        est_local[ind, ] =  out2$par
        var_cov_local=solve(out2$hessian)
        cov_local[[j]]=var_cov_local

      }
      colnames(est_local) <- c("r", "K", "sigma2", "rho")
      result <- list(estimates_local = est_local, cov_matrix_local = cov_local)
      data_local(result)


      show_toast('Progress',
                 text = 'Done, We Appreciate Your Patience',
                 timer = 5000,
                 position = "center",
                 type = 'success')

      output$table2 <- renderTable({
        # Create a data frame from the covariance matrix
        est_local_df <- as.data.frame(data_local()$estimates_local)

        # Rename the columns
        colnames(est_local_df) <- c('r', 'K', 'Sigma', 'rho')

        # Format the numeric columns to display 5 decimal places
        est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

        head(est_local_df,8)
      }, width = "100%", align = 'c')


      output$plot_1 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local()$estimates_local[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')
      })

      output$plot_2 <- renderPlot({
        plot(1:(length(t2)-P+1),  data_local()$estimates_local[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(K(Delta~t))[MLE]),
             ylab = expression(hat(K)),xlab = 'Time Points')
      })

      output$plot_3 <- renderPlot({
        plot(1:(length(t2)-P+1),  data_local()$estimates_local[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')
      })

      output$plot_4 <- renderPlot({
        plot(1:(length(t2)-P+1),  data_local()$estimates_local[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')
      })




    } else if (input$model == 'Exponential Model'){


      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n2 <- input$trajectory
      data2 <- traj_plot_expo(input$expo_in, input$expo_r, input$expo_s, input$expo_rho, input$trajectory, input$t_range)$data
      x02 <- input$expo_in
      q <- length(t2)

      est_local_e <- matrix(data = NA, nrow = length(t2)-P+1, ncol = 3)
      cov_local_e <- list()

      fun_local_likelihood = function(param){
        r = param[1]
        sigma2  = param[2]
        rho = param[3]
        mu = x02*exp(r*t2) # Exponential mean function

        cov_mat = matrix(data = NA, nrow = length(t2), ncol = length(t2))
        for (i in 1:length(t2)) {
          for (j in 1:length(t2)) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }


        # localized estimates
        mu_local = mu[ind:(ind + P-1)]
        data_local = data2[,ind:(ind + P - 1)]
        cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

        exponent = 0
        for(i in 1:n2){
          exponent = exponent + t(data_local[i,] - mu_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_local)
        }





        log_lik = - (n2*length(ind:(ind+P-1))/2)*log(2*pi) - (n2/2)*log(det(cov_mat_local)) - (1/2)*exponent
        return(-log_lik)
      }



      for(j in 1:(length(t2)-P+1)){
        ind = j
        param_0 = c(input$expo_r, 2,0.5)
        out2 = optim(param_0, fun_local_likelihood, method = "L-BFGS-B",
                     hessian = TRUE, lower = c(0.001, 0.0001, -0.001),
                     upper = c(3, 10, 0.999))
        est_local_e[ind, ] =  out2$par
        var_cov_local_e=solve(out2$hessian)
        cov_local_e[[j]]=var_cov_local_e
      }
      colnames(est_local_e) <- c("r", "sigma2", "rho")

      result_expo <- list(estimates_local_expo = est_local_e, cov_matrix_local_expo = cov_local_e)
      data_local_expo(result_expo)



      output$table2 <- renderTable({
        # Create a data frame from the covariance matrix
        est_local_df <- as.data.frame(data_local_expo()$estimates_local_expo)

        # Rename the columns
        colnames(est_local_df) <- c('r', 'Sigma', 'rho')

        # Format the numeric columns to display 5 decimal places
        est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

        head(est_local_df,8)
      }, width = "100%", align = 'c')



      output$plot_1 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_expo()$estimates_local_expo[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')
      })

      output$plot_2 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_expo()$estimates_local_expo[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')
      })

      output$plot_3 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_expo()$estimates_local_expo[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')
      })

      output$plot_4 <- renderPlot({

      })


    } else if (input$model == 'Theta - Logistic Model') {

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n2 <- input$trajectory
      data2 <- traj_plot_theta(input$theta_k, input$theta_in,  input$theta_th,input$theta_r, input$theta_s, input$theta_rho,input$trajectory, input$t_range)$data
      x02 <- input$theta_in
      q <- length(t2)

      est_local = matrix(data = NA, nrow = length(t2)-P+1, ncol = 5)
      cov_local = list()

      fun_local_likelihood = function(param){
        r <- param[1]
        K <- param[2]
        theta <- param[3]
        sigma2  <- param[4]
        rho <- param[5]
        mu_theta <- K/(1 + ((K/x02)^theta - 1)*exp(-r*theta*t2))
        q <- length(t2)

        cov_mat = matrix(data = NA, nrow = length(t2), ncol = length(t2))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }

        # localized estimates
        mu_local = mu_theta[ind:(ind + P-1)]
        data_local = data2[,ind:(ind + P-1)]
        cov_mat_local = cov_mat[ind:(ind+ P-1), ind:(ind+ P-1)]

        exponent = 0
        for(i in 1:n2){
          exponent = exponent + t(data_local[i,] - mu_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_local)
        }

        log_lik = - (n2*length(ind:(ind+ P-1))/2)*log(2*pi) - (n2/2)*log(det(cov_mat_local)) - (1/2)*exponent
        return(-log_lik)
      }

      for(j in 1:(q-P+1)){
        ind = j
        param_0 <- c(0.3, input$theta_k, input$theta_th,2,0.5)
        out1 <- optim(param_0, fun_local_likelihood, method = "L-BFGS-B",
                      hessian = TRUE, lower = c(0.001, input$theta_k-(input$theta_k*0.2), 0.1, 0.0001, -0.001),
                      upper = c(3, input$theta_k+(input$theta_k*0.2), 2, 10, 0.999))
        est_local[j, ] <-  out1$par
        var_cov_local <- solve(out1$hessian)
        cov_local[[j]] <- var_cov_local
      }
      colnames(est_local) = c("r", "K", "theta", "sigma2", "rho")

      result_theta <- list(estimates_local_theta = est_local, cov_matrix_local_theta = cov_local)
      data_local_theta(result_theta)


      output$table2 <- renderTable({
        # Create a data frame from the covariance matrix
        est_local_df <- as.data.frame(data_local_theta()$estimates_local_theta)

        # Rename the columns
        colnames(est_local_df) <- c("r", "K", "theta", "sigma2", "rho")

        # Format the numeric columns to display 5 decimal places
        est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

        head(est_local_df,8)
      }, width = "100%", align = 'c')


      output$plot_1 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')
      })

      output$plot_2 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(K(Delta~t))[MLE]),
             ylab = expression(hat(K)),xlab = 'Time Points')
      })

      output$plot_3 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(theta(Delta~t))[MLE]),
             ylab = expression(hat(theta)),xlab = 'Time Points')
      })


      output$plot_4 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')
      })

      output$plot_5 <- renderPlot({
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,5], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')
      })


    } else if (input$model == 'Von-Bertallanfy Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n2 <- input$trajectory
      data2 <- traj_plot_von(input$von_k, input$von_in, input$von_r, input$von_s, input$von_rho, input$trajectory, input$t_range)$data
      x02 <- input$von_in
      q <- length(t2)

      est_local = matrix(data = NA, nrow = length(t2)-P+1, ncol = 4)
      cov_local = list()

      fun_local_likelihood = function(param){
        r = param[1]
        K = param[2]
        sigma2  = param[3]
        rho = param[4]
        mu = K*((1 + ((x02/K)^(1/3) - 1)*exp(-r*(t2/3)))^3)

        cov_mat = matrix(data = NA, nrow = length(t2), ncol = length(t2))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }

        # localized estimates
        mu_local = mu[ind:(ind + P-1)]
        data_local = data2[,ind:(ind + P - 1)]
        cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

        exponent = 0
        for(i in 1:n2){
          exponent = exponent + t(data_local[i,] - mu_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_local)
        }

        log_lik = - (n2*length(ind:(ind+P-1))/2)*log(2*pi) - (n2/2)*log(det(cov_mat_local)) - (1/2)*exponent
        return(-log_lik)
      }


      for(j in 1:(q-P+1)){
        ind = j
        param_0 = c(0.8, input$von_k, 2,0.5)
        out2 = optim(param_0, fun_local_likelihood, method = "L-BFGS-B",
                     hessian = TRUE, lower = c(0.001, input$von_k-(input$von_k*0.2), 0.0001, -0.001),
                     upper = c(3, input$von_k+(input$von_k*0.2), 10, 0.999))
        est_local[j, ] =  out2$par
        var_cov_local=solve(out2$hessian)
        cov_local[[j]]=var_cov_local
      }
      colnames(est_local) = c("r", "K", "sigma2", "rho")

      result <- list(estimates_local_von = est_local, cov_matrix_local_von = cov_local)
      data_local_von(result)



      output$table2 <- renderTable({
        # Create a data frame from the covariance matrix
        est_local_df <- as.data.frame(data_local_von()$estimates_local_von)

        # Rename the columns
        colnames(est_local_df) <- c('r', 'K', 'Sigma', 'rho')

        # Format the numeric columns to display 5 decimal places
        est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

        head(est_local_df,8)
      }, width = "100%", align = 'c')


      output$plot_1 <- renderPlot({
        plot(1:(q-P+1), data_local_von()$estimates_local_von[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')
      })

      output$plot_2 <- renderPlot({
        plot(1:(q-P+1), data_local_von()$estimates_local_von[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(K(Delta~t))[MLE]),
             ylab = expression(hat(K)),xlab = 'Time Points')
      })

      output$plot_3 <- renderPlot({
        plot(1:(q-P+1), data_local_von()$estimates_local_von[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')
      })

      output$plot_4 <- renderPlot({
        plot(1:(q-P+1), data_local_von()$estimates_local_von[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')
      })





    } else if (input$model == 'Gompertz Model'){


      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      n2 <- input$trajectory
      data2 <- traj_plot_gom(input$gom_b, input$gom_in, input$gom_c, input$gom_s, input$gom_rho, input$trajectory, input$t_range)$data
      x02 <- input$gom_in
      q <- length(t2)

      est_local = matrix(data = NA, nrow = length(t2)-P+1, ncol = 4)
      cov_local = list()

      fun_local_likelihood = function(param){
        c = param[1]
        b = param[2]
        sigma2  = param[3]
        rho = param[4]
        mu = x02*exp((b/c)*(1-exp(-c*t2)))

        cov_mat = matrix(data = NA, nrow = length(t2), ncol = length(t2))
        for (i in 1:q) {
          for (j in 1:q) {
            cov_mat[i,j] = sigma2*rho^(abs(i-j))
          }
        }

        # localized estimates
        mu_local = mu[ind:(ind + P-1)]
        data_local = data2[,ind:(ind + P-1)]
        cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

        exponent = 0
        for(i in 1:n2){
          exponent = exponent + t(data_local[i,] - mu_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_local)
        }

        log_lik = - (n2*length(ind:(ind+P-1))/2)*log(2*pi) - (n2/2)*log(det(cov_mat_local)) - (1/2)*exponent
        return(-log_lik)
      }



      for(j in 1:(q-P+1)){
        ind = j
        param_0 <- c(input$gom_c, input$gom_b, 2,0.5)
        out2 <- optim(param_0, fun_local_likelihood, method = "L-BFGS-B",
                      hessian = TRUE, lower = c(0.05, 0.05, 0.005, -0.001),
                      upper = c(5, 10, 10, 0.999))
        est_local[ind, ] <-  out2$par
        var_cov_local <- solve(out2$hessian)
        cov_local[[j]] <- var_cov_local
      }

      colnames(est_local) <- c("c", "b", "sigma2", "rho")
      result <- list(estimates_local_gom = est_local, cov_matrix_local_gom = cov_local)
      data_local_gom(result)



      output$table2 <- renderTable({
        # Create a data frame from the covariance matrix
        est_local_df <- as.data.frame(data_local_gom()$estimates_local_gom)

        # Rename the columns
        colnames(est_local_df) <-  c("c", "b", "sigma2", "rho")

        # Format the numeric columns to display 5 decimal places
        est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

        head(est_local_df,8)
      }, width = "100%", align = 'c')


      output$plot_1 <- renderPlot({

        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(c(Delta~t))[MLE]),
             ylab = expression(hat(c)),xlab = 'Time Points')

      })

      output$plot_2 <- renderPlot({
        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(b(Delta~t))[MLE]),
             ylab = expression(hat(b)),xlab = 'Time Points')

      })

      output$plot_3 <- renderPlot({
        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')

      })

      output$plot_4 <- renderPlot({
        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')
      })



    }})


  # covariance button


  observeEvent(input$Cov_M,{

    if (input$model == 'Logistic Model'){
      showNotification('Processing...',duration=3,type = 'warning')
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      P <- input$window

      showModal(
        modalDialog(
          renderTable({
            # Create a data frame from the covariance matrix
            cov_local_df <- as.data.frame(data_local()$cov_matrix_local[[length(t2)-P+1]])

            # Rename the columns
            colnames(cov_local_df) <- c('r', 'K', 'Sigma', 'rho')

            # Format the numeric columns to display 5 decimal places
            cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            cov_local_df
          }, width = "100%", align = 'c',height = 10000)
          ,easyClose = TRUE,
          size = 'l',
          footer = downloadButton('cov_data','Covariance Data'))
      )



    } else if (input$model == 'Exponential Model'){
      showNotification('Processing...',duration=3,type = 'warning')

      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      P <- input$window

      showModal(
        modalDialog(
          renderTable({
            # Create a data frame from the covariance matrix
            cov_local_df <- as.data.frame(data_local_expo()$cov_matrix_local_expo[[length(t2)-P+1]])

            # Rename the columns
            colnames(cov_local_df) <- c('r', 'Sigma', 'rho')

            # Format the numeric columns to display 5 decimal places
            cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            cov_local_df
          }, width = "100%", align = 'c',height = 10000)
          ,easyClose = TRUE,
          size = 'l',
          footer = downloadButton('cov_data','Covariance Data'))
      )


    } else if (input$model == 'Theta - Logistic Model'){
      showNotification('Processing...',duration=3,type = 'warning')
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      P <- input$window

      showModal(
        modalDialog(
          renderTable({
            # Create a data frame from the covariance matrix
            cov_local_df <- as.data.frame(data_local_theta()$cov_matrix_local_theta[[length(t2)-P+1]])

            # Rename the columns
            colnames(cov_local_df) <-  c('r', 'K','Theta', 'Sigma', 'rho')

            # Format the numeric columns to display 5 decimal places
            cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            cov_local_df
          }, width = "100%", align = 'c',height = 10000)
          ,easyClose = TRUE,
          size = 'l',
          footer = downloadButton('cov_data','Covariance Data'))
      )

    }
    else if (input$model == 'Von-Bertallanfy Model'){
      showNotification('Processing...',duration=3,type = 'warning')

      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      P <- input$window

      showModal(
        modalDialog(
          renderTable({
            # Create a data frame from the covariance matrix
            cov_local_df <- as.data.frame(data_local_von()$cov_matrix_local_von[[length(t2)-P+1]])

            # Rename the columns
            colnames(cov_local_df) <-  c('r', 'K', 'Sigma', 'rho')

            # Format the numeric columns to display 5 decimal places
            cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            cov_local_df
          }, width = "100%", align = 'c',height = 10000)
          ,easyClose = TRUE,
          size = 'l',
          footer = downloadButton('cov_data','Covariance Data'))
      )

    } else if (input$model == 'Gompertz Model'){
      showNotification('Processing...',duration=3,type = 'warning')
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      P <- input$window

      showModal(
        modalDialog(
          renderTable({
            # Create a data frame from the covariance matrix
            cov_local_df <- as.data.frame(data_local_gom()$cov_matrix_local_gom[[length(t2)-P+1]])

            # Rename the columns
            colnames(cov_local_df) <-  c("c", "b", "Sigma", "Rho")

            # Format the numeric columns to display 5 decimal places
            cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            cov_local_df
          }, width = "100%", align = 'c',height = 10000)
          ,easyClose = TRUE,
          size = 'l',
          footer = downloadButton('cov_data','Covariance Data'))
      )

    }})
  # download handler for covariance

  output$local_est <- downloadHandler(
    filename = function() {
      if (input$model == 'Logistic Model') {
        if (input$format2 == "csv") {
          paste("Local_Estimates_logi-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Exponential Model') {
        if (input$format2 == "csv") {
          paste("Local_Estimates_expo-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format2 == "csv") {
          paste("Local_Estimates_theta-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format2 == "csv") {
          paste("Local_Estimates_von-", Sys.Date(), ".csv", sep = "")
        }
      } else if (input$model == 'Gompertz Model'){
        if (input$format2 == "csv") {
          paste("Local_Estimates_gom-", Sys.Date(), ".csv", sep = "")
        }
      }
    },
    content = function(file) {
      if (input$model == 'Logistic Model') {
        if (input$format2 == "csv") {
          write.csv(data_local()$estimates_local, file)
        }
      } else if (input$model == 'Exponential Model') {
        if (input$format2 == "csv") {
          write.csv(data_local_expo()$estimates_local_expo, file)
        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format2 == "csv") {
          write.csv(data_local_theta()$estimates_local_theta, file)
        }
      } else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format2 == "csv") {
          write.csv(data_local_von()$estimates_local_von, file)
        }
      }  else if (input$model == 'Gompertz Model'){
        if (input$format2 == "csv") {
          write.csv(data_local_gom()$estimates_local_gom, file)
        }
      }
    }
  )


  output$cov_data <- downloadHandler(
    if (input$model == 'Logistic Model'){
      filename = function() {
        paste("covariance_data_logi.csv")
      }} else if (input$model == 'Exponential Model'){
        filename = function() {
          paste("covariance_data_expo.csv")
        }} else if (input$model == 'Theta - Logistic Model'){
          filename = function() {
            paste("covariance_data_theta.csv")
          }
        } else if (input$model == 'Von-Bertallanfy Model'){
          filename = function() {
            paste("covariance_data_von.csv")
          }
        } else if (input$model == 'Gompertz Model'){
          filename = function() {
            paste("covariance_data_von.csv")
          }
        },
    if (input$model == 'Logistic Model'){
      content = function(file) {
        write.csv( data_local()$cov_matrix_local, file)}

    } else if (input$model == 'Exponential Model'){
      content = function(file) {
        write.csv( data_local_expo()$cov_matrix_local_expo, file)}

    } else if (input$model == 'Theta - Logistic Model'){
      content = function(file) {
        write.csv( data_local_theta()$cov_matrix_local_theta, file)}
    }  else if (input$model == 'Von-Bertallanfy Model'){
      content = function(file) {
        write.csv( data_local_von()$cov_matrix_local_von, file)}
    }  else if (input$model == 'Gompertz Model'){
      content = function(file) {
        write.csv( data_local_gom()$cov_matrix_local_gom, file)}
    })

  #zooming table

  observeEvent(input$zoom3,{
    if (input$model == 'Logistic Model'){
      showNotification('Processing...',duration=3,type = 'warning')

      showModal(
        modalDialog(
          paste('Local Estimation Table (Complete)'),
          br(),
          renderTable({
            # Create a data frame from the covariance matrix
            est_local_df <- as.data.frame(data_local()$estimates_local)

            # Rename the columns
            colnames(est_local_df) <- c('r', 'K', 'Sigma', 'rho')

            # Format the numeric columns to display 5 decimal places
            est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            est_local_df
          }, width = "100%", align = 'c',height = 10000)
        )
      )

    } else if (input$model =='Exponential Model'){

      showNotification('Processing...',duration=3,type = 'warning')

      showModal(
        modalDialog(
          paste('Local Estimation Table (Complete)'),
          br(),
          renderTable({
            # Create a data frame from the covariance matrix
            est_local_df <- as.data.frame(data_local_expo()$estimates_local_expo)

            # Rename the columns
            colnames(est_local_df) <- c('r', 'Sigma', 'rho')

            # Format the numeric columns to display 5 decimal places
            est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            est_local_df
          }, width = "100%", align = 'c',height = 10000)
        )
      )


    } else if (input$model == 'Theta - Logistic Model'){
      showNotification('Processing...',duration=3,type = 'warning')

      showModal(
        modalDialog(
          paste('Local Estimation Table (Complete)'),
          br(),
          renderTable({
            # Create a data frame from the covariance matrix
            est_local_df <- as.data.frame(data_local_theta()$estimates_local_theta)

            # Rename the columns
            colnames(est_local_df) <- c("r", "K", "theta", "sigma2", "rho")

            # Format the numeric columns to display 5 decimal places
            est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            est_local_df
          }, width = "100%", align = 'c',height = 10000)
        )
      )
    } else if (input$model == 'Von-Bertallanfy Model'){
      showNotification('Processing...',duration=3,type = 'warning')

      showModal(
        modalDialog(
          paste('Local Estimation Table (Complete)'),
          br(),
          renderTable({
            # Create a data frame from the covariance matrix
            est_local_df <- as.data.frame(data_local_von()$estimates_local_von)

            # Rename the columns
            colnames(est_local_df) <- c('r','K', 'Sigma', 'rho')

            # Format the numeric columns to display 5 decimal places
            est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            est_local_df
          }, width = "100%", align = 'c',height = 10000)
        )
      )
    } else if (input$model == 'Gompertz Model'){
      showNotification('Processing...',type = 'warning',duration=3)

      showModal(
        modalDialog(
          paste('Local Estimation Table (Complete)'),
          br(),
          renderTable({
            # Create a data frame from the covariance matrix
            est_local_df <- as.data.frame(data_local_gom()$estimates_local_gom)

            # Rename the columns
            colnames(est_local_df) <- c("c", "b", "Sigma", "Rho")

            # Format the numeric columns to display 5 decimal places
            est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

            est_local_df
          }, width = "100%", align = 'c',height = 10000)
        )
      )
    }})


  #zooming of estimates plot (local)

  observeEvent(input$zoom4,{

    if (input$model == 'Logistic Model'){
      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)

      showModal(modalDialog(

        renderPlot({
          plot(1:(length(t2)-P+1), data_local()$estimates_local[,1], type = "b", col = "red", lwd = 2,
               main = expression(widehat(r(Delta~t))[MLE]),
               ylab = expression(hat(r)),xlab = 'Time Points')
        }, height = 500),
        renderPlot({
          plot(1:(length(t2)-P+1),data_local()$estimates_local[,2], type = "b", col = "red", lwd = 2,
               main = expression(widehat(K(Delta~t))[MLE]),
               ylab = expression(hat(K)),xlab = 'Time Points')
        },height=500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    } else if (input$model == 'Exponential Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)

      showModal(modalDialog(

        renderPlot({
          plot(1:(length(t2)-P+1), data_local_expo()$estimates_local_expo[,1], type = "b", col = "red", lwd = 2,
               main = expression(widehat(r(Delta~t))[MLE]),
               ylab = expression(hat(r)),xlab = 'Time Points')
        }, height = 500),
        renderPlot({
          plot(1:(length(t2)-P+1), data_local_expo()$estimates_local_expo[,2], type = "b", col = "red", lwd = 2,
               main = expression(widehat(sigma^2)),
               ylab = expression(hat(sigma^2)),xlab = 'Time Points')
        }, height = 500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))

    } else if (input$model == 'Theta - Logistic Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)

      showModal(modalDialog(

        renderPlot({
          plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,1], type = "b", col = "red", lwd = 2,
               main = expression(widehat(r(Delta~t))[MLE]),
               ylab = expression(hat(r)),xlab = 'Time Points')
        },height = 500),

        renderPlot({
          plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,2], type = "b", col = "red", lwd = 2,
               main = expression(widehat(K(Delta~t))[MLE]),
               ylab = expression(hat(K)),xlab = 'Time Points')
        },height = 500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    }  else if (input$model == 'Von-Bertallanfy Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      q <- length(t2)

      showModal(modalDialog(


        renderPlot({
          plot(1:(q-P+1), data_local_von()$estimates_local_von[,1], type = "b", col = "red", lwd = 2,
               main = expression(widehat(r(Delta~t))[MLE]),
               ylab = expression(hat(r)),xlab = 'Time Points')

        },height = 500),

        renderPlot({
          plot(1:(q-P+1), data_local_von()$estimates_local_von[,2], type = "b", col = "red", lwd = 2,
               main = expression(widehat(K(Delta~t))[MLE]),
               ylab = expression(hat(K)),xlab = 'Time Points')

        },height = 500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    } else if (input$model == 'Gompertz Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      q <- length(t2)

      showModal(modalDialog(


        renderPlot({
          plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,1], type = "b", col = "red", lwd = 2,
               main = expression(widehat(c(Delta~t))[MLE]),
               ylab = expression(hat(c)),xlab = 'Time Points')


        },height = 500),

        renderPlot({
          plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,2], type = "b", col = "red", lwd = 2,
               main = expression(widehat(b(Delta~t))[MLE]),
               ylab = expression(hat(b)),xlab = 'Time Points')


        },height = 500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    }})


  observeEvent(input$zoom5,{
    if (input$model == 'Logistic Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)

      showModal(modalDialog(

        renderPlot({
          plot(1:(length(t2)-P+1), data_local()$estimates_local[,3], type = "b", col = "red", lwd = 2,
               main = expression(widehat(sigma^2)),
               ylab = expression(hat(sigma^2)),xlab = 'Time Points')
        }, height = 500),
        renderPlot({
          plot(1:(length(t2)-P+1), data_local()$estimates_local[,4], type = "b", col = "red", lwd = 2,
               main = expression(widehat(rho)),
               ylab = expression(hat(rho)),xlab = 'Time Points')
        },height=500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))
    } else if (input$model == 'Exponential Model'){

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)

      showModal(modalDialog(

        renderPlot({
          plot(1:(length(t2)-P+1), data_local_expo()$estimates_local_expo[,3], type = "b", col = "red", lwd = 2,
               main = expression(widehat(rho)),
               ylab = expression(hat(rho)),xlab = 'Time Points')
        },height=500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))

    } else if (input$model == 'Theta - Logistic Model') {

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)

      showModal(modalDialog(

        renderPlot({
          plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,3], type = "b", col = "red", lwd = 2,
               main = expression(widehat(theta(Delta~t))[MLE]),
               ylab = expression(hat(theta)),xlab = 'Time Points')
        },height = 500),


        renderPlot({
          plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,4], type = "b", col = "red", lwd = 2,
               main = expression(widehat(sigma^2)),
               ylab = expression(hat(sigma^2)),xlab = 'Time Points')
        },height = 500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))

    } else if (input$model == 'Von-Bertallanfy Model') {

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      q <- length(t2)

      showModal(modalDialog(

        renderPlot({
          plot(1:(q-P+1), data_local_von()$estimates_local_von[,3], type = "b", col = "red", lwd = 2,
               main = expression(widehat(sigma^2)),
               ylab = expression(hat(sigma^2)),xlab = 'Time Points')

        },height = 500),


        renderPlot({
          plot(1:(q-P+1), data_local_von()$estimates_local_von[,4], type = "b", col = "red", lwd = 2,
               main = expression(widehat(rho)),
               ylab = expression(hat(rho)),xlab = 'Time Points')
        },height = 500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))

    } else if (input$model == 'Gompertz Model') {

      P <- input$window
      t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
      q <- length(t2)

      showModal(modalDialog(

        renderPlot({
          plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,3], type = "b", col = "red", lwd = 2,
               main = expression(widehat(sigma^2)),
               ylab = expression(hat(sigma^2)),xlab = 'Time Points')


        },height = 500),


        renderPlot({
          plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,4], type = "b", col = "red", lwd = 2,
               main = expression(widehat(rho)),
               ylab = expression(hat(rho)),xlab = 'Time Points')
        },height = 500),
        easyClose = TRUE,
        size = "l",
        footer = NULL
      ))

    }})

  observeEvent(input$zoom6 ,{

    P <- input$window
    t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)

    showModal(modalDialog(

      renderPlot({
        plot(1:(length(t2)-P+1),data_local_theta()$estimates_local_theta[,5], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')
      },height = 500),

      easyClose = TRUE,
      size = "l",
      footer = NULL
    ))


  })

  # download handler for estimate plots (local)

  output$l_plot_d <- downloadHandler(
    filename = function() {
      if (input$model == 'Logistic Model') {
        if (input$format3 == "png") {
          paste("data_plot_local_logi-", Sys.Date(), ".png", sep = "")
        } else if (input$format3 == "pdf") {
          paste("data_plot_local_logi-", Sys.Date(), ".pdf", sep = "")
        }
      } else if (input$model == 'Exponential Model') {
        if (input$format3 == "png") {
          paste("data_plot_local_expo-", Sys.Date(), ".png", sep = "")
        } else if (input$format3 == "pdf") {
          paste("data_plot_local_expo-", Sys.Date(), ".pdf", sep = "")
        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format3 == "png") {
          paste("data_plot_local_theta-", Sys.Date(), ".png", sep = "")
        } else if (input$format3 == "pdf") {
          paste("data_plot_local_theta-", Sys.Date(), ".pdf", sep = "")
        }
      } else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format3 == "png") {
          paste("data_plot_local_von-", Sys.Date(), ".png", sep = "")
        } else if (input$format3 == "pdf") {
          paste("data_plot_local_von-", Sys.Date(), ".pdf", sep = "")
        }
      } else if (input$model == 'Gompertz Model'){
        if (input$format3 == "png") {
          paste("data_plot_local_gom-", Sys.Date(), ".png", sep = "")
        } else if (input$format3 == "pdf") {
          paste("data_plot_local_gom-", Sys.Date(), ".pdf", sep = "")
        }
      }
    },
    content = function(file) {
      if (input$model == 'Logistic Model') {
        # Plot for Logistic Model
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        if (input$format3 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format3 == "pdf") {
          pdf(file,  width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(length(t2) - P + 1), data_local()$estimates_local[, 1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')

        plot(1:(length(t2) - P + 1), data_local()$estimates_local[, 2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(K(Delta~t))[MLE]),
             ylab = expression(hat(K)),xlab = 'Time Points')

        dev.off()
      } else if (input$model == 'Exponential Model') {
        # Plot for Exponential Model
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        if (input$format3 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format3 == "pdf") {
          pdf(file,  width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(length(t2) - P + 1), data_local_expo()$estimates_local_expo[, 1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')

        plot(1:(length(t2) - P + 1), data_local_expo()$estimates_local_expo[, 2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')

        dev.off()
      } else if (input$model == 'Theta - Logistic Model'){
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        q <- length(t2)
        if (input$format3 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format3 == "pdf") {
          pdf(file, width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')

        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(K(Delta~t))[MLE]), ylab = expression(hat(K)),xlab = 'Time Points')

        dev.off()
      } else if (input$model == 'Von-Bertallanfy Model'){
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        q <- length(t2)
        if (input$format3 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format3 == "pdf") {
          pdf(file, width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(q-P+1), data_local_von()$estimates_local_von[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(r(Delta~t))[MLE]),
             ylab = expression(hat(r)),xlab = 'Time Points')

        plot(1:(q-P+1), data_local_von()$estimates_local_von[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(K(Delta~t))[MLE]),
             ylab = expression(hat(K)),xlab = 'Time Points')


        dev.off()
      } else if (input$model == 'Gompertz Model'){
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        q <- length(t2)
        if (input$format3 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format3 == "pdf") {
          pdf(file, width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,1], type = "b", col = "red", lwd = 2,
             main = expression(widehat(c(Delta~t))[MLE]),
             ylab = expression(hat(c)),xlab = 'Time Points')

        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,2], type = "b", col = "red", lwd = 2,
             main = expression(widehat(b(Delta~t))[MLE]),
             ylab = expression(hat(b)),xlab = 'Time Points')



        dev.off()
      }
    }
  )

  output$l_plot_d2 <- downloadHandler(

    filename = function() {
      if (input$model == 'Logistic Model') {
        if (input$format4 == "png") {
          paste("data_plot_local_logi-", Sys.Date(), ".png", sep = "")
        } else if (input$format4 == "pdf") {
          paste("data_plot_local_logi-", Sys.Date(), ".pdf", sep = "")
        }
      } else if (input$model == 'Exponential Model') {
        if (input$format4 == "png") {
          paste("data_plot_local_expo-", Sys.Date(), ".png", sep = "")
        } else if (input$format4 == "pdf") {
          paste("data_plot_local_expo-", Sys.Date(), ".pdf", sep = "")
        }
      } else if (input$model == 'Theta - Logistic Model'){
        if (input$format4 == "png") {
          paste("data_plot_local_theta-", Sys.Date(), ".png", sep = "")
        } else if (input$format4 == "pdf") {
          paste("data_plot_local_theta-", Sys.Date(), ".pdf", sep = "")
        }
      }  else if (input$model == 'Von-Bertallanfy Model'){
        if (input$format4 == "png") {
          paste("data_plot_local_von-", Sys.Date(), ".png", sep = "")
        } else if (input$format4 == "pdf") {
          paste("data_plot_local_von-", Sys.Date(), ".pdf", sep = "")
        }
      }  else if (input$model == 'Gompertz Model'){
        if (input$format4 == "png") {
          paste("data_plot_local_gom-", Sys.Date(), ".png", sep = "")
        } else if (input$format4 == "pdf") {
          paste("data_plot_local_gom-", Sys.Date(), ".pdf", sep = "")
        }
      }
    },
    content = function(file) {
      if (input$model == 'Logistic Model'){
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        if (input$format4 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format4 == "pdf") {
          pdf(file, width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(length(t2)-P+1), data_local()$estimates_local[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')
        plot(1:(length(t2)-P+1), data_local()$estimates_local[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')

        dev.off()

      } else if (input$model == 'Exponential Model'){

        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        if (input$format4 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format4 == "pdf") {
          pdf(file,  width = 11, height = 8.5)
        }

        # Set up a 1x2 grid for two plots side by side
        plot(1:(length(t2)-P+1), data_local_expo()$estimates_local_expo[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')
        dev.off()

      } else if (input$model == 'Theta - Logistic Model'){
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        if (input$format4 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format4 == "pdf") {
          pdf(file,  width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(theta(Delta~t))[MLE]),
             ylab = expression(hat(theta)),xlab = 'Time Points')
        plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)), ylab = expression(hat(sigma^2)),xlab = 'Time Points')

        dev.off()

      } else if (input$model == 'Von-Bertallanfy Model'){
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        q <- length(t2)
        if (input$format4 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format4 == "pdf") {
          pdf(file, width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(q-P+1), data_local_von()$estimates_local_von[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')

        plot(1:(q-P+1), data_local_von()$estimates_local_von[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')

        dev.off()

      } else if (input$model == 'Gompertz Model'){
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        q <- length(t2)
        if (input$format4 == "png") {
          png(file, width = 1100, height = 500)
        } else if (input$format4 == "pdf") {
          pdf(file, width = 11, height = 8.5)
        }

        par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,3], type = "b", col = "red", lwd = 2,
             main = expression(widehat(sigma^2)),
             ylab = expression(hat(sigma^2)),xlab = 'Time Points')

        plot(1:(q-P+1), data_local_gom()$estimates_local_gom[,4], type = "b", col = "red", lwd = 2,
             main = expression(widehat(rho)),
             ylab = expression(hat(rho)),xlab = 'Time Points')

        dev.off()

      }}
  )

  output$l_plot_d3 <- downloadHandler(
    filename = function() {
      if (input$model == 'Theta - Logistic Model') {
        if (input$format5 == "png") {
          paste("data_plot_local_theta-", Sys.Date(), ".png", sep = "")
        } else if (input$format5 == "pdf") {
          paste("data_plot_local_theta-", Sys.Date(), ".pdf", sep = "")
        }
      }
    },
    content = function(file) {
      if (input$model == 'Theta - Logistic Model') {
        P <- input$window
        t2 <- seq(from = input$t_range[1], to = input$t_range[2], by = 1)
        if (input$format5 == "png") {
          png(file, width = 1100, height = 500)
          plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,5], type = "b", col = "red", lwd = 2,
               main = expression(widehat(rho)),
               ylab = expression(hat(rho)),xlab = 'Time Points')
        } else if (input$format5 == "pdf") {
          pdf(file,  width = 11, height = 8.5)
          plot(1:(length(t2)-P+1), data_local_theta()$estimates_local_theta[,5], type = "b", col = "red", lwd = 2,
               main = expression(widehat(rho)),
               ylab = expression(hat(rho)),xlab = 'Time Points')
        }
        dev.off()
      }
    }
  )



  # <------------------------------- SIMULATED DATA (END) -------------------------------------->


  #<--------------------------------- REAL DATA (START) ---------------------------------------->


  observeEvent(input$upload, {
    # Check if a file has been uploaded
    if (is.null(input$upload)) {
      showNotification("Please upload a data file.", type = "error")
    } else {
      # Read the uploaded data


      if (tools::file_ext(input$upload$name)  == "txt") {
        # Read a TXT file
        data <- read.table(input$upload$datapath)
        data <- as.matrix(data)
      } else if (tools::file_ext(input$upload$name)  == "csv") {
        # Read a CSV file
        data <- read.csv(input$upload$datapath)
        data <- as.matrix(data)
      }

      # Check if all columns are numeric
      if (all(sapply(data, is.numeric))) {
        # Data is numeric, perform calculation
        output$data <- renderDataTable({
          if (!is.null(data)) {
            data.table::data.table(data)


          }
        })



        observeEvent(input$button2,{
          show_toast('Important',
                     text = 'Dont click this button twice, Computation May take little time',
                     timer = 11000,
                     position = "center",
                     type = 'success')
        })

        observeEvent(input$button2,{
          output$model_name2 <- renderText({
            selected_model <- input$model1

            if (selected_model == "Logistic Model") {
              return("Estimation using Logistic Model")
            } else if (selected_model == "Von-Bertallanfy Model") {
              return("Estimation using Von-Bertallanfy Model")
            } else if (selected_model == "Gompertz Model") {
              return("Estimation using Gompertz Model")
            } else {
              return("Please select a model.")
            }
          })
        })


        size_pro <- eventReactive(input$button2,{
          matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2)
        })

        output$data_p <- renderPlot({
          size_pro()
        })

        observeEvent(input$zoom7, {

          showModal(modalDialog(

            renderPlot({
              matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2)
            }, height = 600),
            easyClose = TRUE,
            size = "l",
            footer = NULL
          ))

        })

        output$r_data_down <- downloadHandler(
          filename = function() {

            if (input$format6 == "png") {
              paste("data_plot_real-", Sys.Date(), ".png", sep = "")
            } else if (input$format6 == "pdf") {
              paste("data_plot_real-", Sys.Date(), ".pdf", sep = "")
            }

          },
          content = function(file) {
            if (input$format6 == "png") {
              png(file)

              matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2)
              dev.off()
            } else if (input$format6 == "pdf") {
              pdf(file)
              matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2)
              dev.off()
            }}

        )

        # Global estimates Real

        # this place

        r_data_global_logi <- reactiveVal(NULL)
        r_data_global_von <- reactiveVal(NULL)
        r_data_global_gom <- reactiveVal(NULL)

        observeEvent(input$button2, {
          if (input$model1 == 'Logistic Model'){

            q <- ncol(data)
            t <- 1:q
            n <- nrow(data)
            r1 <- input$prob
            x0 <- mean(data[,1])
            K <- mean(data[,q])
            K1 <- K-K*0.3
            K2 <- K+K*0.3




            # Global Likelihood

            log_fun_likelihood = function(param_real){
              r = param_real[1]
              K = param_real[2]
              sigma2  = param_real[3]
              rho = param_real[4]
              mu_log = K/(1+(K/x0 - 1)*exp(-r*t))

              cov_mat = matrix(data = NA, nrow = q, ncol = q)
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }
              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data[i,]-mu_log)%*%(solve(cov_mat))%*%(data[i,]-mu_log)
              }
              log_lik = - (q/2)*log(2*pi) - (1/2)*log(det(cov_mat)) - (1/2)*exponent
              return(-log_lik)
            }
            param_rr <- c(r1, K, 3, 0.3)
            out_log <- optim(param_rr, log_fun_likelihood, method = "L-BFGS-B",
                             hessian = TRUE, lower = c(0.001, K1, 0.0001, -0.001),
                             upper = c(3, K2, 10, 0.999))
            cov_r <- data.frame(solve(out_log$hessian))
            colnames(cov_r) <- c('r','K','Sigma','Rho')

            # Create a data frame with parameter names
            result_df_log <- data.frame(r = out_log$par[1], K = out_log$par[2], sigma = out_log$par[3], rho = out_log$par[4])


            result_list_r <- list(estimates_real_logi = result_df_log, cov_matrix_real_logi = cov_r)
            r_data_global_logi(result_list_r)


            output$table1_r <- renderTable({
              result_df_format <- r_data_global_logi()$estimates_real_logi  # Create a copy of the data frame to format
              # Format numeric columns to display 5 decimal places
              numeric_columns <- sapply(result_df_format, is.numeric)
              result_df_format[, numeric_columns] <- lapply(result_df_format[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
              result_df_format
            }, width = "80%", align = 'c',height = 1000)


            output$covM_r <- renderTable({
              # Create a data frame from the covariance matrix
              cov_df_r <- as.data.frame(r_data_global_logi()$cov_matrix_real_logi)

              # Rename the columns
              colnames(cov_df_r) <- c('r', 'K', 'Sigma', 'rho')

              # Format the numeric columns to display 5 decimal places
              cov_df_r[, sapply(cov_df_r, is.numeric)] <- lapply(cov_df_r[, sapply(cov_df_r, is.numeric)], function(x) format(x, nsmall = 5))

              cov_df_r
            }, width = "80%", align = 'c',height= 1000)}

          else if (input$model1 == 'Von-Bertallanfy Model'){
            q <- ncol(data)
            t <- 1:q
            n <- nrow(data)
            r2 <- input$von_r1
            x0 <- mean(data[,1])
            K <- mean(data[,q])
            K1 <- K-K*0.3
            K2 <- K+K*0.3

            # Global Likelihood

            VB_fun_likelihood = function(param){
              r = param[1]
              K = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_VB = K*((1 + ((x0/K)^(1/3) - 1)*exp(-r*(t/3)))^3)

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }
              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data[i,]-mu_VB)%*%(solve(cov_mat))%*%(data[i,]-mu_VB)
              }


              log_lik = - (n*q/2)*log(2*pi) - (n/2)*log(det(cov_mat)) - (1/2)*exponent
              return(-log_lik)
            }
            param_0 <- c(r2, K, 2, 0.3)
            out_VB <- optim(param_0, VB_fun_likelihood, method = "L-BFGS-B",
                            hessian = TRUE, lower = c(0.001, K1, 0.0001, -0.001),
                            upper = c(3, K2, 3, 0.999))

            cov <- data.frame(solve(out_VB$hessian))
            colnames(cov) <- c('r','K','Sigma','Rho')

            # Create a data frame with parameter names
            result_df <- data.frame(r = out_VB$par[1], K = out_VB$par[2], sigma = out_VB$par[3], rho = out_VB$par[4])


            result_list <- list(estimates_real_von = result_df, cov_matrix_real_von = cov)
            r_data_global_von(result_list)


            output$table1_r <- renderTable({
              result_df_formatted <- r_data_global_von()$estimates_real_von  # Create a copy of the data frame to format
              # Format numeric columns to display 5 decimal places
              numeric_columns <- sapply(result_df_formatted, is.numeric)
              result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
              result_df_formatted
            }, width = "80%", align = 'c',height = 1000)


            output$covM_r <- renderTable({
              # Create a data frame from the covariance matrix
              cov_df <- as.data.frame(r_data_global_von()$cov_matrix_real_von)

              # Rename the columns
              colnames(cov_df) <- c('r', 'K', 'Sigma', 'rho')

              # Format the numeric columns to display 5 decimal places
              cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

              cov_df
            }, width = "80%", align = 'c',height= 1000)
          }     else if (input$model1 == 'Gompertz Model'){
            q <- ncol(data)
            t <- 1:q
            n <- nrow(data)
            b <- input$gom_b1
            c <- input$gom_c1
            x0 <- mean(data[,1])
            K <- mean(data[,q])
            K1 <- K-K*0.3
            K2 <- K+K*0.3

            # Global Likelihood

            gom_fun_likelihood = function(param){
              c = param[1]
              b = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_gom = x0*exp((b/c)*(1-exp(-c*t)))

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }
              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data[i,]-mu_gom)%*%(solve(cov_mat))%*%(data[i,]-mu_gom)
              }


              log_lik = - (n*q/2)*log(2*pi) - (n/2)*log(det(cov_mat)) - (1/2)*exponent
              return(-log_lik)
            }
            param_0 = c(c, b, 2, 0.3)
            out_gom = optim(param_0, gom_fun_likelihood, method = "L-BFGS-B",
                            hessian = TRUE, lower = c(0.01, 0.01, 0.001, -0.001),
                            upper = c(5, 10, 3, 0.999))

            cov <- data.frame(solve(out_gom$hessian))
            colnames(cov) <- c("c", "b", "Sigma", "Rho")

            # Create a data frame with parameter names
            result_gom <- data.frame(c = out_gom$par[1], b = out_gom$par[2], Sigma = out_gom$par[3], Rho = out_gom$par[4])


            result_list <- list(estimates_real_gom = result_gom, cov_matrix_real_gom = cov)
            r_data_global_gom(result_list)


            output$table1_r <- renderTable({
              result_df_formatted <- r_data_global_gom()$estimates_real_gom  # Create a copy of the data frame to format
              # Format numeric columns to display 5 decimal places
              numeric_columns <- sapply(result_df_formatted, is.numeric)
              result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
              result_df_formatted
            }, width = "80%", align = 'c',height = 1000)


            output$covM_r <- renderTable({
              # Create a data frame from the covariance matrix
              cov_df <- as.data.frame(r_data_global_gom()$cov_matrix_real_gom)

              # Rename the columns
              colnames(cov_df) <- c("c", "b", "Sigma", "Rho")

              # Format the numeric columns to display 5 decimal places
              cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

              cov_df
            }, width = "80%", align = 'c',height= 1000)
          }


          show_toast('Progress',
                     text = 'Global Estimates are Ready',
                     timer = 5000,
                     position = "center",
                     type = 'success')
        })




        observeEvent(input$zoom8, {
          if (input$model1 == 'Logistic Model'){
            showModal(modalDialog(
              paste('Global Estimation Table'),
              br(),
              renderTable({
                result_df_formatted <- r_data_global_logi()$estimates_real_logi  # Create a copy of the data frame to format
                # Format numeric columns to display 5 decimal places
                numeric_columns <- sapply(result_df_formatted, is.numeric)
                result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
                result_df_formatted
              }, width = "100%", align = 'c',height = 1000),
              br(),
              br(),
              paste('Covariance Matrix'),
              br(),
              renderTable({
                # Create a data frame from the covariance matrix
                cov_df <- as.data.frame(r_data_global_logi()$cov_matrix_real_logi)

                # Rename the columns
                colnames(cov_df) <- c('r', 'K', 'Sigma', 'rho')

                # Format the numeric columns to display 5 decimal places
                cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

                cov_df
              }, width = "100%", align = 'c',height= 1000),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          }


          else if (input$model1 == 'Von-Bertallanfy Model'){

            showModal(modalDialog(
              paste('Global Estimation Table'),
              br(),
              renderTable({
                result_df_formatted <- r_data_global_von()$estimates_real_von  # Create a copy of the data frame to format
                # Format numeric columns to display 5 decimal places
                numeric_columns <- sapply(result_df_formatted, is.numeric)
                result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
                result_df_formatted
              }, width = "100%", align = 'c',height = 1000),
              br(),
              br(),
              paste('Covariance Matrix'),
              br(),
              renderTable({
                # Create a data frame from the covariance matrix
                cov_df <- as.data.frame(r_data_global_von()$cov_matrix_real_von)

                # Rename the columns
                colnames(cov_df) <- c('r', 'K', 'Sigma', 'rho')

                # Format the numeric columns to display 5 decimal places
                cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

                cov_df
              }, width = "100%", align = 'c',height= 1000),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))


          } else if (input$model1 == 'Gompertz Model'){

            showModal(modalDialog(
              paste('Global Estimation Table'),
              br(),
              renderTable({
                result_df_formatted <- r_data_global_gom()$estimates_real_gom  # Create a copy of the data frame to format
                # Format numeric columns to display 5 decimal places
                numeric_columns <- sapply(result_df_formatted, is.numeric)
                result_df_formatted[, numeric_columns] <- lapply(result_df_formatted[, numeric_columns], function(x) formatC(x, format = "f", digits = 5))
                result_df_formatted
              }, width = "100%", align = 'c',height = 1000),
              br(),
              br(),
              paste('Covariance Matrix'),
              br(),
              renderTable({
                # Create a data frame from the covariance matrix
                cov_df <- as.data.frame(r_data_global_gom()$cov_matrix_real_gom)

                # Rename the columns
                colnames(cov_df) <- c("c", "b", "Sigma", "Rho")

                # Format the numeric columns to display 5 decimal places
                cov_df[, sapply(cov_df, is.numeric)] <- lapply(cov_df[, sapply(cov_df, is.numeric)], function(x) format(x, nsmall = 5))

                cov_df
              }, width = "100%", align = 'c',height= 1000),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))


          }
        })

        output$c_mat_r <- downloadHandler(
          filename = function() {
            if (input$model1 == 'Logistic Model' && input$format7_1 == 'csv') {
              paste("Cov_Matrix_logi-", Sys.Date(), ".csv", sep = "")
            }  else if (input$model1 == 'Von-Bertallanfy Model' && input$format7_1 == 'csv'){
              paste("Cov_Matrix_von-", Sys.Date(), ".csv", sep = "")
            }  else if (input$model1 == 'Gompertz Model' && input$format7_1 == 'csv'){
              paste("Cov_Matrix_gom-", Sys.Date(), ".csv", sep = "")
            }
          },
          content = function(file) {
            if (input$model1 == 'Logistic Model' && input$format7_1 == 'csv') {
              write.csv(r_data_global_logi()$cov_matrix_real_logi, file)
            } else if (input$model1 == 'Von-Bertallanfy Model' && input$format7_1 == 'csv'){
              write.csv(r_data_global_von()$cov_matrix_real_von, file)
            }  else if (input$model1 == 'Gompertz Model' && input$format7_1 == 'csv'){
              write.csv( r_data_global_gom()$cov_matrix_real_gom, file)
            }
          }
        )

        output$down_global_r <- downloadHandler(
          filename = function() {
            if (input$model1 == 'Logistic Model') {
              if (input$format7_2 == "csv") {
                paste("Global_Estimates_logi-", Sys.Date(), ".csv", sep = "")
              }
            }  else if (input$model1 == 'Von-Bertallanfy Model'){
              if (input$format7_2 == "csv") {
                paste("Global_Estimates_von-", Sys.Date(), ".csv", sep = "")
              }
            }   else if (input$model1 == 'Gompertz Model'){
              if (input$format7_2 == "csv") {
                paste("Global_Estimates_gom-", Sys.Date(), ".csv", sep = "")
              }
            }
          },
          content = function(file) {
            if (input$model1 == 'Logistic Model') {
              if (input$format7_2 == "csv") {
                write.csv(r_data_global_logi()$estimates_real_logi, file)
              }
            } else if (input$model1 == 'Von-Bertallanfy Model'){
              if (input$format7_2 == "csv") {
                write.csv(r_data_global_von()$estimates_real_von, file)
              }
            } else if (input$model1 == 'Gompertz Model') {
              if (input$format7_2 == "csv") {
                write.csv(r_data_global_gom()$estimates_real_gom, file)
              }
            }

          }
        )
        generatePlot2 <- eventReactive(input$button2, {
          if (input$model1 == 'Logistic Model'){
            #com_data <- melt(data.frame(x=1:ncol(data),t(data)), id='x')
            K <- r_data_global_logi()$estimates_real_logi[,2]
            x0 <- mean(data[,1])
            r <- r_data_global_logi()$estimates_real_logi[,1]
            q <- ncol(data)
            t <- 1:q
            mu_log <- K/(1+(K/x0 - 1)*exp(-r*t))
            lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")


          } else if (input$model1 == 'Von-Bertallanfy Model'){

            K <- r_data_global_von()$estimates_real_von[,2]
            x0 <- mean(data[,1])
            r <- r_data_global_von()$estimates_real_von[,1]
            q <- ncol(data)
            t <- 1:q
            mu_log <- K*((1 + ((x0/K)^(1/3) - 1)*exp(-r*(t/3)))^3)
            lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")

          } else if (input$model1 == 'Gompertz Model'){


            b <- r_data_global_gom()$estimates_real_gom[,2]
            x0 <- mean(data[,1])
            c <- r_data_global_gom()$estimates_real_gom[,1]
            q <- ncol(data)
            t <- 1:q
            mu_log <- x0*exp((b/c)*(1-exp(-c*t)))
            lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")

          }
        })

        output$com_plot <- renderPlot({
          matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2,col = 'grey')
          if (input$model1 == 'Logistic Model'){
            generatePlot2()

          } else if (input$model1 == 'Von-Bertallanfy Model'){
            generatePlot2()
          } else if (input$model1 == 'Gompertz Model'){
            generatePlot2()
          }
        })

        observeEvent(input$zoom9,{

          if (input$model1 == 'Logistic Model'){
            P <- input$window2
            q <- ncol(data)
            t <- 1:q

            showModal(modalDialog(

              renderPlot({
                matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2,col = 'grey')
                K <- r_data_global_logi()$estimates_real_logi[,2]
                x0 <- mean(data[,1])
                r <- r_data_global_logi()$estimates_real_logi[,1]
                q <- ncol(data)
                t <- 1:q
                mu_log <- K/(1+(K/x0 - 1)*exp(-r*t))
                lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")
              }, height = 500),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          }   else if (input$model1 == 'Von-Bertallanfy Model'){

            P <- input$window2
            q <- ncol(data)
            t <- 1:q

            showModal(modalDialog(


              renderPlot({
                matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2,col = 'grey')
                K <- r_data_global_logi()$estimates_real_logi[,2]
                x0 <- mean(data[,1])
                r <- r_data_global_logi()$estimates_real_logi[,1]
                q <- ncol(data)
                t <- 1:q
                mu_log <- K/(1+(K/x0 - 1)*exp(-r*t))
                lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")
              },height = 500),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          } else if (input$model1 == 'Gompertz Model'){

            P <- input$window2
            q <- ncol(data)
            t <- 1:q

            showModal(modalDialog(


              renderPlot({
                matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2,col = 'grey')
                K <- r_data_global_logi()$estimates_real_logi[,2]
                x0 <- mean(data[,1])
                r <- r_data_global_logi()$estimates_real_logi[,1]
                q <- ncol(data)
                t <- 1:q
                mu_log <- K/(1+(K/x0 - 1)*exp(-r*t))
                lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")
              },height = 500),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          }


        })

        output$com_plot_r <- downloadHandler(

          filename = function() {
            if (input$model1 == 'Logistic Model') {
              if (input$format8 == "png") {
                paste("com_logi-", Sys.Date(), ".png", sep = "")
              } else if (input$format8 == "pdf") {
                paste("com_logi-", Sys.Date(), ".pdf", sep = "")
              }
            } else if (input$model1 == 'Von-Bertallanfy Model'){
              if (input$format8 == "png") {
                paste("com_von-", Sys.Date(), ".png", sep = "")
              } else if (input$format8 == "pdf") {
                paste("com_von-", Sys.Date(), ".pdf", sep = "")
              }
            } else if (input$model == 'Gompertz Model'){
              if (input$format8 == "png") {
                paste("com_gom-", Sys.Date(), ".png", sep = "")
              } else if (input$format8 == "pdf") {
                paste("com_gom-", Sys.Date(), ".pdf", sep = "")
              }
            }
          },
          content = function(file) {
            if (input$model1 == 'Logistic Model') {


              if (input$format8 == "png") {
                png(file, width = 1100, height = 500)
              } else if (input$format8 == "pdf") {
                pdf(file,  width = 11, height = 8.5)
              }
              matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2,col = 'grey')
              K <- r_data_global_logi()$estimates_real_logi[,2]
              x0 <- mean(data[,1])
              r <- r_data_global_logi()$estimates_real_logi[,1]
              q <- ncol(data)
              t <- 1:q
              mu_log <- K/(1+(K/x0 - 1)*exp(-r*t))
              lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")
              dev.off()
            } else if (input$model1 == 'Von-Bertallanfy Model'){

              if (input$format8 == "png") {
                png(file, width = 1100, height = 500)
              } else if (input$format3 == "pdf") {
                pdf(file, width = 8, height = 8.5)
              }
              matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2,col = 'grey')
              K <- r_data_global_logi()$estimates_real_logi[,2]
              x0 <- mean(data[,1])
              r <- r_data_global_logi()$estimates_real_logi[,1]
              q <- ncol(data)
              t <- 1:q
              mu_log <- K/(1+(K/x0 - 1)*exp(-r*t))
              lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")
              dev.off()
            } else if (input$model1 == 'Gompertz Model'){

              if (input$format8 == "png") {
                png(file, width = 1100, height = 500)
              } else if (input$format8 == "pdf") {
                pdf(file, width = 11, height = 8.5)
              }
              matplot(t(data), type = "l",main = 'Size Profile Of Uplaoded Data',cex.lab = 1.2,col = 'grey')
              K <- r_data_global_logi()$estimates_real_logi[,2]
              x0 <- mean(data[,1])
              r <- r_data_global_logi()$estimates_real_logi[,1]
              q <- ncol(data)
              t <- 1:q
              mu_log <- K/(1+(K/x0 - 1)*exp(-r*t))
              lines(t, mu_log, col="black", lwd=3, lty=2, type = "l")
              dev.off()
            }
          }
        )



        r_data_local_logi <- reactiveVal(NULL)
        r_data_local_von <- reactiveVal(NULL)
        r_data_local_gom <- reactiveVal(NULL)

        observeEvent(input$button2,{

          if (input$model1 == 'Logistic Model'){

            P <- input$window2
            q <- ncol(data)
            t <- 1:q
            n <- nrow(data)
            r1 <- input$prob
            x0 <- mean(data[,1])
            K <- mean(data[,q])
            K1 <- K-K*0.3
            K2 <- K+K*0.3

            r_est_local_log = matrix(data = NA, nrow = q - P+1, ncol = 4)  # store local estimates
            r_cov_local_log = list()

            log_fun_local_likelihood = function(param){
              r = param[1]
              K = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_log = K/(1 + (K/x0 - 1)*exp(-r*t))

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }

              # localized estimates
              mu_log_local = mu_log[ind:(ind + P-1)]
              data_local = data[,ind:(ind + P - 1)]
              cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data_local[i,] - mu_log_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_log_local)
              }

              log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
              return(-log_lik)
            }

            show_toast('Progress',
                       text = 'Half Way Through, We Appreciate Your Patience',
                       timer = 5000,
                       position = "center",
                       type = 'success')

            for(j in 1:(q - P+1)){
              ind = j
              param_0 = c(r1, K, 3, 0.3)

              r_local_out_log = optim(param_0, log_fun_local_likelihood, method = "L-BFGS-B",
                                      hessian = TRUE, lower = c(0.001, K1, 0.0001, -0.001),
                                      upper = c(3, K2, 10, 0.999))
              r_est_local_log[j, ] =  r_local_out_log$par       # compute estimate
              r_cov_local_log[[j]] = solve(r_local_out_log$hessian)
            }
            colnames(r_est_local_log) <- c("r", "K", "sigma2", "rho")
            result <- list(estimates_local_r = r_est_local_log, cov_matrix_local_r = r_cov_local_log)
            r_data_local_logi(result)


            show_toast('Progress',
                       text = 'Done, We Appreciate Your Patience',
                       timer = 5000,
                       position = "center",
                       type = 'success')

            output$table2_r <- renderTable({
              # Create a data frame from the covariance matrix
              est_local_df <- as.data.frame(r_data_local_logi()$estimates_local_r)

              # Rename the columns
              colnames(est_local_df) <- c('r', 'K', 'Sigma', 'rho')

              # Format the numeric columns to display 5 decimal places
              est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

              head(est_local_df,8)
            }, width = "100%", align = 'c')


            output$plot_1_r <- renderPlot({
              plot(1:(length(t)-P+1), r_data_local_logi()$estimates_local_r[,1], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(r(Delta~t))[MLE]),
                   ylab = expression(hat(r)),xlab = 'Time Points')
            })

            output$plot_2_r <- renderPlot({
              plot(1:(length(t)-P+1),  r_data_local_logi()$estimates_local_r[,2], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(K(Delta~t))[MLE]),
                   ylab = expression(hat(K)),xlab = 'Time Points')
            })

          }
          else if (input$model1 == 'Von-Bertallanfy Model'){

            P <- input$window2
            q <- ncol(data)
            t <- 1:q
            n <- nrow(data)
            r <- input$von_r1
            x0 <- mean(data[,1])
            K <- mean(data[,q])
            K1 <- K-K*0.3
            K2 <- K+K*0.3

            r_est_local_v = matrix(data = NA, nrow = length(t)-P+1, ncol = 4)
            r_cov_local_v = list()

            VB_fun_local_likelihood = function(param){
              r = param[1]
              K = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_VB = K*((1 + ((x0/K)^(1/3) - 1)*exp(-r*(t/3)))^3)

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }

              # localized estimates
              mu_VB_local = mu_VB[ind:(ind + P-1)]
              data_local = data[,ind:(ind + P - 1)]
              cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data_local[i,] - mu_VB_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_VB_local)
              }

              log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
              return(-log_lik)
            }


            for(j in 1:(q - P+1)){
              ind = j
              param_0 = c(r, K, 3, 0.3)

              local_out_VB = optim(param_0, VB_fun_local_likelihood, method = "L-BFGS-B",
                                   hessian = TRUE, lower = c(0.001, K1, 0.0001, -0.001),
                                   upper = c(3, K2, 10, 0.999))
              r_est_local_v[j, ] =  local_out_VB$par       # compute estimate
              r_cov_local_v[[j]] = solve(local_out_VB$hessian)
            }
            colnames(r_est_local_v) = c("r", "K", "sigma2", "rho")

            result_vb <- list(r_estimates_local_von = r_est_local_v, r_cov_matrix_local_von = r_cov_local_v)
            r_data_local_von(result_vb)



            output$table2_r <- renderTable({
              # Create a data frame from the covariance matrix
              est_local_df <- as.data.frame(r_data_local_von()$r_estimates_local_von)

              # Rename the columns
              colnames(est_local_df) <- c('r', 'K', 'Sigma', 'rho')

              # Format the numeric columns to display 5 decimal places
              est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

              head(est_local_df,8)
            }, width = "100%", align = 'c')


            output$plot_1_r <- renderPlot({
              plot(1:(q-P+1), r_data_local_von()$r_estimates_local_von[,1], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(r(Delta~t))[MLE]),
                   ylab = expression(hat(r)),xlab = 'Time Points')
            })

            output$plot_2_r <- renderPlot({
              plot(1:(q-P+1), r_data_local_von()$r_estimates_local_von[,2], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(K(Delta~t))[MLE]),
                   ylab = expression(hat(K)),xlab = 'Time Points')
            })



          }
          else if (input$model1 == 'Gompertz Model'){


            P <- input$window2
            q <- ncol(data)
            t <- 1:q
            n <- nrow(data)
            b <- input$gom_b1
            c <- input$gom_c1
            x0 <- mean(data[,1])
            K <- mean(data[,q])
            K1 <- K-K*0.3
            K2 <- K+K*0.3


            r_est_local_gom = matrix(data = NA, nrow = q-P+1, ncol = 4)# to store the parameter estimates
            r_cov_local_gom=list() # to store the local variance-covariance matrix.

            gom_fun_local_likelihood = function(param){
              c = param[1]
              b = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_gom = x0*exp((b/c)*(1-exp(-c*t)))

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }

              # localized estimates
              mu_local_gom = mu_gom[ind:(ind + P-1)]
              data_local = data[,ind:(ind + P-1)]
              cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data_local[i,] - mu_local_gom)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_local_gom)
              }

              log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
              return(-log_lik)
            }



            for(j in 1:(q-P+1)){
              ind = j
              param_0 = c(c, b, 2, 0.3)
              local_out_gom = optim(param_0, gom_fun_local_likelihood, method = "L-BFGS-B",
                                    hessian = TRUE, lower = c(0.01, 0.01, 0.001, -0.001),
                                    upper = c(5, 10, 3, 0.999))
              r_est_local_gom[j, ] = local_out_gom$par
              r_cov_local_gom[[j]] = solve(local_out_gom$hessian)
            }
            colnames(r_est_local_gom) = c("c", "b", "sigma2", "rho")
            result <- list(r_estimates_local_gom = r_est_local_gom, r_cov_matrix_local_gom = r_cov_local_gom)
            r_data_local_gom(result)



            output$table2_r <- renderTable({
              # Create a data frame from the covariance matrix
              est_local_df <- as.data.frame(r_data_local_gom()$r_estimates_local_gom)

              # Rename the columns
              colnames(est_local_df) <-  c("c", "b", "sigma2", "rho")

              # Format the numeric columns to display 5 decimal places
              est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

              head(est_local_df,8)
            }, width = "100%", align = 'c')


            output$plot_1_r <- renderPlot({

              plot(1:(q-P+1), r_data_local_gom()$r_estimates_local_gom[,1], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(c(Delta~t))[MLE]),
                   ylab = expression(hat(c)),xlab = 'Time Points')

            })

            output$plot_2_r <- renderPlot({
              plot(1:(q-P+1), r_data_local_gom()$r_estimates_local_gom[,2], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(b(Delta~t))[MLE]),
                   ylab = expression(hat(b)),xlab = 'Time Points')

            })

          }})

        observeEvent(input$Cov_M_r,{

          if (input$model1 == 'Logistic Model'){
            showNotification('Processing...',duration=3,type = 'warning')
            q <- ncol(data)
            t <- 1:q
            P <- input$window2

            showModal(
              modalDialog(
                renderTable({
                  # Create a data frame from the covariance matrix
                  cov_local_df <- as.data.frame(r_data_local_logi()$cov_matrix_local_r[[length(t)-P+1]])

                  # Rename the columns
                  colnames(cov_local_df) <- c('r', 'K', 'Sigma', 'rho')

                  # Format the numeric columns to display 5 decimal places
                  cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

                  cov_local_df
                }, width = "100%", align = 'c',height = 10000)
                ,easyClose = TRUE,
                size = 'l',
                footer = downloadButton('cov_data_r','Covariance Data'))
            )



          }
          else if (input$model1 == 'Von-Bertallanfy Model'){
            showNotification('Processing...',duration=3,type = 'warning')

            q <- ncol(data)
            t <- 1:q
            P <- input$window2

            showModal(
              modalDialog(
                renderTable({
                  # Create a data frame from the covariance matrix
                  cov_local_df <- as.data.frame(r_data_local_von()$r_cov_matrix_local_von[[length(t)-P+1]])

                  # Rename the columns
                  colnames(cov_local_df) <-  c('r', 'K', 'Sigma', 'rho')

                  # Format the numeric columns to display 5 decimal places
                  cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

                  cov_local_df
                }, width = "100%", align = 'c',height = 10000)
                ,easyClose = TRUE,
                size = 'l',
                footer = downloadButton('cov_data_r','Covariance Data'))
            )

          } else if (input$model1 == 'Gompertz Model'){
            showNotification('Processing...',duration=3,type = 'warning')
            q <- ncol(data)
            t <- 1:q
            P <- input$window2

            showModal(
              modalDialog(
                renderTable({
                  # Create a data frame from the covariance matrix
                  cov_local_df <- as.data.frame(r_data_local_gom()$r_cov_matrix_local_gom[[length(t)-P+1]])

                  # Rename the columns
                  colnames(cov_local_df) <-  c("c", "b", "Sigma", "Rho")

                  # Format the numeric columns to display 5 decimal places
                  cov_local_df[, sapply(cov_local_df, is.numeric)] <- lapply(cov_local_df[, sapply(cov_local_df, is.numeric)], function(x) format(x, nsmall = 5))

                  cov_local_df
                }, width = "100%", align = 'c',height = 10000)
                ,easyClose = TRUE,
                size = 'l',
                footer = downloadButton('cov_data_r','Covariance Data'))
            )

          }})

        output$local_est_real <- downloadHandler(
          filename = function() {
            if (input$model1 == 'Logistic Model') {
              if (input$format9 == "csv") {
                paste("Local_Estimates_logi-", Sys.Date(), ".csv", sep = "")
              }
            } else if (input$model1 == 'Von-Bertallanfy Model'){
              if (input$format9 == "csv") {
                paste("Local_Estimates_von-", Sys.Date(), ".csv", sep = "")
              }
            } else if (input$model1 == 'Gompertz Model'){
              if (input$format9 == "csv") {
                paste("Local_Estimates_gom-", Sys.Date(), ".csv", sep = "")
              }
            }
          },
          content = function(file) {
            if (input$model1 == 'Logistic Model') {
              if (input$format9 == "csv") {
                write.csv(r_data_local_logi()$estimates_local_r, file)
              }
            } else if (input$model1 == 'Von-Bertallanfy Model'){
              if (input$format9 == "csv") {
                write.csv(r_data_local_von()$r_estimates_local_von, file)
              }
            }  else if (input$model1 == 'Gompertz Model'){
              if (input$format9 == "csv") {
                write.csv(r_data_local_gom()$r_estimates_local_gom, file)
              }
            }
          }
        )


        output$cov_data_r <- downloadHandler(
          if (input$model1 == 'Logistic Model'){
            filename = function() {
              paste("covariance_data_logi.csv")
            }} else if (input$model1 == 'Von-Bertallanfy Model'){
              filename = function() {
                paste("covariance_data_von.csv")
              }
            } else if (input$model1 == 'Gompertz Model'){
              filename = function() {
                paste("covariance_data_von.csv")
              }
            },
          if (input$model1 == 'Logistic Model'){
            content = function(file) {
              write.csv( r_data_local_logi()$cov_matrix_local_r, file)}

          }   else if (input$model1 == 'Von-Bertallanfy Model'){
            content = function(file) {
              write.csv(r_data_local_von()$r_cov_matrix_local_von, file)}
          }  else if (input$model1 == 'Gompertz Model'){
            content = function(file) {
              write.csv(r_data_local_gom()$r_cov_matrix_local_gom, file)}
          })


        observeEvent(input$zoom10,{
          if (input$model1 == 'Logistic Model'){
            showNotification('Processing...',duration=3,type = 'warning')

            showModal(
              modalDialog(
                paste('Local Estimation Table (Complete)'),
                br(),
                renderTable({
                  # Create a data frame from the covariance matrix
                  est_local_df <- as.data.frame(r_data_local_logi()$estimates_local_r)

                  # Rename the columns
                  colnames(est_local_df) <- c('r', 'K', 'Sigma', 'rho')

                  # Format the numeric columns to display 5 decimal places
                  est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

                  est_local_df
                }, width = "100%", align = 'c',height = 10000)
              )
            )

          }  else if (input$model1 == 'Von-Bertallanfy Model'){
            showNotification('Processing...',duration=3,type = 'warning')

            showModal(
              modalDialog(
                paste('Local Estimation Table (Complete)'),
                br(),
                renderTable({
                  # Create a data frame from the covariance matrix
                  est_local_df <- as.data.frame(r_data_local_von()$r_estimates_local_von)

                  # Rename the columns
                  colnames(est_local_df) <- c('r','K', 'Sigma', 'rho')

                  # Format the numeric columns to display 5 decimal places
                  est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

                  est_local_df
                }, width = "100%", align = 'c',height = 10000)
              )
            )
          } else if (input$model == 'Gompertz Model'){
            showNotification('Processing...',type = 'warning',duration=3)

            showModal(
              modalDialog(
                paste('Local Estimation Table (Complete)'),
                br(),
                renderTable({
                  # Create a data frame from the covariance matrix
                  est_local_df <- as.data.frame(r_data_local_gom()$r_estimates_local_gom)

                  # Rename the columns
                  colnames(est_local_df) <- c("c", "b", "Sigma", "Rho")

                  # Format the numeric columns to display 5 decimal places
                  est_local_df[, sapply(est_local_df, is.numeric)] <- lapply(est_local_df[, sapply(est_local_df, is.numeric)], function(x) format(x, nsmall = 5))

                  est_local_df
                }, width = "100%", align = 'c',height = 10000)
              )
            )
          }})


        #zooming of estimates plot (local)

        observeEvent(input$zoom11,{

          if (input$model1 == 'Logistic Model'){
            P <- input$window2
            q <- ncol(data)
            t <- 1:q

            showModal(modalDialog(

              renderPlot({
                plot(1:(length(t)-P+1), r_data_local_logi()$estimates_local_r[,1], type = "b", col = "red", lwd = 2,
                     main = expression(widehat(r(Delta~t))[MLE]),
                     ylab = expression(hat(r)),xlab = 'Time Points')
              }, height = 500),

              renderPlot({
                plot(1:(length(t)-P+1),  r_data_local_logi()$estimates_local_r[,2], type = "b", col = "red", lwd = 2,
                     main = expression(widehat(K(Delta~t))[MLE]),
                     ylab = expression(hat(K)),xlab = 'Time Points')
              },height=500),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          }   else if (input$model1 == 'Von-Bertallanfy Model'){

            P <- input$window2
            q <- ncol(data)
            t <- 1:q

            showModal(modalDialog(


              renderPlot({
                plot(1:(q-P+1), r_data_local_von()$r_estimates_local_von[,1], type = "b", col = "red", lwd = 2,
                     main = expression(widehat(r(Delta~t))[MLE]),
                     ylab = expression(hat(r)),xlab = 'Time Points')

              },height = 500),

              renderPlot({
                plot(1:(q-P+1), r_data_local_von()$r_estimates_local_von[,2], type = "b", col = "red", lwd = 2,
                     main = expression(widehat(K(Delta~t))[MLE]),
                     ylab = expression(hat(K)),xlab = 'Time Points')

              },height = 500),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          } else if (input$model1 == 'Gompertz Model'){

            P <- input$window2
            q <- ncol(data)


            showModal(modalDialog(


              renderPlot({
                plot(1:(q-P+1), r_data_local_gom()$r_estimates_local_gom[,1], type = "b", col = "red", lwd = 2,
                     main = expression(widehat(c(Delta~t))[MLE]),
                     ylab = expression(hat(c)),xlab = 'Time Points')


              },height = 500),

              renderPlot({
                plot(1:(q-P+1), r_data_local_gom()$r_estimates_local_gom[,2], type = "b", col = "red", lwd = 2,
                     main = expression(widehat(b(Delta~t))[MLE]),
                     ylab = expression(hat(b)),xlab = 'Time Points')


              },height = 500),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          }})


        output$l_plot_d_real <- downloadHandler(
          filename = function() {
            if (input$model1 == 'Logistic Model') {
              if (input$format10 == "png") {
                paste("data_plot_local_logi-", Sys.Date(), ".png", sep = "")
              } else if (input$format10 == "pdf") {
                paste("data_plot_local_logi-", Sys.Date(), ".pdf", sep = "")
              }
            }  else if (input$model1 == 'Von-Bertallanfy Model'){
              if (input$format10 == "png") {
                paste("data_plot_local_von-", Sys.Date(), ".png", sep = "")
              } else if (input$format10 == "pdf") {
                paste("data_plot_local_von-", Sys.Date(), ".pdf", sep = "")
              }
            } else if (input$model1 == 'Gompertz Model'){
              if (input$format10 == "png") {
                paste("data_plot_local_gom-", Sys.Date(), ".png", sep = "")
              } else if (input$format10 == "pdf") {
                paste("data_plot_local_gom-", Sys.Date(), ".pdf", sep = "")
              }
            }
          },
          content = function(file) {
            if (input$model1 == 'Logistic Model') {
              # Plot for Logistic Model
              P <- input$window2
              q <- ncol(data)
              t <- 1:q

              if (input$format10 == "png") {
                png(file, width = 1100, height = 500)
              } else if (input$format10 == "pdf") {
                pdf(file,  width = 11, height = 8.5)
              }

              par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
              plot(1:(length(t)-P+1), r_data_local_logi()$estimates_local_r[,1], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(r(Delta~t))[MLE]),
                   ylab = expression(hat(r)),xlab = 'Time Points')

              plot(1:(length(t)-P+1),  r_data_local_logi()$estimates_local_r[,2], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(K(Delta~t))[MLE]),
                   ylab = expression(hat(K)),xlab = 'Time Points')

              dev.off()
            } else if (input$model1 == 'Von-Bertallanfy Model'){

              P <- input$window2
              q <- ncol(data)
              t <- 1:q

              if (input$format10 == "png") {
                png(file, width = 1100, height = 500)
              } else if (input$format10 == "pdf") {
                pdf(file, width = 11, height = 8.5)
              }

              par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
              plot(1:(q-P+1), r_data_local_von()$r_estimates_local_von[,1], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(r(Delta~t))[MLE]),
                   ylab = expression(hat(r)),xlab = 'Time Points')

              plot(1:(q-P+1), r_data_local_von()$r_estimates_local_von[,2], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(K(Delta~t))[MLE]),
                   ylab = expression(hat(K)),xlab = 'Time Points')


              dev.off()
            } else if (input$model1 == 'Gompertz Model'){
              P <- input$window2
              q <- ncol(data)
              t <- 1:q

              if (input$format10 == "png") {
                png(file, width = 1100, height = 500)
              } else if (input$format10 == "pdf") {
                pdf(file, width = 11, height = 8.5)
              }

              par(mfrow = c(2, 1))  # Set up a 1x2 grid for two plots side by side
              plot(1:(q-P+1), r_data_local_gom()$r_estimates_local_gom[,1], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(c(Delta~t))[MLE]),
                   ylab = expression(hat(c)),xlab = 'Time Points')

              plot(1:(q-P+1), r_data_local_gom()$r_estimates_local_gom[,2], type = "b", col = "red", lwd = 2,
                   main = expression(widehat(b(Delta~t))[MLE]),
                   ylab = expression(hat(b)),xlab = 'Time Points')



              dev.off()
            }
          }
        )


        #-------------------------------- P value--------------------
        p_val_log <- reactiveVal(NULL)
        p_val_von <- reactiveVal(NULL)
        p_val_gom <- reactiveVal(NULL)

        observeEvent(input$button2,{

          if (input$model1 == 'Logistic Model'){
            q <- ncol(data)
            P <- input$window2


            RGR_hat_log <- numeric(length = q-P+1)
            for(i in 1:(q-P+1)){
              RGR_hat_log[i] <- r_data_local_logi()$estimates_local_r[i,1]*(1- mean(data[,i])/r_data_local_logi()$estimates_local_r[i,2])
            }

            #for(i in 1:length(r_data_local_logi()$cov_matrix_local_r)){
            # colnames(r_data_local_logi()$cov_matrix_local_r[[i]]) <- c("r", "K", "sigma2", "rho")
            #rownames(r_data_local_logi()$cov_matrix_local_r[[i]]) = c("r", "K", "sigma2", "rho")
            #}

            var_RGR_hat_log <- numeric(length = q-P+1)
            for(i in 1:(q-2)){
              A = r_data_local_logi()$cov_matrix_local_r[[i]][1,1]*(1-mean(data[,i])/r_data_local_logi()$estimates_local_r[i,2])^2
              B = r_data_local_logi()$cov_matrix_local_r[[i]][2,2]*(r_data_local_logi()$estimates_local_r[i,1]*mean(data[,i])/r_data_local_logi()$estimates_local_r[i,2]^2)^2
              C = 2*(1-mean(data[,i])/r_data_local_logi()$estimates_local_r[i,2])*(r_data_local_logi()$estimates_local_r[i,1]*mean(data[,i])/r_data_local_logi()$estimates_local_r[i,2]^2)*r_data_local_logi()$cov_matrix_local_r[[i]][1,2]
              var_RGR_hat_log[i] = A+B+C
            }

            V_log <- numeric(length = q-2)
            for(i in 1:(q-2)){
              V_log[i] = log(1/RGR_hat_log[i] - 1/r_data_local_logi()$estimates_local_r[i,1])
            }

            var_V_log <- numeric(length = q-2)
            for(i in 1:(q-2)){
              A = (-r_data_local_logi()$estimates_local_r[i,1]/(RGR_hat_log[i]*(r_data_local_logi()$estimates_local_r[i,1]-RGR_hat_log[i])))^2*var_RGR_hat_log[i]
              B = (RGR_hat_log[i]/(r_data_local_logi()$estimates_local_r[i,1]*(r_data_local_logi()$estimates_local_r[i,1]-RGR_hat_log[i])))^2*r_data_local_logi()$cov_matrix_local_r[[i]][1,1]
              C = 2*(-1/(r_data_local_logi()$estimates_local_r[i,1]-RGR_hat_log[i])^2)*(r_data_local_logi()$cov_matrix_local_r[[i]][1,1]*(1-mean(data[,i])/r_data_local_logi()$estimates_local_r[i,2]) + (r_data_local_logi()$estimates_local_r[i,1]*mean(data[,i])/(r_data_local_logi()$estimates_local_r[i,2])^2)*r_data_local_logi()$cov_matrix_local_r[[i]][1,2])
              var_V_log[i] = A + B + C
            }

            test_stat_log <- numeric(length = q-P-1)
            var_test_stat_log <- numeric(length = q-P-1)
            for(i in 1:(q-P-1)){
              var_test_stat_log[i] = var_V_log[i+2] + 4*var_V_log[i+1] + var_V_log[i]
              test_stat_log[i] = (V_log[i+2]-2*V_log[i+1] + V_log[i])
            }

            p_value_local_log_l <- numeric(q-P-1)
            for(i in 1:(q-P-1)){
              p_value_local_log_l[i] = 1- pnorm(abs(test_stat_log[i]/sqrt(var_test_stat_log[i])))
            }

            p_val_log(p_value_local_log_l)
            output$p_plot <- renderPlot(
              {

                plot(1:(q-P-1), p_val_log(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                points(1:(q-P-1), p_val_log(), col = "red", pch = 19)
                segments(1:(q-P-1), 0, 1:(q-P-1), p_val_log(), col = "black", lwd = 2)
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

              }
            )

          } else if (input$model1 == 'Von-Bertallanfy Model'){
            q <- ncol(data)
            P <- input$window2

            RGR_hat_VB = numeric(length = q-P+1)
            for(i in 1:(q-P+1)){
              RGR_hat_VB[i] = r_data_local_von()$r_estimates_local_von[i,1]*((r_data_local_von()$r_estimates_local_von[i,2]/mean(data[,i]))^(1/3))*(1- (mean(data[,i])/r_data_local_von()$r_estimates_local_von[i,2])^(1/3))
            }


            var_RGR_hat_VB <- numeric(length = q-P+1)
            for(i in 1:(q-2)){
              A = r_data_local_von()$r_cov_matrix_local_von[[i]][1,1]*((r_data_local_von()$r_estimates_local_von[i,2]/mean(data[,i]))^(1/3)-1)^2
              B = r_data_local_von()$r_cov_matrix_local_von[[i]][2,2]*(r_data_local_von()$r_estimates_local_von[i,1]/(3*(mean(data[,i])*r_data_local_von()$r_estimates_local_von[i,2]^2)^(1/3)))^2
              C = 2*((r_data_local_von()$r_estimates_local_von[i,2]/mean(data[,i]))^(1/3)-1)*(r_data_local_von()$r_estimates_local_von[i,1]/(3*(mean(data[,i])*r_data_local_von()$r_estimates_local_von[i,2]^2)^(1/3)))*r_data_local_von()$r_cov_matrix_local_von[[i]][1,2]
              var_RGR_hat_VB[i] = A+B+C
            }

            # Computation of V
            V_VB <- numeric(length = q-2)
            for(i in 1:(q-2)){
              V_VB[i] = log(1/RGR_hat_VB[i] + 1/r_data_local_von()$r_estimates_local_von[i,1])
            }
            print(V_VB)

            # Computation of Variance of V
            var_V_VB <- numeric(length = q-2)
            for(i in 1:(q-2)){
              A = (-r_data_local_von()$r_estimates_local_von[i,1]/(RGR_hat_VB[i]*(r_data_local_von()$r_estimates_local_von[i,1]+RGR_hat_VB[i])))^2*var_RGR_hat_VB[i]
              B = (-RGR_hat_VB[i]/(r_data_local_von()$r_estimates_local_von[i,1]*(r_data_local_von()$r_estimates_local_von[i,1]+RGR_hat_VB[i])))^2*r_data_local_von()$r_cov_matrix_local_von[[i]][1,1]
              C = 2*(1/(r_data_local_von()$r_estimates_local_von[i,1]+RGR_hat_VB[i])^2)*(r_data_local_von()$r_cov_matrix_local_von[[i]][1,1]*(1-mean(data[,i])/r_data_local_von()$r_estimates_local_von[i,2]) + (r_data_local_von()$r_estimates_local_von[i,1]*mean(data[,i])/(r_data_local_von()$r_estimates_local_von[i,2])^2)*r_data_local_von()$r_cov_matrix_local_von[[i]][1,2])
              var_V_VB[i] = A + B + C
            }
            var_V_VB

            # Computation of the test statistic
            test_stat_VB <- numeric(length = q-P-1)
            var_test_stat_VB = numeric(length = q-P-1)
            for(i in 1:(q-P-1)){
              var_test_stat_VB[i] = var_V_VB[i+2] + 4*var_V_VB[i+1] + var_V_VB[i]
              test_stat_VB[i] = (V_VB[i+2]-2*V_VB[i+1] + V_VB[i])
            }

            abs(test_stat_VB/sqrt(var_test_stat_VB))

            # Computation of the p-value
            p_value_local_VB <- numeric(q-P-1)
            for(i in 1:(q-P-1)){
              p_value_local_VB[i] = 1- pnorm(test_stat_VB[i]/sqrt(var_test_stat_VB[i]))
            }

            p_val_von(p_value_local_VB)
            output$p_plot <- renderPlot(
              {

                plot(1:(q-P-1), p_val_von(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot')
                points(1:(q-P-1), p_val_von(), col = "red", pch = 19)
                segments(1:(q-P-1), 0, 1:(q-P-1), p_val_von(), col = "black", lwd = 2)
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)}
            )

          } else if (input$model1 == 'Gompertz Model'){
            q <- ncol(data)
            P <- input$window2

            RGR_hat_gom <- numeric(length = q-P+1)
            for(i in 1:(q-P+1)){
              RGR_hat_gom[i] = r_data_local_gom()$r_estimates_local_gom[i,2]*exp(-r_data_local_gom()$r_estimates_local_gom[i,1]*t[i])
            }



            # Computation of variance of Modified RGR using Delta method
            var_RGR_hat_gom <- numeric(length = q-P+1)
            for(i in 1:(q-2)){
              A = r_data_local_gom()$r_cov_matrix_local_gom[[i]][1,1]*(-r_data_local_gom()$r_estimates_local_gom[i,2]*t[i]*exp(-r_data_local_gom()$r_estimates_local_gom[i,1]*t[i]))^2
              B = r_data_local_gom()$r_cov_matrix_local_gom[[i]][2,2]*(exp(-r_data_local_gom()$r_estimates_local_gom[i,1]*t[i]))^2
              C = 2*(-r_data_local_gom()$r_estimates_local_gom[i,2]*t[i]*exp(-r_data_local_gom()$r_estimates_local_gom[i,1]*t[i]))*(exp(-r_data_local_gom()$r_estimates_local_gom[i,1]*t[i]))*r_data_local_gom()$r_cov_matrix_local_gom[[i]][1,2]
              var_RGR_hat_gom[i] = A+B+C
            }

            # Computation of V
            V_gom = numeric(length = q-2)
            for(i in 1:(q-2)){
              V_gom[i] = log(1/RGR_hat_gom[i])
            }


            # Computation of Variance of V
            var_V_gom <- numeric(length = q-2)
            for(i in 1:(q-2)){
              A = r_data_local_gom()$r_cov_matrix_local_gom[[i]][1,1]*(-t[i])^2
              B = r_data_local_gom()$r_cov_matrix_local_gom[[i]][2,2]*(1/r_data_local_gom()$r_estimates_local_gom[i,2])^2
              C = 2*(-t[i]/r_data_local_gom()$r_estimates_local_gom[i,2])*r_data_local_gom()$r_cov_matrix_local_gom[[i]][1,2]
              var_V_gom[i] = A + B + C
            }


            # Computation of the test statistic
            test_stat_gom <- numeric(length = q-P-1)
            var_test_stat_gom <- numeric(length = q-P-1)
            for(i in 1:(q-P-1)){
              var_test_stat_gom[i] = var_V_gom[i+2] + 4*var_V_gom[i+1] + var_V_gom[i]
              test_stat_gom[i] = (V_gom[i+2]-2*V_gom[i+1] + V_gom[i])
            }



            # Computation of the p-value
            p_value_local_gom <- numeric(q-P-1)
            for(i in 1:(q-P-1)){
              p_value_local_gom[i] = 1- pnorm(abs(test_stat_gom[i]/sqrt(var_test_stat_gom[i])))
            }
            p_val_gom(p_value_local_gom)

            output$p_plot <- renderPlot(
              {

                plot(1:(q-P-1), p_val_gom(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                points(1:(q-P-1), p_val_gom(), col = "red", pch = 19)
                segments(1:(q-P-1), 0, 1:(q-P-1), p_val_gom(), col = "black", lwd = 2)
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)
              }
            )

          }


          observeEvent(input$zoom12,{
            if (input$model1 == 'Logistic Model'){
              P <- input$window2
              q <- ncol(data)
              t <- 1:q

              showModal(modalDialog(

                renderPlot({

                  plot(1:(q-P-1), p_val_log(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                  points(1:(q-P-1), p_val_log(), col = "red", pch = 19)
                  segments(1:(q-P-1), 0, 1:(q-P-1), p_val_log(), col = "black", lwd = 2)
                  abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)
                }, height = 500),
                easyClose = TRUE,
                size = "l",
                footer = NULL
              ))
            }   else if (input$model1 == 'Von-Bertallanfy Model'){

              P <- input$window2
              q <- ncol(data)
              t <- 1:q

              showModal(modalDialog(


                renderPlot({
                  plot(1:(q-P-1), p_val_von(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                  points(1:(q-P-1), p_val_von(), col = "red", pch = 19)
                  segments(1:(q-P-1), 0, 1:(q-P-1), p_val_von(), col = "black", lwd = 2)
                  abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

                },height = 500),
                easyClose = TRUE,
                size = "l",
                footer = NULL
              ))
            } else if (input$model1 == 'Gompertz Model'){

              P <- input$window2
              q <- ncol(data)
              t <- 1:q

              showModal(modalDialog(


                renderPlot({

                  plot(1:(q-P-1), p_val_gom(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                  points(1:(q-P-1), p_val_gom(), col = "red", pch = 19)
                  segments(1:(q-P-1), 0, 1:(q-P-1), p_val_gom(), col = "black", lwd = 2)
                  abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)


                },height = 500),
                easyClose = TRUE,
                size = "l",
                footer = NULL
              ))
            }
          })


          output$p_plot_r <- downloadHandler(
            filename = function() {
              if (input$model1 == 'Logistic Model') {
                if (input$format11 == "png") {
                  paste("p_val_logi-", Sys.Date(), ".png", sep = "")
                } else if (input$format11 == "pdf") {
                  paste("p_val_logi-", Sys.Date(), ".pdf", sep = "")
                }
              } else if (input$model1 == 'Von-Bertallanfy Model'){
                if (input$format11 == "png") {
                  paste("p_val_von-", Sys.Date(), ".png", sep = "")
                } else if (input$format11 == "pdf") {
                  paste("p_val_von-", Sys.Date(), ".pdf", sep = "")
                }
              } else if (input$model == 'Gompertz Model'){
                if (input$format11 == "png") {
                  paste("p_val_gom-", Sys.Date(), ".png", sep = "")
                } else if (input$format11 == "pdf") {
                  paste("p_val_gom-", Sys.Date(), ".pdf", sep = "")
                }
              }
            },
            content = function(file) {
              if (input$model1 == 'Logistic Model') {
                # Plot for Logistic Model
                q <- ncol(data)
                P <- input$window2
                t <- 1:q

                if (input$format11 == "png") {
                  png(file, width = 1100, height = 500)
                } else if (input$format11 == "pdf") {
                  pdf(file,  width = 11, height = 8.5)
                }

                plot(1:(q-P-1), p_val_log(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                points(1:(q-P-1), p_val_log(), col = "red", pch = 19)
                segments(1:(q-P-1), 0, 1:(q-P-1), p_val_log(), col = "black", lwd = 2)
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)




                dev.off()
              } else if (input$model1 == 'Von-Bertallanfy Model'){
                q <- ncol(data)
                P <- input$window2
                t <- 1:q
                if (input$format11 == "png") {
                  png(file, width = 1100, height = 500)
                } else if (input$format3 == "pdf") {
                  pdf(file, width = 11, height = 8.5)
                }


                plot(1:(q-P-1), p_val_von(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                points(1:(q-P-1), p_val_von(), col = "red", pch = 19)
                segments(1:(q-P-1), 0, 1:(q-P-1), p_val_von(), col = "black", lwd = 2)
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

                dev.off()
              } else if (input$model1 == 'Gompertz Model'){
                q <- ncol(data)
                P <- input$window2
                t <- 1:q
                if (input$format11 == "png") {
                  png(file, width = 1100, height = 500)
                } else if (input$format11 == "pdf") {
                  pdf(file, width = 11, height = 8.5)
                }


                plot(1:(q-P-1), p_val_gom(), type = "n", ylim = c(0, 1), ylab = "p-value", lwd = 2, pch = 19, main = 'P Value Plot',xlab = 'Time Points')
                points(1:(q-P-1), p_val_gom(), col = "red", pch = 19)
                segments(1:(q-P-1), 0, 1:(q-P-1), p_val_gom(), col = "black", lwd = 2)
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)



                dev.off()
              }
            }
          )
        })

        # comparison plot

        observeEvent(input$button3,{




          show_toast('Important',
                     text = 'Dont click this button twice, Computation May take little time',
                     timer = 11000,
                     position = "center",
                     type = 'success')

          P <- input$window2
          q <- ncol(data)
          t <- 1:q
          n <- nrow(data)
          r1 <- input$prob
          b <- input$gom_b1
          c <- input$gom_c1
          r <- input$von_r1
          x0 <- mean(data[,1])
          x01 <- x0-x0*0.1
          x02 <- x0+x0*0.1
          K <- mean(data[,q])
          K1 <- K-K*0.3
          K2 <- K+K*0.3
          B_ <- 20


          ## -------------- ,Log Param likelihood below---------------

          r_est_local_log = matrix(data = NA, nrow = q - P+1, ncol = 4)  # store local estimates
          r_cov_local_log = list()

          log_fun_local_likelihood = function(param){
            r = param[1]
            K = param[2]
            sigma2  = param[3]
            rho = param[4]
            mu_log = K/(1 + (K/x0 - 1)*exp(-r*t))

            cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
            for (i in 1:q) {
              for (j in 1:q) {
                cov_mat[i,j] = sigma2*rho^(abs(i-j))
              }
            }

            # localized estimates
            mu_log_local = mu_log[ind:(ind + P-1)]
            data_local = data[,ind:(ind + P - 1)]
            cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

            exponent = 0
            for(i in 1:n){
              exponent = exponent + t(data_local[i,] - mu_log_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_log_local)
            }

            log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
            return(-log_lik)
          }


          for(j in 1:(q - P+1)){
            ind = j
            param_0 = c(r, K, 3, 0.3)

            r_local_out_log = optim(param_0, log_fun_local_likelihood, method = "L-BFGS-B",
                                    hessian = TRUE, lower = c(0.001, K1, 0.0001, -0.001),
                                    upper = c(3, K2, 10, 0.999))
            r_est_local_log[j, ] =  r_local_out_log$par       # compute estimate
            r_cov_local_log[[j]] = solve(r_local_out_log$hessian)
          }
          colnames(r_est_local_log) <- c("r", "K", "sigma2", "rho")

          ## -------------- ,Log Param likelihood above---------------


          ## -------------- ,Log Non Param likelihood below---------------
          est_local_log_list_np = list()
          x0_mean_vals_log_np = numeric(length = B_)
          data_np_boot_list = list()

          for(l in 1:B_){
            set.seed(l)
            print(l)
            boot = sample(1:n, size = n, replace = TRUE)
            data_np_boot  = data[boot, ]
            data_np_boot_list[[l]] = data_np_boot

            est_local_log_np_boot = matrix(data = NA, nrow = q-P+1, ncol = 4)
            x0_mean_vals_log_np = mean(data_np_boot[,1])

            log_fun_local_likelihood_np_boot = function(param){
              r = param[1]
              K = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_log = K/(1 + (K/x0 - 1)*exp(-r*t))

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }

              # localized estimates
              mu_log_local = mu_log[ind:(ind + P-1)]
              data_local = data_np_boot[,ind:(ind + P - 1)]
              cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data_local[i,] - mu_log_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_log_local)
              }

              log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
              return(-log_lik)
            }

            for(j in 1:(q - P+1)){
              ind = j
              param_0 = c(r, K, 3, 0.3)

              local_out_log_np_boot = optim(param_0, log_fun_local_likelihood_np_boot, method = "L-BFGS-B",
                                            hessian = TRUE, lower = c(0.001, K1, 0.5, -0.001),
                                            upper = c(3, K2, 10, 0.9))
              est_local_log_np_boot[j, ] =  local_out_log_np_boot$par       # compute estimate
            }

            colnames(est_local_log_np_boot) = c("r", "K", "sigma2", "rho")
            est_local_log_list_np[[l]] = est_local_log_np_boot
          }



          ## -------------- ,Log Non Param likelihood above---------------




          show_toast('Processing...',
                     text = 'Half way through...',
                     timer = 8000,
                     position = "center",
                     type = 'success')

          ## -------------- VB Param likelihood below---------------

          r_est_local_v = matrix(data = NA, nrow = length(t)-P+1, ncol = 4)
          r_cov_local_v = list()

          VB_fun_local_likelihood = function(param){
            r = param[1]
            K = param[2]
            sigma2  = param[3]
            rho = param[4]
            mu_VB = K*((1 + ((x0/K)^(1/3) - 1)*exp(-r*(t/3)))^3)

            cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
            for (i in 1:q) {
              for (j in 1:q) {
                cov_mat[i,j] = sigma2*rho^(abs(i-j))
              }
            }

            # localized estimates
            mu_VB_local = mu_VB[ind:(ind + P-1)]
            data_local = data[,ind:(ind + P - 1)]
            cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

            exponent = 0
            for(i in 1:n){
              exponent = exponent + t(data_local[i,] - mu_VB_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_VB_local)
            }

            log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
            return(-log_lik)
          }


          for(j in 1:(q - P+1)){
            ind = j
            param_0 = c(r, K, 3, 0.3)

            local_out_VB = optim(param_0, VB_fun_local_likelihood, method = "L-BFGS-B",
                                 hessian = TRUE, lower = c(0.001, K1, 0.0001, -0.001),
                                 upper = c(3, K2, 10, 0.999))
            r_est_local_v[j, ] =  local_out_VB$par       # compute estimate
            r_cov_local_v[[j]] = solve(local_out_VB$hessian)
          }
          colnames(r_est_local_v) = c("r", "K", "sigma2", "rho")

          ## -------------- VB Param likelihood above---------------



          ## -------------- VB Non Param likelihood below---------------

          est_local_VB_list_np = list()
          x0_mean_vals_VB_np = numeric(length = B_)
          data_np_boot_VB_list = list()
          for(l in 1:B_){
            set.seed(l)
            print(l)
            boot = sample(1:n, size = n, replace = TRUE)
            data_np_boot  = data[boot, ]
            data_np_boot_VB_list[[l]] = data_np_boot

            est_local_VB_np_boot = matrix(data = NA, nrow = q-P+1, ncol = 4)
            x0_mean_vals_VB_np = mean(data_np_boot[,1])

            VB_fun_local_likelihood_np_boot = function(param){
              r = param[1]
              K = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_VB = K*((1 + ((x0/K)^(1/3) - 1)*exp(-r*(t/3)))^3)

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }

              # localized estimates
              mu_VB_local = mu_VB[ind:(ind + P-1)]
              data_local = data_np_boot[,ind:(ind + P - 1)]
              cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data_local[i,] - mu_VB_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_VB_local)
              }

              log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
              return(-log_lik)
            }
            for(j in 1:(q - P+1)){
              ind = j
              param_0 = c(r, K, 3, 0.3)

              local_out_VB_np_boot = optim(param_0, VB_fun_local_likelihood_np_boot, method = "L-BFGS-B",
                                           hessian = TRUE, lower = c(0.001, K1, 0.5, -0.001),
                                           upper = c(3, K2, 10, 0.9))
              est_local_VB_np_boot[j, ] =  local_out_VB_np_boot$par       # compute estimate
            }
            colnames(est_local_VB_np_boot) = c("r", "K", "sigma2", "rho")
            est_local_VB_list_np[[l]] = est_local_VB_np_boot
          }


          ## -------------- VB NOn Param likelihood above---------------



          # --------------------- Gom Param Likelihood below ----------------

          r_est_local_gom <- matrix(data = NA, nrow = q-P+1, ncol = 4)# to store the parameter estimates
          r_cov_local_gom <- list() # to store the local variance-covariance matrix.

          gom_fun_local_likelihood = function(param){
            c = param[1]
            b = param[2]
            sigma2  = param[3]
            rho = param[4]
            mu_gom = x0*exp((b/c)*(1-exp(-c*t)))

            cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
            for (i in 1:q) {
              for (j in 1:q) {
                cov_mat[i,j] = sigma2*rho^(abs(i-j))
              }
            }

            # localized estimates
            mu_local_gom = mu_gom[ind:(ind + P-1)]
            data_local = data[,ind:(ind + P-1)]
            cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

            exponent = 0
            for(i in 1:n){
              exponent = exponent + t(data_local[i,] - mu_local_gom)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_local_gom)
            }

            log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
            return(-log_lik)
          }



          for(j in 1:(q-P+1)){
            ind = j
            param_0 = c(c, b, 2, 0.3)
            local_out_gom = optim(param_0, gom_fun_local_likelihood, method = "L-BFGS-B",
                                  hessian = TRUE, lower = c(0.01, 0.01, 0.001, -0.001),
                                  upper = c(5, 10, 3, 0.999))
            r_est_local_gom[j, ] = local_out_gom$par
            r_cov_local_gom[[j]] = solve(local_out_gom$hessian)
          }
          colnames(r_est_local_gom) <- c("c", "b", "sigma2", "rho")

          # --------------------- Gom Param Likelihood above ----------------


          # --------------------- Gom Non Param Likelihood below ----------------

          est_local_gom_list_np = list()
          x0_mean_vals_gom_np = numeric(length = B_)
          data_np_boot_list = list()
          for(l in 1:B_){
            set.seed(l)
            print(l)
            boot = sample(1:n, size = n, replace = TRUE)
            data_np_boot  = data[boot, ]
            data_np_boot_list[[l]] = data_np_boot

            est_local_gom_np_boot = matrix(data = NA, nrow = q-P+1, ncol = 4)
            x0_mean_vals_gom_np = mean(data_np_boot[,1])

            gom_fun_local_likelihood_np_boot = function(param){
              b = param[1]
              c = param[2]
              sigma2  = param[3]
              rho = param[4]
              mu_gom = x0*exp((b/c)*(1-exp(-c*t)))

              cov_mat = matrix(data = NA, nrow = length(t), ncol = length(t))
              for (i in 1:q) {
                for (j in 1:q) {
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }

              # localized estimates
              mu_gom_local = mu_gom[ind:(ind + P-1)]
              data_local = data_np_boot[,ind:(ind + P - 1)]
              cov_mat_local = cov_mat[ind:(ind+P-1), ind:(ind+P-1)]

              exponent = 0
              for(i in 1:n){
                exponent = exponent + t(data_local[i,] - mu_gom_local)%*%(solve(cov_mat_local))%*%(data_local[i,] - mu_gom_local)
              }

              log_lik = - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
              return(-log_lik)
            }
            for(j in 1:(q - P+1)){
              ind = j
              param_0 = c(b, c, 3, 0.3)

              local_out_gom_np_boot = optim(param_0, gom_fun_local_likelihood_np_boot, method = "L-BFGS-B",
                                            hessian = TRUE, lower = c(0.01, 0.01, 0.5, -0.001),
                                            upper = c(10, 10, 10, 0.9))
              est_local_gom_np_boot[j, ] =  local_out_gom_np_boot$par       # compute estimate
            }
            colnames(est_local_gom_np_boot) = c("b", "c", "sigma2", "rho")
            est_local_gom_list_np[[l]] = est_local_gom_np_boot
          }


          # --------------------- Gom Non Param Likelihood above ----------------


          show_toast('Only a bit more',
                     text = 'A little more..',
                     timer = 8000,
                     position = "center",
                     type = 'success')


          # ----- RGR log param ----
          RGR_hat_log <- numeric(length = q-P+1)
          for(i in 1:(q-P+1)){
            RGR_hat_log[i] <- r_est_local_log[i,1]*(1- mean(data[,i])/r_est_local_log[i,2])
          }


          V_log <- numeric(length = q-2)
          for(i in 1:(q-2)){
            V_log[i] <- log(1/RGR_hat_log[i] - 1/r_est_local_log[i,1])
          }


          test_stat_log <- numeric(length = q-P-1)
          #var_test_stat_log <- numeric(length = q-P-1)
          for(i in 1:(q-P-1)){
            #var_test_stat_log[i] <- var_V_log[i+2] + 4*var_V_log[i+1] + var_V_log[i]
            test_stat_log[i] <- (V_log[i+2]-2*V_log[i+1] + V_log[i])

          }



          #--------------------------------------------
          # BOOTSTRAPPING LOGISTIC

          mu_log = K/(1 + (K/x0 - 1)*exp(-r1*t))

          cov_mat_logi = matrix(data = NA, nrow = length(t), ncol = length(t))
          for (i in 1:q) {
            for (j in 1:q) {
              cov_mat_logi[i,j] = sigma2*rho^(abs(i-j))
            }
          }



          B <- 20

          est_local_log_list_p <- list()
          x0_mean_vals_log_p <-  numeric(length = B)
          data_boot_list_p <- list()

          for (l in 1:B){
            data_boot <- mvtnorm::rmvnorm(n=1,mean=mu_log,sigma=cov_mat_logi)
            data_boot_list_p[[l]] <- data_boot
            est_local_log_boot <- matrix(data = NA, nrow = q-P+1, ncol = 4)
            x0_mean_vals_log_p <- mean(data_boot[,1])

            log_fun_local_likelihood_boot <- function(param){
              r = param[1]
              K = param[2]
              sigma2 = param[3]
              rhp = param[4]
              mu_log <- K/(1+(K/x0-1)*exp(-r*t))

              cov_mat <- matrix(data=NA , nrow = length(t),ncol = length(t))
              for (i in 1:q){
                for (j in 1:q){
                  cov_mat[i,j] = sigma2*rho^(abs(i-j))
                }
              }

              mu_log_local <- mu_log[ind:(ind+P-1)]
              data_local <- data_boot[,ind:(ind+P-1)]
              cov_mat_local <- cov_mat[ind:(ind+P-1),ind:(ind+P-1)]

              exponent <- 0
              exponent <- exponent + t(data_local-mu_log_local)%*%(solve(cov_mat_local))%*%(data_local - mu_log_local)
              log_lik <- - (n*length(ind:(ind+P-1))/2)*log(2*pi) - (n/2)*log(det(cov_mat_local)) - (1/2)*exponent
              return(-log_lik)

            }

            for (j in 1:(q-P+1)){
              ind <- j
              param_0 <- c(r,K,3,0.3)
              local_out_log_boot = optim(param_0, log_fun_local_likelihood_boot, method = "L-BFGS-B",
                                         hessian = TRUE, lower = c(0.001, K1, 0.5, -0.001),
                                         upper = c(3, K2, 10, 0.9))
              est_local_log_boot[j, ] =  local_out_log_boot$par
            }

            colnames(est_local_log_boot) = c("r", "K", "sigma2", "rho")
            est_local_log_list_p[[l]] = est_local_log_boot

          }


          # Estimates
          mat_est_r_MLE_log_p <- matrix(data = NA, nrow = q-P+1, ncol = B)
          for (j in 1:B) {
            mat_est_r_MLE_log_p[,j] = est_local_log_list_p[[j]][,1]
          }

          mat_est_K_MLE_log_p <- matrix(data = NA, nrow = q-P+1, ncol = B)
          for (j in 1:B) {
            mat_est_K_MLE_log_p[,j] = est_local_log_list_p[[j]][,2]
          }

          modified_RGR_log_p <- matrix(data = NA, nrow = q-P+1, ncol = B)
          for (i in 1:(q-P+1)) {
            numerator = mat_est_r_MLE_log_p[i,]*exp(-mat_est_r_MLE_log_p[i,]*t[i])*(mat_est_K_MLE_log_p[i,] - x0_mean_vals_log_p)
            denominator = x0_mean_vals_log_p + exp(-mat_est_r_MLE_log_p[i,]*t[i])*(mat_est_K_MLE_log_p[i,] - x0_mean_vals_log_p)
            modified_RGR_log_p[i,] = numerator/denominator
          }



          V_boot_log_p <- matrix(data = NA, nrow = q-P+1, ncol = B)
          for (i in 1:(q-P+1)) {
            V_boot_log_p[i, ] = log(1/modified_RGR_log_p[i,] - 1/ mat_est_r_MLE_log_p[i,])
          }

          test_stat_boot_log_p <- matrix(data = NA, nrow = q-P-1, ncol = B)
          for (i in 1:(q-P-1)) {
            test_stat_boot_log_p[i,] = V_boot_log_p[i+2,] - 2*V_boot_log_p[i+1,] + V_boot_log_p[i,]
          }

          sd_test_stat_log_p = numeric(q-P-1)
          for (i in 1:(q-P- 1)) {
            sd_test_stat_log_p[i] = sd(test_stat_boot_log_p[i,])
          }


          p_value_log_p = numeric(q-P-1)
          for (i in 1:(q-P- 1)) {
            p_value_log_p[i] = 1- pnorm(abs(test_stat_boot_log_p[i]/sd_test_stat_log_p[i]))
          }

          #------------------LOGISTIC BOOTSTRAP END--------------



          # ------------ Log non param below ------------------


          mat_est_r_MLE_log_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (j in 1:B_) {
            mat_est_r_MLE_log_np[,j] = est_local_log_list_np[[j]][,1]
          }

          mat_est_K_MLE_log_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (j in 1:B_) {
            mat_est_K_MLE_log_np[,j] = est_local_log_list_np[[j]][,2]
          }

          # Estimation of modified RGR and its distribution
          modified_RGR_log_np= matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (i in 1:(q-P+1)) {
            numerator = mat_est_r_MLE_log_np[i,]*exp(-mat_est_r_MLE_log_np[i,]*t[i])*(mat_est_K_MLE_log_np[i,] - x0_mean_vals_log_np)
            denominator = x0_mean_vals_log_np + exp(-mat_est_r_MLE_log_np[i,]*t[i])*(mat_est_K_MLE_log_np[i,] - x0_mean_vals_log_np)
            modified_RGR_log_np[i,] = numerator/denominator
          }

          # Computation of V
          V_boot_log_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (i in 1:(q-P+1)) {
            V_boot_log_np[i, ] = log(1/modified_RGR_log_np[i,] - 1/ mat_est_r_MLE_log_np[i,])
          }

          # Computation of test stat
          test_stat_boot_log_np = matrix(data = NA, nrow = q-P-1, ncol = B_)
          for (i in 1:(q-P-1)) {
            test_stat_boot_log_np[i,] = V_boot_log_np[i+2,] - 2*V_boot_log_np[i+1,] + V_boot_log_np[i,]
          }

          sd_test_stat_log_np = numeric(q-P-1)
          for (i in 1:(q-P- 1)) {
            sd_test_stat_log_np[i] = sd(test_stat_boot_log_np[i,])
          }

          p_value_log_np = numeric(q-P-1)
          for (i in 1:(q-P-1)) {
            p_value_log_np[i] = 1- pnorm(abs(test_stat_log[i]/sd_test_stat_log_np[i]))
          }

          # --------------------- Log non param above -----------------

          # von

          show_toast('Only a bit more',
                     text = 'Almost there..',
                     timer = 8000,
                     position = "center",
                     type = 'success')

          RGR_hat_VB = numeric(length = q-P+1)
          for(i in 1:(q-P+1)){
            RGR_hat_VB[i] = r_est_local_v[i,1]*((r_est_local_v[i,2]/mean(data[,i]))^(1/3))*(1- (mean(data[,i])/r_est_local_v[i,2])^(1/3))
          }



          var_RGR_hat_VB <- numeric(length = q-P+1)
          for(i in 1:(q-2)){
            A = r_cov_local_v[[i]][1,1]*((r_est_local_v[i,2]/mean(data[,i]))^(1/3)-1)^2
            B = r_cov_local_v[[i]][2,2]*(r_est_local_v[i,1]/(3*(mean(data[,i])*r_est_local_v[i,2]^2)^(1/3)))^2
            C = 2*((r_est_local_v[i,2]/mean(data[,i]))^(1/3)-1)*(r_est_local_v[i,1]/(3*(mean(data[,i])*r_est_local_v[i,2]^2)^(1/3)))*r_cov_local_v[[i]][1,2]
            var_RGR_hat_VB[i] = A+B+C
          }


          # Computation of V
          V_VB <- numeric(length = q-2)
          for(i in 1:(q-2)){
            V_VB[i] = log(1/RGR_hat_VB[i] + 1/r_est_local_v[i,1])
          }



          # Computation of Variance of V
          var_V_VB <- numeric(length = q-2)
          for(i in 1:(q-2)){
            A = (-r_est_local_v[i,1]/(RGR_hat_VB[i]*(r_est_local_v[i,1]+RGR_hat_VB[i])))^2*var_RGR_hat_VB[i]
            B = (-RGR_hat_VB[i]/(r_est_local_v[i,1]*(r_est_local_v[i,1]+RGR_hat_VB[i])))^2*r_cov_local_v[[i]][1,1]
            C = 2*(1/(r_est_local_v[i,1]+RGR_hat_VB[i])^2)*(r_cov_local_v[[i]][1,1]*(1-mean(data[,i])/r_est_local_v[i,2]) + (r_est_local_v[i,1]*mean(data[,i])/(r_est_local_v[i,2])^2)*r_cov_local_v[[i]][1,2])
            var_V_VB[i] = A + B + C
          }



          # Computation of the test statistic
          test_stat_VB <- numeric(length = q-P-1)
          var_test_stat_VB = numeric(length = q-P-1)
          for(i in 1:(q-P-1)){
            var_test_stat_VB[i] = var_V_VB[i+2] + 4*var_V_VB[i+1] + var_V_VB[i]
            test_stat_VB[i] = (V_VB[i+2]-2*V_VB[i+1] + V_VB[i])

          }





          # Computation of the p-value
          p_value_local_VB <- numeric(q-P-1)
          for(i in 1:(q-P-1)){
            p_value_local_VB[i] = 1- pnorm(test_stat_VB[i]/sqrt(var_test_stat_VB[i]))
          }



          # ------------ VB non param below ---------------------



          mat_est_r_MLE_VB_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (j in 1:B_) {
            mat_est_r_MLE_VB_np[,j] = est_local_VB_list_np[[j]][,1]
          }

          mat_est_K_MLE_VB_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (j in 1:B_) {
            mat_est_K_MLE_VB_np[,j] = est_local_VB_list_np[[j]][,2]
          }

          # Estimation of modified RGR and its distribution
          modified_RGR_VB_np= matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (i in 1:(q-P+1)) {
            numerator = mat_est_r_MLE_VB_np[i,]*exp(-mat_est_r_MLE_VB_np[i,]*(t[i]/3))*
              ((mat_est_K_MLE_VB_np[i,])^(1/3) - (x0_mean_vals_VB_np)^(1/3))
            denominator = (mat_est_K_MLE_VB_np[i,])^(1/3) + exp(-mat_est_r_MLE_VB_np[i,]*(t[i]/3))*
              ((x0_mean_vals_VB_np)^(1/3) - (mat_est_K_MLE_VB_np[i,])^(1/3))
            modified_RGR_VB_np[i,] = numerator/denominator
          }

          # Computation of V
          V_boot_VB_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (i in 1:(q-P+1)) {
            V_boot_VB_np[i, ] = log(1/modified_RGR_VB_np[i,] + 1/ mat_est_r_MLE_VB_np[i,])
          }

          # Computation of test stat
          test_stat_boot_VB_np = matrix(data = NA, nrow = q-P-1, ncol = B_)
          for (i in 1:(q-P-1)) {
            test_stat_boot_VB_np[i,] = V_boot_VB_np[i+2,] - 2*V_boot_VB_np[i+1,] + V_boot_VB_np[i,]
          }


          sd_test_stat_VB_np = numeric(q-P-1)
          for (i in 1:(q-P- 1)) {
            sd_test_stat_VB_np[i] = sd(test_stat_boot_VB_np[i,])
          }

          # Computation of the p-value
          p_value_VB_np = numeric(q-P-1)
          for (i in 1:(q-P- 1)) {
            p_value_VB_np[i] = 1- pnorm(abs(test_stat_VB[i]/sd_test_stat_VB_np[i]))
          }

          # --------------------- VB non param above -----------------


          show_toast('Only a bit more',
                     text = 'Getting your results..',
                     timer = 8000,
                     position = "center",
                     type = 'success')


          #gompertz p_val

          RGR_hat_gom = numeric(length = q-P+1)
          for(i in 1:(q-P+1)){
            RGR_hat_gom[i] = r_est_local_gom[i,2]*exp(-r_est_local_gom[i,1]*t[i])
          }

          # Computation of variance of Modified RGR using Delta method
          var_RGR_hat_gom = numeric(length = q-P+1)
          for(i in 1:(q-2)){
            A = r_cov_local_gom[[i]][1,1]*(-r_est_local_gom[i,2]*t[i]*exp(-r_est_local_gom[i,1]*t[i]))^2
            B = r_cov_local_gom[[i]][2,2]*(exp(-r_est_local_gom[i,1]*t[i]))^2
            C = 2*(-r_est_local_gom[i,2]*t[i]*exp(-r_est_local_gom[i,1]*t[i]))*(exp(-r_est_local_gom[i,1]*t[i]))*r_cov_local_gom[[i]][1,2]
            var_RGR_hat_gom[i] = A+B+C
          }

          # Computation of V
          V_gom = numeric(length = q-2)
          for(i in 1:(q-2)){
            V_gom[i] = log(1/RGR_hat_gom[i])
          }


          # Computation of Variance of V
          var_V_gom = numeric(length = q-2)
          for(i in 1:(q-2)){
            A = r_cov_local_gom[[i]][1,1]*(-t[i])^2
            B = r_cov_local_gom[[i]][2,2]*(1/r_est_local_gom[i,2])^2
            C = 2*(-t[i]/r_est_local_gom[i,2])*r_cov_local_gom[[i]][1,2]
            var_V_gom[i] = A + B + C
          }

          # Computation of the test statistic
          test_stat_gom = numeric(length = q-P-1)
          var_test_stat_gom = numeric(length = q-P-1)
          for(i in 1:(q-P-1)){
            var_test_stat_gom[i] = var_V_gom[i+2] + 4*var_V_gom[i+1] + var_V_gom[i]
            test_stat_gom[i] = (V_gom[i+2]-2*V_gom[i+1] + V_gom[i])
          }

          # Computation of the p-value
          p_value_local_gom = numeric(q-P-1)
          for(i in 1:(q-P-1)){
            p_value_local_gom[i] = 1- pnorm(abs(test_stat_gom[i]/sqrt(var_test_stat_gom[i])))
          }


          show_toast('Only a bit more',
                     text = 'Getting your plots..',
                     timer = 8000,
                     position = "center",
                     type = 'success')

          # ------------------ GOM non param below  ------------------

          mat_est_b_MLE_gom_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (j in 1:B_) {
            mat_est_b_MLE_gom_np[,j] = est_local_gom_list_np[[j]][,1]
          }

          mat_est_c_MLE_gom_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (j in 1:B_) {
            mat_est_c_MLE_gom_np[,j] = est_local_gom_list_np[[j]][,2]
          }

          # Estimation of modified RGR and its distribution
          modified_RGR_gom_np= matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (i in 1:(q-P+1)) {
            modified_RGR_gom_np[i,] = mat_est_b_MLE_gom_np[i,]*exp(-mat_est_c_MLE_gom_np[i,]*t[i])
          }

          # Computation of V
          V_boot_gom_np = matrix(data = NA, nrow = q-P+1, ncol = B_)
          for (i in 1:(q-P+1)) {
            V_boot_gom_np[i, ] = log(1/modified_RGR_gom_np[i,])
          }

          # Computation of test stat
          test_stat_boot_gom_np = matrix(data = NA, nrow = q-P-1, ncol = B_)
          for (i in 1:(q-P-1)) {
            test_stat_boot_gom_np[i,] = V_boot_gom_np[i+2,] - 2*V_boot_gom_np[i+1,] + V_boot_gom_np[i,]
          }


          sd_test_stat_gom_np = numeric(q-P-1)
          for (i in 1:(q-P- 1)) {
            sd_test_stat_gom_np[i] = sd(test_stat_boot_gom_np[i,])
          }

          # Computation of the p-value
          p_value_gom_np = numeric(q-P-1)
          for (i in 1:(q-P- 1)) {
            p_value_gom_np[i] = 1- pnorm(abs(test_stat_gom[i]/sd_test_stat_gom_np[i]))
          }


          # ------------------- Gom non param above ----------------------

          show_toast('Progress',
                     text = 'Done, We Appreciate Your Patience',
                     timer = 3000,
                     position = "center",
                     type = 'success')

          if (input$method_p == 'Parametric Bootstrap'){

            output$ap_plot <- renderPlot({
              plot(1:(q-P-1),p_value_log_p, col = "red", type = "b",
                   ylim = c(0,1), ylab = "p-value", lwd = 2, pch = 19,xlab = 'Time Points')
              lines(1:(q - P - 1), p_value_local_VB, col = "green", type = "b")
              lines(1:(q - P - 1), p_value_local_gom, col = "magenta", type = "b")

              # Add a single legend for all lines
              legend("topright", legend = c("Logistic P value Plot", "VB P value Plot", "Gompertz P value Plot"),
                     col = c("red", "green", "magenta"), lty = 1, lwd = 2, pch = 19)

              # Add a horizontal line at the specified tolerance level
              abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

            })

            observeEvent(input$zoom13,{

              showModal(modalDialog(

                renderPlot({
                  plot(1:(q-P-1),p_value_log_p, col = "red", type = "b",
                       ylim = c(0,1), ylab = "p-value", lwd = 2, pch = 19,xlab = 'Time Points')
                  lines(1:(q - P - 1), p_value_local_VB, col = "green", type = "b")
                  lines(1:(q - P - 1), p_value_local_gom, col = "magenta", type = "b")

                  # Add a single legend for all lines
                  legend("topright", legend = c("Logistic P value Plot", "VB P value Plot", "Gompertz P value Plot"),
                         col = c("red", "green", "magenta"), lty = 1, lwd = 2, pch = 19)

                  # Add a horizontal line at the specified tolerance level
                  abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)
                }, height = 500),
                easyClose = TRUE,
                size = "l",
                footer = NULL
              ))
            })

            output$ap_plot_r <- downloadHandler(
              filename = function() {

                if (input$format12 == "png") {
                  paste("p_val_comp-", Sys.Date(), ".png", sep = "")
                } else if (input$format12 == "pdf") {
                  paste("p_val_comp-", Sys.Date(), ".pdf", sep = "")
                }

              },
              content = function(file) {


                if (input$format12 == "png") {
                  png(file, width = 1100, height = 500)
                } else if (input$format12 == "pdf") {
                  pdf(file,  width = 11, height = 8.5)
                }

                plot(1:(q-P-1),p_value_log_p, col = "red", type = "b",
                     ylim = c(0,1), ylab = "p-value", lwd = 2, pch = 19,xlab = 'Time Points')
                lines(1:(q - P - 1), p_value_local_VB, col = "green", type = "b")
                lines(1:(q - P - 1), p_value_local_gom, col = "magenta", type = "b")

                # Add a single legend for all lines
                legend("topright", legend = c("Logistic P value Plot", "VB P value Plot", "Gompertz P value Plot"),
                       col = c("red", "green", "magenta"), lty = 1, lwd = 2, pch = 19)

                # Add a horizontal line at the specified tolerance level
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

                dev.off()
              } )}

          else if (input$method_p == 'Non-parametric Bootstrap'){
            output$ap_plot <- renderPlot({
              plot(1:(q-P-1), p_value_log_np, col = "red", type = "b", ylim = c(0,1), ylab = "p-value",
                   xlab= "Time Points", lwd = 2, pch = 19)
              lines(1:(q-P-1), p_value_VB_np, col = "green", type = "b")
              lines(1:(q-P-1), p_value_gom_np, col = "magenta", type = "b")

              # Add a single legend for all lines
              legend("topright", legend = c("Logistic P value Non - param Plot", "VB P value Non - param Plot", "Gompertz P value Non - param Plot"),
                     col = c("red", "green", "magenta"), lty = 1, lwd = 2, pch = 19)

              # Add a horizontal line at the specified tolerance level
              abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

            })

            observeEvent(input$zoom13,{

              showModal(modalDialog(

                renderPlot({
                  plot(1:(q-P-1), p_value_log_np, col = "red", type = "b", ylim = c(0,1), ylab = "p-value",
                       xlab= "Time Points", lwd = 2, pch = 19)
                  lines(1:(q-P-1), p_value_VB_np, col = "green", type = "b")
                  lines(1:(q-P-1), p_value_gom_np, col = "magenta", type = "b")

                  # Add a single legend for all lines
                  legend("topright", legend = c("Logistic P value Non - param Plot", "VB P value Non - param Plot", "Gompertz P value Non - param Plot"),
                         col = c("red", "green", "magenta"), lty = 1, lwd = 2, pch = 19)

                  # Add a horizontal line at the specified tolerance level
                  abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

                }, height = 500),
                easyClose = TRUE,
                size = "l",
                footer = NULL
              ))
            })

            output$ap_plot_r <- downloadHandler(
              filename = function() {

                if (input$format12 == "png") {
                  paste("p_val_comp-", Sys.Date(), ".png", sep = "")
                } else if (input$format12 == "pdf") {
                  paste("p_val_comp-", Sys.Date(), ".pdf", sep = "")
                }

              },
              content = function(file) {


                if (input$format12 == "png") {
                  png(file, width = 1100, height = 500)
                } else if (input$format12 == "pdf") {
                  pdf(file,  width = 11, height = 8.5)
                }

                plot(1:(q-P-1), p_value_log_np, col = "red", type = "b", ylim = c(0,1), ylab = "p-value",
                     xlab= "Time Points", lwd = 2, pch = 19)
                lines(1:(q-P-1), p_value_VB_np, col = "green", type = "b")
                lines(1:(q-P-1), p_value_gom_np, col = "magenta", type = "b")

                # Add a single legend for all lines
                legend("topright", legend = c("Logistic P value Non - param Plot", "VB P value Non - paramPlot", "Gompertz P value Non - param Plot"),
                       col = c("red", "green", "magenta"), lty = 1, lwd = 2, pch = 19)

                # Add a horizontal line at the specified tolerance level
                abline(h = input$tolerance, lwd = 2, col = "black", lty = 2)

                dev.off()
              } )

          }


          # here now












          # Your calculation code here
        })} else {
          # Data is not numeric, show a toast message
          show_toast('Caution ⚠️',
                     text = 'Data Must be all Numeric type',
                     timer = 5000,
                     position = "center",
                     type = 'info')
        }
    }})




  #data <- reactive({
  #  req(input$upload)
  # read.csv(input$upload$datapath)
  #})


  # <---------------------- Template Code --------------------------------------------->
  observeEvent(input$navbarID, {
    if(input$navbarID == "Home"){
      session$sendCustomMessage("background-color", "#2e4052")
    } else {
      session$sendCustomMessage("background-color", "#e6f4f1")
    }
  })



  observe({
    if (input$navbarID == 'Simulation Study') {
      showModal(modalDialog(
        includeHTML("www/intro_1.html"),
        easyClose = TRUE
      ))
    }
  })

  observe({
    if (input$navbarID == 'Real Data') {
      showModal(modalDialog(
        includeHTML("www/intro_2.html"),
        easyClose = TRUE
      ))
    }
  })

  #Download Handler



}
