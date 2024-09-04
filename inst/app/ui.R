
ui <- tagList(tags$head(
  tags$script("
      Shiny.addCustomMessageHandler('background-color', function(color) {
        document.body.style.backgroundColor = color;
      });
    ")
),
tags$head(
  tags$style(
    HTML(".custom-table {
             width: 80%; /* Adjust the width as needed */
             height: 300px; /* Adjust the height as needed */
           }")
  )
),


navbarPage(
  windowTitle = 'GEM-R',

  title = HTML('<span style="color: white;font-size = 10px;">G<sub>p</sub>EM-R</span>'),

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "stylefinal.css")
  ),
  id = "navbarID",




  tabPanel('Home',
           tags$head(
             tags$link(rel = "stylesheet", type = "text/css", href = "www/stylehome.css")
           ),
           tags$div(
             includeHTML("www/Home.html")
           ),
           # Include the JavaScript file
           tags$script(type = "text/javascript", src = "www/homeapp.js")),
  tabPanel(
    'Simulation Study',
    id = 'Tab_1',
    sidebarLayout(
      sidebarPanel(

        actionButton('button1',div(id='btntext','📉 Get Estimates'), style = "font-size: 20px; padding: 10px;,width : 100%",class='custom-btn'),
        br(),
        br(),
        br(),

        sliderInput(
          inputId = "t_range",
          label = div(id = 'label',"⏲️ Time Points:"),
          min = 0,
          max = 50,
          value = c(1,10),
          step = 1
        ),
        helpText('Select Time Points',class = "help-text" ),
        br(),

        numericInput("trajectory", div(id='label',"🔢 No. of Independent Trajectories (n):"), value = 10,min = 1),
        helpText('Select No. of Independent Trajectories (n)',class = "help-text"),

        br(),

        numericInput("window",div(id='label',"🔢 Window Size :"),value = 3,min = 0),
        helpText("Select a Window size for Local Estimation",class = "help-text"),

        br(),

        selectInput(
          "model",
          div(id='label',"⚙️ Choose a Model :"),
          choices = c("Logistic Model","Exponential Model", "Theta - Logistic Model", "Von-Bertallanfy Model", "Gompertz Model"),selected = "Logistic Model"),
        helpText('Select a Model You Want to Use',class = "help-text"),
        class = "none",

        br(),


        conditionalPanel(
          condition = "input.model == 'Logistic Model'",
          numericInput("logi_r", div(id = 'label'," r :"), value = 0.2,min = 0.001),
          helpText('Select a Value of r ',class = "help-text"),

          br(),

          numericInput("logi_k", div(id = 'label'," K :"), value = 100,min = 0.001),
          helpText('Select a K Value',class = "help-text"),

          br(),

          numericInput("logi_s", div(id = 'label'," Sigma (σ):"), value = 2,min = 0.001),
          helpText('Select a Sigma (σ) Value',class = "help-text"),

          br(),

          numericInput("logi_rho", div(id = 'label'," Rho (ρ) :"), value = 0.5),
          helpText('Select a Rho (ρ) Value',class = "help-text"),


          br(),

          numericInput("logi_in", div(id = 'label'," Initial Value (x0) :"), value = 10,min = 0.001),
          helpText('Select a Theta Value',class = "help-text")),


        conditionalPanel(
          condition = "input.model == 'Exponential Model'",
          numericInput("expo_r", div(id = 'label'," r :"), value = 0.1,min = 0.001),
          helpText('Select Value of Exponential Parameter r ',class = "help-text"),

          br(),


          numericInput("expo_rho",div(id = 'label'," Rho (ρ) :"), value = 0.5),
          helpText('Select a Rho (ρ) Value',class = "help-text"),

          br(),

          numericInput("expo_s",div(id = 'label'," Sigma (σ) :"), value = 2,min = 0.001),
          helpText('Select a Sigma Value',class = "help-text"),

          br(),


          numericInput("expo_in", div(id = 'label'," Initial Point (x0) :"), value = 1,min = 0.001),
          helpText('Set an Initial Point',class = "help-text")
        ),


        conditionalPanel(
          condition = "input.model == 'Theta - Logistic Model'",

          numericInput("theta_r", div(id = 'label'," r :"), value = 0.2,min = 0.001),
          helpText('Select Value of r ',class = "help-text"),
          br(),

          numericInput("theta_th",div(id = 'label'," Theta :"), value = 1.2,min = 0.001),
          helpText('Select a Theta Value between (0.1 to 2) only',class = "help-text"),

          br(),


          numericInput("theta_k",div(id = 'label'," K :"), value = 100,min = 0.001),
          helpText('Select a K Value',class = "help-text"),
          br(),

          numericInput("theta_rho", div(id = 'label'," Rho (ρ) :"), value = 0.5),
          helpText('Select a Rho (ρ) Value ',class = "help-text"),

          br(),

          numericInput("theta_in",div(id = 'label'," Initial Value (x0) :"), value = 10,min = 0.001),
          helpText('Select a Theta Value',class = "help-text"),

          br(),

          numericInput("theta_s", div(id = 'label'," Sigma (σ) :"), value = 2,min = 0.001),
          helpText('Select a Sigma (σ) Value',class = "help-text")),


        conditionalPanel(
          condition = "input.model == 'Von-Bertallanfy Model'",

          numericInput("von_r", div(id = 'label'," r :"), value = 0.3,min = 0.001),
          helpText('Select Value of r ',class = "help-text"),

          br(),

          numericInput("von_k", div(id = 'label'," K :"), value = 100,min = 0.001),
          helpText('Select a Value of K',class = "help-text"),

          br(),

          numericInput("von_rho", div(id = 'label'," Rho (ρ) :"), value = 0.5),
          helpText('Select a Rho (ρ) Value',class = "help-text"),

          br(),


          numericInput("von_in", div(id = 'label'," Initial Value (x0) :"), value = 10,min = 0.001),
          helpText('Select a Theta Value',class = "help-text"),

          br(),

          numericInput("von_s", div(id = 'label'," Sigma (σ) :"), value = 2,min = 0.001),
          helpText('Select a Sigma (σ) Value',class = "help-text")
        ),



        conditionalPanel(
          condition = "input.model == 'Gompertz Model'",

          numericInput("gom_b", div(id = 'label'," b :"), value = 0.3,min = 0.001),
          helpText('Select a Value of b ',class = "help-text"),

          br(),

          numericInput("gom_c", div(id = 'label'," c :"), value = 0.2,min = 0.001),
          helpText('Select a c Value',class = "help-text"),

          br(),

          numericInput("gom_rho", div(id = 'label'," Rho (ρ) :"), value = 0.5),
          helpText('Select a Rho (ρ) Value',class = "help-text"),

          br(),

          numericInput("gom_in", div(id = 'label'," Initial Value (x0) :"), value = 10,min = 0.001),
          helpText('Select a Theta Value',class = "help-text"),

          br(),

          numericInput("gom_s", div(id = 'label'," Sigma (σ) :"), value = 2,min = 0.001),
          helpText('Select a Sigma (σ) Value',class = "help-text")

        ),

        prettyRadioButtons("radio_button", div(id = 'label',"👆 Method Selected to Compute ISRP (Default) : "), choices = c("Local Maximizaton"),animation = 'pulse',selected ="Local Maximizaton" ),
        helpText('A method is selected by Default',class = "help-text"),class = 'sidebar'),



      mainPanel(
        fluidRow(
          column(12,div(id = 'title',textOutput('model_name')),class = 'name_m')
        ),
        br(),
        fluidRow(
          column(9,div(id='plot',plotOutput('plot1')
          ),class = 'columns'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','This plot is Size Profile of Data'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format",
                         choiceNames = c("png",'pdf',"csv (data)"),
                         choiceValues = c('png','pdf','csv'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('plot_down','📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),

        br(),
        fluidRow(
          column(5,'Global Estimates Table',br(),
                 tableOutput("table1"),class = 'colu1'),
          column(5,'Covariance Matrix',
                 conditionalPanel(
                   condition = "input.model == 'Theta - Logistic Model'",
                   paste(' '),
                   paste(' '),
                   HTML('<span style="color: white; font-size = medium;">(Zoom to check all Estimates)</span>'),
                   br()
                 ),
                 tableOutput('covM'),div(id = 'downbtn',
                                         dropdown(
                                           radioGroupButtons(
                                             inputId = "format1_1",
                                             choiceNames = c('csv'),
                                             choiceValues = c('csv'),
                                             selected = "csv",
                                             direction = "vertical"
                                           ),
                                           downloadButton('c_mat','📟'),
                                           size = "s",
                                           icon = ("Download Covariance Matrix Here"),
                                           up = TRUE
                                         )),class = 'colu'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom2",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','Global Estimates of Parameters and its Covariance Matrix'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format1_2",
                         label = "Select type",
                         choiceNames = c('csv'),
                         choiceValues = c('csv'),
                         selected = "csv",
                         direction = "vertical"
                       ),
                       downloadButton('down_global','📟'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,'Local Estimates Table (Head Values) (Zoom to check full Table)',tableOutput('table2'),paste(' '),paste('Click Here to get Covariance Matrix of Local Estimates 👉'),
                 actionButton('Cov_M','Covariance table'),class = 'columns'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom3",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','Local Estimates of Parameters and its Covariance Matrix'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format2",
                         label = "Select type",
                         choiceNames = c('csv'),
                         choiceValues = c('csv'),
                         selected = "csv",
                         direction = "vertical"
                       ),
                       downloadButton('local_est','📟'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,'Local Estimates Plot',
                 div(
                   style = "display: flex; flex-wrap: wrap;",
                   div(
                     style = "flex: 1;",
                     plotOutput("plot_1")
                   ),
                   div(
                     style = "flex: 1;",
                     plotOutput("plot_2")
                   )
                 ),class = 'columns_plot'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom4",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext_plot','Plot(s) of Local Estimates'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format3",
                         choiceNames = c("png","pdf"),
                         choiceValues = c('png','pdf'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('l_plot_d','📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol_plot',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,'Local Estimates Plot',
                 br(),
                 div(
                   style = "display: flex; flex-wrap: wrap;",
                   div(
                     style = "flex: 1;",
                     plotOutput("plot_3")
                   ),
                   div(
                     style = "flex: 1;",
                     plotOutput("plot_4")
                   )
                 ),class = 'columns_plot'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom5",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext_plot','Plot(s) of Local Estimates'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format4",
                         choiceNames = c("png","pdf"),
                         choiceValues = c('png','pdf'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('l_plot_d2','📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol_plot',align = 'center')
        ),
        br(),
        conditionalPanel(
          condition = "input.model == 'Theta - Logistic Model'",
          fluidRow(
            column(9, 'Local Estimates Plot', br(),
                   div(
                     style = "display: flex; flex-wrap: wrap;",
                     div(
                       style = "flex: 1;",
                       plotOutput("plot_5")
                     )
                   ),
                   class = 'columns_plot'
            ),
            column(1,
                   div(id = 'zoombtn',
                       actionBttn(
                         inputId = "zoom6",
                         label = 'Zoom',
                         icon = icon("search-plus", class = "opt"),
                         style = "fill",
                         color = "success",
                         size = "s"
                       )
                   ),
                   div(id = 'desctext_plot', 'Plot of Local Estimates'),
                   div(id = 'downbtn',
                       dropdown(
                         radioGroupButtons(
                           inputId = "format5",
                           choiceNames = c("png", "pdf"),
                           choiceValues = c('png', 'pdf'),
                           selected = 'png',
                           direction = "vertical"
                         ),
                         downloadButton('l_plot_d3', '📈'),
                         size = "s",
                         icon = icon("download", class = "opt"),
                         up = TRUE
                       )
                   ),
                   class = 'smallCol_plot',
                   align = 'center'
            )
          )
        )

      )
    )
  ),
  tabPanel(
    'Real Data',
    sidebarLayout(
      sidebarPanel(

        actionButton('button2',div(id='btntext','📉 Get Estimates'), style = "font-size: 20px; padding: 10px;,width : 100%",class='custom-btn'),
        br(),
        br(),
        br(),


        fileInput("upload", div(id = 'label'," 📅 Upload File"), accept = c(".csv", ".txt")),
        helpText('Upload Your File (csv, txt format only)',class = "help-text"),
        class = "",
        br(),

        numericInput("window2",div(id='label',"🔢 Window Size :"),value = 3),
        helpText("Select a Window size for Local Estimation",class = "help-text"),

        br(),


        selectInput(
          "model1",
          div(id='label',"⚙️ Choose a Model :"),
          choices = c("Logistic Model", "Von-Bertallanfy Model", "Gompertz Model"),selected = "Logistic Model"
        ),
        helpText('Select a Model You Want to Use',class = "help-text"),
        class = "",

        conditionalPanel(
          condition = "input.model1 == 'Logistic Model'",
          numericInput(inputId = "prob", div(id = 'label'," r :"), value = 0.7,min = 0.001),

          helpText('Select a Value of r ',class = "help-text")),

        conditionalPanel(
          condition = "input.model1 == 'Von-Bertallanfy Model'",

          numericInput("von_r1", div(id = 'label'," r :"), value = 0.8,min = 0.001),
          helpText('Select Value of r ',class = "help-text")),

        conditionalPanel(
          condition = "input.model1 == 'Gompertz Model'",

          numericInput("gom_b1", div(id = 'label'," b :"), value = 0.3,min = 0.001),
          helpText('Select a Value of b ',class = "help-text"),

          br(),

          numericInput("gom_c1", div(id = 'label'," c :"), value = 0.2,min = 0.001),
          helpText('Select a c Value',class = "help-text")),

        prettyRadioButtons('tolerance',div(id = 'label',"🔢 Select a Tolerance level for P value"),choices = c('0.01','0.05','0.1'),animation = 'pulse',selected = '0.05'),
        helpText('Select a Tolerance level for P value calculation',class = 'help-text'),


        prettyRadioButtons("radio_button", div(id = 'label',"👆 Method Selected to Compute ISRP (Default) : "), choices = c("Local Maximizaton"),animation = 'pulse',selected ="Local Maximizaton" ),
        helpText('A method is selected by Default',class = "help-text"),

        br(),
        prettyRadioButtons('method_p',div(id='label','👆 Select a p value calculation method'),choices = c('Parametric Bootstrap','Non-parametric Bootstrap'),animation = 'pulse',selected = 'Parametric Bootstrap'),
        helpText('Select a method',class = "help-text"),
        actionButton('button3',div(id='btntext_2','📉 Compare p value Plots'), style = "font-size: 20px; padding: 10px;,width : 100%",class='custom-btn'),
        br(),
        helpText('To Compare P value plot generated by all model',class = 'help-text'),class = 'sidebar'),


      mainPanel(
        fluidRow(
          column(12,div(id = 'title',textOutput('model_name2')),class = 'name_m')
        ),
        br(),

        fluidRow(
          column(8,div(id = 'datatable',dataTableOutput('data')), class = 'columns'
          ),
          column(1,div(id = 'desctext','Your Uploaded Data will be displayed Here'),class='smallCol',align = "center")),

        br(),
        fluidRow(
          column(9,
                 div(id = 'plot',plotOutput('data_p')),class = 'columns'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom7",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','Size Profile of Uploaded Data'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format6",
                         choiceNames = c("png", "pdf"),
                         choiceValues = c('png', 'pdf'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('r_data_down', '📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),
        br(),
        fluidRow(
          column(5,'Global Estimates Table',br(),
                 tableOutput("table1_r"),class = 'colu1'),
          column(5,'Covariance Matrix',
                 tableOutput('covM_r'),div(id = 'downbtn',
                                           dropdown(
                                             radioGroupButtons(
                                               inputId = "format7_1",
                                               choiceNames = c('csv'),
                                               choiceValues = c('csv'),
                                               selected = "csv",
                                               direction = "vertical"
                                             ),
                                             downloadButton('c_mat_r','📟'),
                                             size = "s",
                                             icon = ("Download Covariance Matrix Here"),
                                             up = TRUE
                                           )),class = 'colu'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom8",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','Global Estimates of Parameters and its Covariance Matrix'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format7_2",
                         label = "Select type",
                         choiceNames = c('csv'),
                         choiceValues = c('csv'),
                         selected = "csv",
                         direction = "vertical"
                       ),
                       downloadButton('down_global_r','📟'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,div(id = 'plot',plotOutput('com_plot')),class = 'columns'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom9",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','Visualization of Fitted Curve'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format8",
                         choiceNames = c("png", "pdf"),
                         choiceValues = c('png', 'pdf'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('com_plot_r', '📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,'Local Estimates Table (Head Values) (Zoom to check full Table)',tableOutput('table2_r'),paste(' '),paste('Click Here to get Covariance Matrix of Local Estimates 👉'),
                 actionButton('Cov_M_r','Covariance table'),class = 'columns'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom10",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','Local Estimates of Parameters and its Covariance Matrix'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format9",
                         label = "Select type",
                         choiceNames = c('csv'),
                         choiceValues = c('csv'),
                         selected = "csv",
                         direction = "vertical"
                       ),
                       downloadButton('local_est_real','📟'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,'Local Estimates Plot',
                 div(
                   style = "display: flex; flex-wrap: wrap;",
                   div(
                     style = "flex: 1;",
                     plotOutput("plot_1_r")
                   ),
                   div(
                     style = "flex: 1;",
                     plotOutput("plot_2_r")
                   )
                 ),class = 'columns_plot'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom11",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext_plot','Plot(s) of Local Estimates'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format10",
                         choiceNames = c("png","pdf"),
                         choiceValues = c('png','pdf'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('l_plot_d_real','📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol_plot',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,div(id = 'plot',plotOutput('p_plot')),class = 'columns'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom12",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','P-value Plot of Selected Model'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format11",
                         choiceNames = c("png", "pdf"),
                         choiceValues = c('png', 'pdf'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('p_plot_r', '📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        ),
        br(),
        fluidRow(
          column(9,div(id = 'plot',plotOutput('ap_plot')),class = 'columns'),
          column(1,
                 div(id = 'zoombtn',
                     actionBttn(
                       inputId = "zoom13",
                       label = 'Zoom',
                       icon = icon("search-plus", class = "opt"),
                       style = "fill",
                       color = "success",
                       size = "s"
                     )),
                 div(id = 'desctext','Comparison Plot of P-value of all Models'),
                 div(id = 'downbtn',
                     dropdown(
                       radioGroupButtons(
                         inputId = "format12",
                         choiceNames = c("png", "pdf"),
                         choiceValues = c('png', 'pdf'),
                         selected = 'png',
                         direction = "vertical"
                       ),
                       downloadButton('ap_plot_r', '📈'),
                       size = "s",
                       icon = icon("download", class = "opt"),
                       up = TRUE
                     ))
                 ,class='smallCol',align = 'center')
        )



      )))

))



