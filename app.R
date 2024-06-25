#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)
library(dplyr)

ui <- navbarPage(
  "Exploring Change in Estimate Criterion",   
  tabPanel("Summary", 
           
           h2("Background"),
           
           p("There are multiple methods to select potential covariates when building a statistical models. Change in Estimate Criterion (CEC) are a data-driven method to select covariates (using forward or backward deletion) based on their impact on the accuracy of the model."),
           
           p("In his 2014 publication, Dr. Paul H. Lee performed a series of simulations to assess the validity of the commonly-used 10% cut-off for CECs He performed simulations for both linear and logistic regressions, finding a 95% significance cut-off for the change in estimate criterion."),
           
           p("For a more detailed discussion on CECs and further background see Dr. Lee's publication. "),
           
           h2("Purpose"),
           
           p("This app is meant to be a simple demonstration of how sample size, variance, and effect size can influence the `optimal` threshold to detect a change in estimate. The Linear Regression and Logistic Regression pages demonstrate changes in CECs for each technique."),
           
           p("Note that the original publication used 10,000 simulations to estimate thresholds. To minimize buffering I have lowered the number of simulations run in this simple demonstration."),
           
           h2("Reference"),
           
           p("Lee P. H. (2014). Is a cutoff of 10% appropriate for the change-in-estimate criterion of confounder identification?.", em("Journal of epidemiology, 24"), "(2), 161–167. https://doi.org/10.2188/jea.je20130062"),
           
           p("For a full copy of the files used to create this app see: https://github.com/Sebastian-Santana-Ort/Shiny_Simulation")
           ),
  tabPanel("Linear Regression", 
           sidebarLayout(
             sidebarPanel(
               sliderInput("sample_size",
                           "Sample Size:",
                           min = 1,
                           max = 1000,
                           value = 500,
                           step = 100),
               
               sliderInput("var_e",
                           "Variance of Error:",
                           min = 0.5,
                           max = 4,
                           value = 1,
                           step = 0.1),
               
               sliderInput("effect",
                           "Effect Size:",
                           min = 0,
                           max = 1,
                           value = 0.5,
                           step = 0.10),
               
               sliderInput("corr",
                           "Exposure-confounder correlation:",
                           min = 0,
                           max = 1,
                           value = 0.1,
                           step = 0.01),
               
               numericInput("sim_size",
                            "Simulation Size:",
                            value = 1000,
                            min = 1000)
             ),
             
             mainPanel(
               verbatimTextOutput("result_linear"),
               plotOutput("histogram_linear")
             )
           )),
  tabPanel("Logistic Regression", 
           sidebarLayout(
             sidebarPanel(
               sliderInput("sample_size_log",
                           "Sample Size:",
                           min = 1,
                           max = 1000,
                           value = 500,
                           step = 100),
               
               sliderInput("corr_log",
                           "Correlation Between Predictor and Confounder:",
                           min = .1,
                           max = 1,
                           value = .2,
                           step = 0.1),
               
               sliderInput("effect_log",
                           "Odds Ratio:",
                           min = 1,
                           max = 4,
                           value = 1.5,
                           step = .5),
               
               numericInput("sim_size_log",
                            "Simulation Size:",
                            value = 1000,
                            min = 1000)
             ),
             
             mainPanel(
               verbatimTextOutput("result_log"),
               plotOutput("histogram_log")
             )
           ))
)
server <- function(input, output, session) {

estimate_lin <- reactive({
    sim_size <- input$sim_size
    sample_size <- input$sample_size
    effect <- input$effect
    var_e <- input$var_e
    cor_x_z <- input$corr
    
    estimate_lin <- numeric(sim_size)
    set.seed(42)
    for (i in 1:sim_size) {
      x <- rnorm(sample_size)
      z <- cor_x_z * x + sqrt(1 - cor_x_z^2) * rnorm(sample_size)
      e <- rnorm(sample_size) * var_e
      y <- effect * x + z + e
      
      reg1 <- lm(y ~ x)
      test1 <- reg1$coefficients[2]
      
      reg2 <- lm(y ~ x + z)
      test2 <- reg2$coefficients[2]
      
      ratio <- test2 / test1
      if (ratio < 1) ratio <- 2 - ratio
      estimate_lin[i] <- ratio
    }
    estimate_lin
  })
  
output$result_linear <- renderPrint({
  estimate_lin <- estimate_lin()
  paste("95% Cut Off for Change in Estimate Criterion is: ", round((quantile(estimate_lin, 0.95) - 1) * 100, 3), "%")
})


output$histogram_linear <- renderPlot({
  estimate_lin <- estimate_lin()
  df <- data.frame(estimate_lin = estimate_lin)
  df_filtered <- df %>% filter(estimate_lin <= quantile(estimate_lin, 0.99))
  ggplot(df_filtered, aes(x = estimate_lin)) +
    geom_histogram(bins = 200, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
    ggtitle("Distribution of Changes in Estimates Across Simulations") +
    labs(caption = "Estimates above 99% were removed for the sake of simplicity in this visual.",
         x = "CEC Value",
         y = "Frequency") +
    theme_classic(base_size = 18)
})


estimate_log <- reactive({
  sim_size <- input$sim_size_log
  sample_size <- input$sample_size_log
  OR <- input$effect_log
  cor_x_z <- input$corr_log ## Exposure-confounder correlation
  
  estimate <- 1:sim_size
  ## Simulation starts
  set.seed(42)
  for (i in 1:sim_size){
    x <- rnorm(sample_size)
    z <- cor_x_z*x+sqrt(1-cor_x_z*cor_x_z)*rnorm(sample_size)
    p <- exp(log(OR)*x) / (1+exp(log(OR)*x))
    
    y <- 1:sample_size
    for (j in 1:sample_size){
      y[j] <- sample(0:1,1,rep=TRUE,prob=c(1-p[j],p[j]))
    }
    
    reg <- glm(y~x, family = binomial)
    test1 <- exp(reg$coefficients[2])
    
    reg <- glm(y~x+z, family = binomial)
    test2 <- exp(reg$coefficients[2])
    
    ratio <- test2/test1
    if (ratio<1) ratio = 2-ratio
    estimate[i] <- ratio
  }
  estimate
})

output$result_log <- renderPrint({
  estimate_log <- estimate_log()
  paste("95% Cut Off for Change in Estimate Criterion is: ", round((quantile(estimate_log, 0.95) - 1) * 100, 3), "%")
})


output$histogram_log = renderPlot({
  estimate_log = estimate_log()
  df <- data.frame(estimate_log = estimate_log)
  df_filtered <- df %>% filter(estimate_log <= quantile(estimate_log, 0.99))
  ggplot(df_filtered, aes(x = estimate_log)) +
    geom_histogram(bins = 200, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
    ggtitle("Distribution of Changes in Estimates Across Simulations") +
    labs(caption = "Estimates above 99% were removed for the sake of simplicity in this visual.",
         x = "CEC Value",
         y = "Frequency") +
    theme_classic(base_size = 18)
  })
  
}

shinyApp(ui = ui, server = server)
