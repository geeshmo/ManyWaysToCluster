library(shiny)

shinyUI(
  
  pageWithSidebar(
    headerPanel("Clustering"),
    
    sidebarPanel(
      
      # The number of clusters to generate  
      numericInput("num.clusters", "Number of clusters (1-10)", 5, min = 1, max = 10),
      
      # The total number of points to generate  
      numericInput("num.points", "Number of points (100-2000)", 1000, min = 100, max = 2000),
      
      checkboxInput("uncor.dims", "Uncorrelated dimensions", value = FALSE),
      br(),
      # Generate data
      actionButton("genBtn", label="Generate Data"),
      br(),
      br(),
      selectInput(
        inputId = "algorithm.to.run",
        label   = "Select an algorithm to run",
        choices = list(
          "K-Means" = "k_means",
          "GMM - Gibbs" = "gmm_gibbs",
          "GMM - Collapsed Gibbs" = "gmm_collapsed_gibbs",
          "GMM - Infinite Mixture" = "gmm_infinite_mixture",
          "Subspace clustering" = "subspace_clustering"
        )),
      
      uiOutput(outputId = "algParams"),
      actionButton("runBtn", label="Run"),
      
      uiOutput(outputId = "xLimSlider"),
      tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",
           function(message) {
              console.log(message)
              eval(message.code);
           }
          );
         ')))
    ),
    
    mainPanel(
      tabsetPanel(id ="theTabs",
                  tabPanel("Input data",
                           plotOutput(outputId = "dataPlot"), value = "dataTab"),
                  tabPanel("Algorithm output", uiOutput(outputId = "algAnimation"), value = "runTab")
      )
    )
    
  )
)
