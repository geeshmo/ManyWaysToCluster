library(ggplot2)
library(shiny)

source("mixture.R")

enableActionButton <- function(id,session, val) {
  session$sendCustomMessage(type="jsCode",
                            list(code= paste("$('#",id,"').prop('disabled',",tolower(toString(!val)) ,")"
                                             ,sep="")))
}

shinyServer(
  
  function(input, output, session) {
    
    observe({
      if (input$genBtn == 0) {
        enableActionButton("runBtn", session, FALSE)
      } else {
        enableActionButton("runBtn", session, TRUE)
      }
    })
    
    observe({
      if (input$runBtn  > 0) {
        updateTabsetPanel(session, inputId = "theTabs", selected = "runTab")
      }
    })
    
    observe({
      if (input$genBtn  > 0) {
        updateTabsetPanel(session, inputId = "theTabs", selected = "dataTab")
      }
    })
    
    run.algorithm <- reactive({
      if (input$runBtn ==0) {
        return(NULL)
      }
      
      updateTabsetPanel(session, inputId = "theTabs", selected = "runTab")
      
      #graph.update$cur.graph <-NULL
      g.d <- isolate(generated.data())
      num.clusters <- isolate(input$num.clusters)
      num.iters <- isolate(input$num.iters)
      alg.name <- isolate(getAlgName())
      
      alg.out <- switch(alg.name,
                        k_means = withProgress(message = 'Running K-means', value = 0, run.kmeans(g.d, isolate(input$k_means.num.clusters), 
                                                                                                  isolate(input$k_means.num.iters), incProgress)),
                        gmm_gibbs = withProgress(message = 'Running Block Gibbs', value = 0,  run.block.gibbs(g.d, isolate(input$gibbs.num.clusters),
                                                                                                              isolate(input$gibbs.num.iters), incProgress)),
                        gmm_collapsed_gibbs = withProgress(message = 'Running Collapsed Gibbs', value = 0, 
                                                           run.collapsed.gibbs(g.d, isolate(input$collapsed.gibbs.num.clusters),
                                                                               isolate(input$collapsed.gibbs.num.iters), incProgress)),
                        gmm_infinite_mixture = withProgress(message = 'Running Infinite Mixture Collapsed Gibbs', value = 0, 
                                                            run.infinite.mixture.collapsed.gibbs(g.d, isolate(input$infinite.mixture.num.iters), incProgress)),
                        subspace_clustering = withProgress(message = 'Running Subspace Clustering', value = 0, run.lac(g.d, isolate(input$lac.num.clusters),
                                                                                                                       isolate(input$lac.num.iters), incProgress))
      )
      return(alg.out)})
    
    
    output$algAnimation <- renderUI({
      alg.run.res <- run.algorithm()
      file.name <- generate.animation(alg.run.res)
      #iframe.tag <- paste("<iframe src=\"", file.name,"\" width=\"100%\"  frameborder=\"0\" height=\"100%\"></iframe>", sep="")
      tags$iframe(src=file.name, style="width: 100%; height: 100vh; margin: 0; padding: 0; overflow-y: hidden; border: 0")
      #HTML(iframe.tag)
      
      
    })
    
    generated.data <- reactive({
      if (input$genBtn == 0) {
        return(NULL)
      }
      
      boundary <- 20
      isolate(generate.data(input$num.clusters, input$num.points, boundary, input$uncor.dims))
    })
    
    # Update data plot every time generated data changes
    dataPlot <- reactive({
      g.d <- generated.data()
      
      if (!is.null(g.d)) {
        out <- plot.generated.data(g.d)
        
      } else { 
        out <- NULL
      }
      return(out)
    })
    
    # Render the input data
    output$dataPlot <- renderPlot({
      print(dataPlot())
    })
    
    getAlgName <- reactive({
      input$algorithm.to.run
    })
    
    output$algParams <- renderUI({
      algName <- getAlgName()
      out <- switch(algName,
                    k_means = list(
                      numericInput(inputId = "k_means.num.clusters", label = "Suggest number of clusters", 5, min = 1, max = 10),
                      numericInput(inputId = "k_means.num.iters", label = "Number of iterations", 20, min = 1, max = 100)
                    ),
                    gmm_gibbs = list(
                      numericInput(inputId = "gibbs.num.clusters", label = "Suggest number of clusters", 5, min = 1, max = 10),
                      numericInput(inputId = "gibbs.num.iters", label = "Number of iterations", 40, min = 1, max = 100)
                    ),
                    gmm_collapsed_gibbs = list(
                      numericInput(inputId = "collapsed.gibbs.num.clusters", label = "Suggest number of clusters", 5, min = 1, max = 10),
                      numericInput(inputId = "collapsed.gibbs.num.iters", label = "Number of iterations", 40, min = 1, max = 100)
                    ),
                    gmm_infinite_mixture = list(
                      numericInput(inputId = "infinite.mixture.num.iters", label = "Number of iterations", 40, min = 1, max = 100)
                    ),
                    subspace_clustering = list(
                      numericInput(inputId = "lac.num.clusters", label = "Suggest number of clusters", 5, min = 1, max = 10),
                      numericInput(inputId = "lac.num.iters", label = "Number of iterations", 40, min = 1, max = 100)
                    )
      )
      return(out)
    })
  }
)
