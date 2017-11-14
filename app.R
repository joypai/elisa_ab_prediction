library(shiny)
library(shinythemes)
library(ggplot2)
library(nplr)
library(dplyr)
library(plotly)
library(openxlsx)
library(reshape2)
library(DT)

# UI
ui <- shinyUI(fluidPage(
  theme = shinytheme("simplex"),
  # Application title
  titlePanel("Eliza Curve Fitting"),
  
  # Sidebar with parameters
  sidebarLayout(
    sidebarPanel(width=3,
      fileInput("input_file", "Import data:", accept=c(".xlsx"))
    ),
    
    # Show neutralization plot
    mainPanel(width=9,
              tabsetPanel(
                tabPanel("curve plot",
                  h3(textOutput("no_data_message")),
                  plotlyOutput("curve_plot"),
                  h3(textOutput("new_values_header")),
                  dataTableOutput("od_table"),
                  h3(textOutput("new_preds_header")),
                  dataTableOutput("pred_conc_table")
                ),
                tabPanel("results",
                  dataTableOutput("full_table")
                  )
              )
    )
  )
))

# predict new concentrations
predict_y = function(params, x) {
  params <- as.numeric(params)
  y <- params[1] + (params[2] - params[1]) / (1+ 10^(params[4]*(params[3]-x))) ^ params[5]
  return(y)
}

predict_x = function(params, y) {
  params <- as.numeric(params)
  #c - log10(((d-a)/(y-a))-1)/b
  
  x <- params[3] - log10(((params[2]-params[1])/(y-params[1]))-1)/params[4]
  #x <- params[3] * (((params[2]-params[1])/(y-params[1]))-1)^(1/params[4])
  return(x)
}

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
  fit_model_and_predict <- reactive({
    results <- list()
    print(input$input_file$datapath)
    contents <- read.xlsx(input$input_file$datapath)
    
    block_starts <- grep("Raw Data", contents[,2])
    sample_ids <- as.character(contents[(block_starts+10),2:13])
    
    data <- apply(contents[(block_starts+2):(block_starts+9),2:13], 2, as.numeric)
    
    sample_dilutions <- as.numeric(contents[(block_starts+3):(block_starts+8),15])
    sc_concentrations <- as.numeric(contents[(block_starts+3):(block_starts+14),16])
    standard_curve_data <- data.frame(conc=sc_concentrations+0.0001, od=apply(data[1:2,], 2, mean))
    
    new_data <- data[3:8,]
    colnames(new_data) <- sample_ids
    rownames(new_data) <- paste("dilution",sample_dilutions,sep="")
    head(new_data)
    new_data_long <- melt(new_data, varnames=c("dilution","sample_id"), value.name="OD")

    # fit model using standard curve data
    np1 <- nplr(standard_curve_data$conc, standard_curve_data$od, useLog = TRUE, npars=4)
    par <- getPar(np1)
    
    #pred_conc <- predict_y(par$params, log10(new_data+0.0001))
    #pred_conc <- 10^predict_x(par$params, new_data)

    predicted_x <- getEstimates(np1, new_data_long$OD)$x
    pred_conc <- as.data.frame(matrix(predicted_x, nrow=6, ncol=12))
    rownames(pred_conc) <- rownames(new_data)
    colnames(pred_conc) <- colnames(new_data)
    
    #new_data_long$pred_conc_log <- predict_x(par$params, new_data_long$OD)
    #new_data_long$pred_conc <- 10^new_data_long$pred_conc_log
    
    new_data_long$pred_conc <- predicted_x
    new_data_long$pred_conc_log <- log10(new_data_long$pred_conc)
    
    results[['model']] <- np1
    results[['pred_long']] <- new_data_long
    results[['pred']] <- pred_conc
    results[['new_data']] <- new_data
    return(results)
  })
  
  output$no_data_message <- renderText({
    if (is.null(input$input_file)) {
      return("upload file with eliza data to get started")
    }
  })
  
  output$new_values_header <- renderText({
    if (!is.null(input$input_file)) {
      return("new OD values")
    }
  })
  
  output$new_preds_header <- renderText({
    if (!is.null(input$input_file)) {
      return("Predicted concentration based on new OD values")
    }
  })
  
  
  output$curve_plot <- renderPlotly ({
    req(input$input_file)
    results <- fit_model_and_predict()
    print(names(results))
    x1 <- getX(results$model); y1 <- getY(results$model)
    x2 <- getXcurve(results$model); y2 <- getYcurve(results$model)
    
    p <- ggplot() + theme_bw() + geom_point(data=data.frame(logconc=x1, OD=y1, conc=10^x1), aes(x=logconc, y=OD, color="standard_curve", label=conc), size=1) + 
      geom_line(data=data.frame(x=x2, y=y2), aes(x=x, y=y), color="blue") + #coord_cartesian(ylim=c(min(y1),1.1), xlim=c(min(x1)-0.1,max(x1+0.1)), expand=F) +
      xlab("log10 concentration (ug/ml)") + ylab("OD") +# ggtitle("elisa curve") +
      geom_point(data=results$pred_long,aes(x=pred_conc_log,y=OD, color="predicted", label=pred_conc), size=1) +
      scale_color_manual(name="legend",values=c("standard_curve"="black","predicted"="red"))#, size=6, shape=13) +
    
    ggplotly(p, tooltip = c("label","x","y"))
  })

  output$pred_conc_table <- renderDataTable({ 
    req(input$input_file)
    results <- fit_model_and_predict()
    datatable(round(results$pred, digits=5), extensions = c('Buttons'), options=list(dom = 'Bfrtip', buttons=c('copy','csv','print')))

  })
  
  output$od_table <- renderDataTable({ 
    req(input$input_file)
    results <- fit_model_and_predict()
    datatable(results$new_data, extensions = c('Buttons'), options=list(dom = 'Bfrtip', buttons=c('copy','csv','print')))
  
  })
  
  output$full_table <- renderDataTable({ 
    req(input$input_file)
    results <- fit_model_and_predict()
    datatable(results$pred_long, extensions=c('Buttons', 'Scroller'), 
              options=list(dom = 'Bfrtip', buttons=c('copy','csv','print'), scrollY=800, scroller=TRUE))
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

