library(shiny)
library(shinyMatrix)

ui <- fluidPage(
  pageWithSidebar(
      titlePanel("CMSC 150 - Project"),
      sidebarPanel(
        selectInput(inputId = "Method", "Please select your choice.",
                    choices = c("Polynomial Regression",
                                "Spline Interpolation",
                                "Simplex Method")
                    ),
        conditionalPanel(condition = "input.Method == 'Polynomial Regression'",
                         fileInput(inputId="csvfile", label = "Attach CSV file",
                                   multiple = FALSE,
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv")),
                         numericInput(inputId="order", label = "Input polynomial order", value = 3),
                         numericInput(inputId="realX", label = "Input a real number", value =1))),
      mainPanel(conditionalPanel(condition = "input.Method == 'Polynomial Regression'",
                                h2("Test data points table"),
                                DT::dataTableOutput("dataPoints")))
  )
)

server<-function(input, output){
  
  setwd('~/Desktop/CMSC150/Project')
  
  source('Gauss.R')
  mydata <- reactive({
    
    inFile <- input$csvfile
    
    if (is.null(inFile))
      return(NULL)
    
    mainTable <- read.csv(inFile$datapath, header=input$header, sep=input$sep,  dec = input$dec)
    
    return(mainTable)
  })

  degree = input$order

  x = mainTable[,1]
  y = mainTable[,2]

  mainList = list(x,y)
  
  #create the main labeled list
  mainData = list(augcoeffmatrix = matrix(0L, byrow = TRUE, nrow = degree+1, ncol = degree+2), unknowns = c(), polynomialString = "function(x) ", polynomialFunction = "function(x) ")
  
  
  #Set up aug coeff matrix
  for(row in 1:(degree+1)){
    for(col in 1:(degree + 2)){
      if (col != (degree + 2))
        mainData$augcoeffmatrix[row,col] = sum(mainList[[1]]^((row-1) + (col-1))) #if the coeff is not in the right hand side
      else
        mainData$augcoeffmatrix[row,col] = sum(mainList[[2]] * mainList[[1]]^(row-1)) #if the coeff is in the RHS
    }
  }
  
  #get the unknowns
  mainData$unknowns = GaussJordan(mainData$augcoeffmatrix)
  
  #getting the qModel (polynomial Model)
  if(degree == 1)
    qModel = lm(y~x)
  else
    qModel = lm(y~poly(x,degree,raw = T))
  
  #plot the function
  plot(mainList[[1]] , mainList[[2]], pch = 20, col = "red", main = "Function plot", xlab = "X", ylab = "Y")
  
  #set lines
  coordinates = lines(x,predict(qModel), col="blue")
  
  #manipulate polynomial String from right to left, if already in the left most, do not add any string manipulated variables and operations
  
  for(degCount in (degree+1):1){
    if(degCount == 1)
      mainData$polynomialString = paste(mainData$polynomialString, qModel$coefficients[[degCount]])
    else{
      mainData$polynomialString = paste(mainData$polynomialString, qModel$coefficients[[degCount]], " * x ^ ", (degCount-1) ," + ", sep = "", collapse = NULL)
    }
  }
  
  #create the function through evaluation and parsing functions
  mainData$polynomialFunction = eval(parse(text = mainData$polynomialString))
  
  
  output$dataPoints = DT::renderDataTable({ mainData$augcoeffmatrix })
  output$trend = renderPlot({ plot(mainList[[1]] , mainList[[2]], pch = 20, col = "red", main = "Function plot", xlab = "X", ylab = "Y") } )
}

shinyApp(ui= ui, server= server)