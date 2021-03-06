library(shiny)
library(shinyMatrix)

ui <- fluidPage(
  
  #changes specific error validation font color to green
  tags$head(
    tags$style(HTML("
      .shiny-output-error-validation {
        color: green;
      }
    "))
  ),
  
  #Main layout of the GUI
  pageWithSidebar(
    titlePanel("CMSC 150 - Project"),
    
    #sidebar panel with multiple conditional panels depending on the generic solver to be selected by the user
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
                       numericInput(inputId="realX", label = "Input a real number", value =1)),
      conditionalPanel(condition = "input.Method == 'Spline Interpolation'",
                       fileInput(inputId="csvfile2", label = "Attach CSV file",
                       multiple = FALSE, accept = c("text/csv", "text/comma-separated-values, text/plain",
                                                    ".csv")),
                       numericInput(inputId="realX2", label = "Input a real number", value=1)),
      conditionalPanel(condition = "input.Method == 'Simplex Method'",
                       selectInput(inputId = "simplexSolution", "Display tableau",
                                   choices = c("Show tableau",
                                               "Hide tableau"), selected = "Hide tableau",
                       ),
                       #INPUT CODE HERE
                       h3("Minimum cost: ")
                       )
      ),
    
    #mainpanel of the layout with conditional panels that changes depending on the generic solver chosen by the user
    mainPanel(conditionalPanel(condition = "input.Method == 'Polynomial Regression'",
                               h3("Polynomial Function: "),
                               textOutput("polynomial_function"),
                               h3("Estimate of f(x): "),
                               textOutput("estimate"),
                               h3("Plot of function"), 
                               plotOutput("trend")),
              conditionalPanel(condition = "input.Method == 'Spline Interpolation'",
                               #SPLINE INTERPOLATION
                               h3("Interval functions: "),
                               textOutput("intervalFunctions"),
                               h3("Estimate of f(x): "),
                               textOutput("intervalEstimate")
                               #plotOutput if needed
                               ),
              conditionalPanel(condition = "input.Method == 'Simplex Method'",
                               #SIMPLEX METHOD
                               matrixInput(inputId = "inputMatrix1", value = matrix(0L, nrow = 1, ncol = 5, dimnames = list(c("Demands"),c("Sacramento","Saltlake", "Albuquerque", "Chicago", "New York"))),
                                           rows = list(names = TRUE, editableNames = FALSE),
                                           cols = list(names = TRUE, editableNames = FALSE),
                                           ),
                               matrixInput(inputId = "inputMatrix2", value = matrix(0L, nrow = 3, ncol = 6, dimnames = list(c("Denver", "Phoenix", "Dallas"),c("Supply", "Shipping costs (Sacramento)","Shipping costs (Saltlake)","Shipping costs (Albuquerque)","Shipping costs (Chicago)", "Shipping costs (New York)"))),
                                           rows = list(names = TRUE, editableNames = FALSE),
                                           cols = list(names = TRUE, editableNames = FALSE),
                                            ),
                               conditionalPanel(condition = "input.simplexSolution == 'Show tableau'",
                                                tableOutput("finalTableau"),
                               )
                            ),
              )
  )
)

server<-function(input, output){
  
  setwd('~/Desktop/CMSC150/Project')

  #outputs text and validates the csvfile 
  output$polynomial_function = renderText({
                    validate(need(input$csvfile != "", "Please select a csv file."))
                    mainFile <- input$csvfile
                    validate(need(input$order > 0, "Please input a valid polynomial order degree"))
                    PolynomialRegression(mainFile$datapath, input$order, input$realX, 1) })
  
  #outputs the estimated value of the real number in the polynomial function
  output$estimate = renderText({                    
                    validate(need(input$csvfile != "", "Please select a csv file."))
                    validate(need(input$order > 0, "Please input a valid polynomial order degree"))
                    mainFile <- input$csvfile
                    PolynomialRegression(mainFile$datapath, input$order, input$realX, 2) })
#                    eval(parse(text= output$polynomial_function )) })
  
  #outputs optional plot of the function
  output$trend = renderPlot({ 
                validate(need(input$csvfile != "", "Please select a csv file."))
                validate(need(input$order > 0, "Please input a valid polynomial order degree"))
                mainFile <- input$csvfile
                PolynomialRegression(mainFile$datapath, input$order, input$realX, 3) } )
  
  #outputs interval function of the spline interpolation intervals
  output$intervalFunctions = renderUI({
    mainFile <-input$csvfile2
    validate(need(input$csvfile2, "Please select a csv file."))
    validate(need(!is.null(QuadSpline(mainFile$datapath, input$realX2, 2)), "Please select a real number inside the bounds of the spline functions"))
    QuadSpline(mainFile$datapath, input$realX2, 1)
  })
  
  #outputs the estimated value of the real number according to the function of the assigned plot
  output$intervalEstimate = renderText({
    mainFile <-input$csvfile2
    validate(need(input$csvfile2, "Please select a csv file."))
    validate(need(!is.null(QuadSpline(mainFile$datapath, input$realX2, 2)), "Please select a real number inside the bounds of the spline functions"))
    QuadSpline(mainFile$datapath, input$realX2, 2)
  })
  
  output$finalTableau = renderTable({
    demand <- input$inputMatrix1
    cost_supply <-input$inputMatrix2
    demandVector <- c(as.numeric(demand))
    costVector <- matrix(as.numeric(cost_supply), nrow = 3)
    mainVector <-c(costVector[,1], demandVector)
    print(costVector)
    print(mainVector)
    simplexSetup(mainVector,costVector)
  })
}
#GaussJordan Function for List parameters
GaussJordan<-function(mainMatrix){
  
  unknowns = length(mainMatrix[1,]) -1
  variableValues = rep(0, (unknowns)) #creates a vector of zeroes corresponding to variables
  
  for(pivotIteration in 1:(unknowns)){ #loops through the rows of the system
    max = 0
    pivotRow = 0
    PivotElement = 0;
    
    for (i in pivotIteration:(unknowns)){#loop for the number of equations
      #get the maximum value on the ith column
      if (max <= abs(mainMatrix[i,pivotIteration])){
        max = abs(mainMatrix[i,pivotIteration])
        pivotRow = i
      }
    }
    
    #swap first row and the new pivot row with max
    if (pivotRow != pivotIteration){
      temp = mainMatrix[pivotIteration,];
      mainMatrix[pivotIteration,] = mainMatrix[pivotRow,]
      mainMatrix[pivotRow,] = temp
    }
    
    #Pivot element before normalization
    PivotElement = mainMatrix[pivotIteration,pivotIteration] #pivot element in main diagonal
    
    if(PivotElement == 0){
      print("Division by zero error. No solution.")
      return(variableValues)
    }
    
    mainMatrix[pivotIteration,] = mainMatrix[pivotIteration,]/PivotElement
    
    #pivot element after normalization
    PivotElement = mainMatrix[pivotIteration,pivotIteration]
    for (upperTriangle in 1:(unknowns)){ #looping all values to be eliminated to create an identity matrix
      if(upperTriangle != pivotIteration){ #ignore the values in the main diagonal and RHS
        VTBE = mainMatrix[upperTriangle,pivotIteration]
        multiplier = VTBE/PivotElement
        
        #temporary storage of pivot row
        tempVector = c(mainMatrix[pivotIteration,])
        tempVector = tempVector*multiplier #M x PivotRow
        
        #row - vector(MxPivotRow) loop
        for(element in 1:(unknowns+1)){
          mainMatrix[upperTriangle, element] = mainMatrix[upperTriangle, element] - (tempVector[element])
        }
      }
    }
  }
  
  for(i in 1:(unknowns)) #input results to a vector
    variableValues[i] = mainMatrix[i,(unknowns+1)]
  
  #create new list for results
  results = list(augcoeffmatrix = mainMatrix, values = variableValues)
  
  return(results$values)
}

QuadSpline<-function(csvfile, realX, returnType){
  
  mydata = read.csv(file = csvfile)
  
  #data points
  x = mydata[,1]
  y = mydata[,2]
  
  #intervals and knots
  intervals = length(x)-1
  interiorKnots = (intervals-1)*2
  endPoints = 2
  equivalenceEqs = intervals-1
  RHS = intervals*3
  
  
  RHSvalues = c()
  
  #equations
  mainMatrix = matrix(0L, nrow = (intervals*3)-1,ncol =(intervals*3)-1 , dimnames = list(c(),c()))
  
  #CONDITION 1 FILL UP MATRIX
  xi_input = 1
  colInput = 1
  for(i in 1:(interiorKnots)){
    for(j in 1:3){
      if(i == 1){
        mainMatrix[i,colInput] = x[i+1]
        mainMatrix[i, colInput+1] = 1
        colInput = colInput + 2
        break;
      }
      else{
        if(j == 1)
          mainMatrix[i,colInput] = x[xi_input+1]^2
        else if (j==2)
          mainMatrix[i,colInput] = x[xi_input+1]
        else if (j==3)
          mainMatrix[i,colInput] = 1
        colInput = colInput + 1
      }
    }
    RHSvalues = c(RHSvalues, y[xi_input+1])
    if(i%%2 == 0){
      colInput = colInput-3
      xi_input = xi_input + 1
    }
  }
  
  
  #CONDITION 2
  xi_input = 1
  colInput = 1
  for(rowInput in (interiorKnots+1):(interiorKnots+2)){
    for(j in 1:3){
      if(rowInput == (interiorKnots+1)){
        mainMatrix[rowInput,j] = x[xi_input]
        mainMatrix[rowInput,j+1] = 1
        colInput = (colInput+2)+(3*((interiorKnots/2)-1))
        xi_input = length(x)
        break;
      }
      else{
        mainMatrix[rowInput, colInput] = x[xi_input]^(3-j)
      }
      colInput = colInput + 1
    }
    rowInput = rowInput + 1
  }
  RHSvalues = c(RHSvalues, y[1])
  RHSvalues = c(RHSvalues, y[xi_input])
  
  #CONDITION 3
  rowInput = 1
  colInput = 1
  xi_input = 2
  for(rowInput in ((intervals*2)+1):(RHS-1)){
    for(j in 1:3){
      if(j == 1){
        if(rowInput != ((intervals*2)+1)){
          mainMatrix[rowInput,colInput] = x[xi_input]*2
          colInput = colInput +1
        }
      }
      else if(j==2){ 
        mainMatrix[rowInput,colInput] = 1
        colInput = colInput + 1
      }
      else
        colInput = colInput + 1
    }
    
    for(j in 1:3){
      if(j == 1){
        mainMatrix[rowInput,colInput] = x[xi_input]*2*(-1)
        colInput = colInput + 1
      }
      else if(j==2){ 
        mainMatrix[rowInput,colInput] = (-1)
        colInput = colInput + 1
      }
    }
    RHSvalues = c(RHSvalues,0)
    colInput = colInput - 2
    rowInput = rowInput + 1
    xi_input = xi_input + 1
  }
  
  #print(RHSvalues)
  #print(mainMatrix)
  mainMatrix = cbind(mainMatrix,RHSvalues)
  unknownValues = GaussJordan(mainMatrix)
  #unknownValues = MatrixGaussJordan(mainMatrix, RHSvalues)
  
  #print(unknownValues)
  
  functionList <- c()
  functionsAdded = 1
  for(func in 1:(intervals)){
    functionString = "function(x) "
    if(func == 1){
      functionString = paste(functionString, unknownValues[1],"*x +" ,unknownValues[2], sep="", collapse = NULL)
      functionList[functionsAdded]<-functionString
      functionsAdded = functionsAdded + 1
      unknownsUsed = 3
    }
    else{
      functionString = "function(x) "
      functionString = paste(functionString, unknownValues[unknownsUsed],"*x^2 +", sep="", collapse = NULL)
      unknownsUsed = unknownsUsed + 1
      functionString = paste(functionString, unknownValues[unknownsUsed],"*x +", sep="", collapse = NULL)
      unknownsUsed = unknownsUsed + 1
      functionString = paste(functionString, unknownValues[unknownsUsed], sep="", collapse = NULL)
      unknownsUsed = unknownsUsed + 1
      functionList[functionsAdded]<- functionString
      #functionList[functionsAdded]<- eval(parse(text = functionString))
      functionsAdded = functionsAdded + 1
    }
  }
  #print(functionList)
  
  interval = NULL

  #print(x)
  for(value in 1:(length(x))){
    if(realX >= x[value] && realX <= x[value+1]){
      #interval 1
      interval = value
    }
  }
  
  #realX value in interval
  #print(functionList[[interval]](realX))
  #print(solve(mainMatrix, RHSvalues))
  
  if(returnType == 1){
    if(is.null(interval))
      return (NULL)
    return(functionList)
  } else if(returnType == 2){
    if(is.null(interval))
      return (NULL)
    function_fx<-eval(parse(text = functionList[interval]))
    return(function_fx(realX))
  }
}

GaussMod<-function(mainMatrix){
  
  highestNeg = 0
  lowestPos = 999999
  isNeg = 0
  #check the number of negative numbers in the objective function
  for(num in mainMatrix[length(mainMatrix[,1]),]){
    if(num < 0)
      isNeg = isNeg + 1
  }
  
  while(isNeg > 0){
    isNeg <- 0
    
    #find maximum pivot column
    for(PC in 1:(length(mainMatrix[1,])-1)){
      #check bottom side for highest negative integer
      if( mainMatrix[length(mainMatrix[,1]), PC] < highestNeg){
        highestNeg = mainMatrix[length(mainMatrix[,1]), PC]
        maxPC = PC
      }
    }
    
    #find minimum pivot row
    for(PR in 1:(length(mainMatrix[,1])-1)){
      if(mainMatrix[PR, maxPC] != 0 && mainMatrix[PR, length(mainMatrix[1,])]/mainMatrix[PR, maxPC] > 0){
        if( mainMatrix[PR, length(mainMatrix[1,])]/mainMatrix[PR, maxPC] < lowestPos ){
          lowestPos = mainMatrix[PR, length(mainMatrix[1,])]/mainMatrix[PR, maxPC]
          minPR = PR
        }
      }
    }
    
    
    print(maxPC)
    print(minPR)    
    #normalize
    mainMatrix[minPR,] =mainMatrix[minPR,]/mainMatrix[minPR, maxPC];
    
    #choose pivot element    
    PE = mainMatrix[minPR, maxPC]
    
    cat("Pivot element :")
    cat(PE)
    cat("\n")
    
    cat("Result of normalization:\n")
    print(mainMatrix)
    cat("\n")
    
    
    for(VTBE in 1:(length(mainMatrix[,PC]))){
      if(VTBE != minPR && mainMatrix[VTBE,maxPC] != 0){
        
        #get multiplier by value to be eliminated over PE
        multiplier = mainMatrix[VTBE, maxPC]/PE
        cat("Multiplier")
        cat(multiplier)
        cat("\n")
        
        vectorToSubtract = multiplier*mainMatrix[minPR,]
        
        mainMatrix[VTBE,] = mainMatrix[VTBE,] - vectorToSubtract
      }
    }
    
    print(mainMatrix)
    cat("\n")
    cat("-------")
    cat("\n")
    highestNeg = 0
    lowestPos = 999999
    for(num in mainMatrix[length(mainMatrix[,1]),]){
      if(num < 0)
        isNeg = isNeg + 1
    }
    
    cat("isNeg")
    cat(isNeg)
    cat("\n")
  }
  
  return(mainMatrix)
}

PolynomialRegression<-function(csvfile, degree, realX, int){
  
  mydata = read.csv(file = csvfile)
  x = mydata[,1]
  y = mydata[,2]
  
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
  plot_f <- plot(mainList[[1]] , mainList[[2]], pch = 20, col = "red", main = "Function plot", xlab = "X", ylab = "Y")
  
  #set lines
  coordinates = lines(x,predict(qModel), col="blue")
  
  #manipulate polynomial String from right to left, if already in the left most, do not add any string manipulated variables and operations
  
  for(degCount in (degree+1):1){
    if(degCount == 1)
      mainData$polynomialString = paste(mainData$polynomialString, mainData$unknowns[[degCount]])
    else{
      mainData$polynomialString = paste(mainData$polynomialString, mainData$unknowns[[degCount]], " * x ^ ", (degCount-1) ," + ", sep = "", collapse = NULL)
    }
  }
  
  #create the function through evaluation and parsing functions
  #print(mainData$polynomialString)
  mainData$polynomialFunction = eval(parse(text = mainData$polynomialString))
  #print(mainData$polynomialFunction)
  print(mainData)

    if(int == 1)
    return(mainData$polynomialString)
  else if(int == 2)
    return(mainData$polynomialFunction(realX))
  else
    return(plot_f)
}

simplexSetup<-function(demand_supply_vector, costMatrix){
  
  #Sets up the initial matrix based from the initial positions of the constraints and equations
  mainMatrix = matrix(0L, nrow = 9, ncol = 16)
  
  #sets up the supply constraints
  mainMatrix[1,] = c(-1,-1,-1,-1,-1, 0,0,0,0,0,0,0,0,0,0, demand_supply_vector[1]*(1))
  mainMatrix[2,] = c(0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0, demand_supply_vector[2]*(1))
  mainMatrix[3,] = c(0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1, demand_supply_vector[3]*(1))
  
  
  #Sets up the demand constraints to the matrix
  col = 4
  i = 1
  demandInput = 4
  while(col < 9){
    while(i <= 15){
      mainMatrix[col, i] = 1
      i = i + 5
    }
    mainMatrix[col, 16] = demand_supply_vector[demandInput]*(-1)
    demandInput = demandInput + 1
    col = col + 1
    i = col - 3
  }
  
  #mainMatrix[9,] = c(5,6,7,8,9,6,7,8,9,10,3,5,7,11,13, 0)
  
  #Gets the initial cost values per shipping area to the matrix
  newMatrix <- matrix(costMatrix[,-1], nrow = 3, byrow = FALSE)
  costVector <- c()
  
  input = 1
  for(i in 1:(length(newMatrix[,1]))){
    for(j in 1:(length(newMatrix[1,]))){
      costVector[input] = newMatrix[i, j]
      input = input + 1
    }
  }
  
  #print(costVector)
  mainMatrix[9,] = c(costVector, 0L)
  
  print(mainMatrix)
  
  #Tranposes the initial matrix
  transposedMatrix = t(mainMatrix)
  
  #creates a matrix of slack variables for the identification of unknowns
  unknownsMatrix = matrix(0L, nrow = 16, ncol = 16)
  for(i in 1:(length(unknownsMatrix[,1]))){
    for(j in 1:(length(unknownsMatrix[1,]))){
      if(i == j)
        if(i > 3)
          unknownsMatrix[i,j] = 1
        else
          unknownsMatrix[i,j] = 1
    }
  }
  
  #binds the matrix of slack variables to the initial matrix
  transposedMatrix = cbind(transposedMatrix, unknownsMatrix)
  
  #Sort the merged matrix to adjust RHS to the rightmost side of the matrix
  for(int in 9:(length(transposedMatrix[1,])-1)){
    temp = transposedMatrix[,int]
    transposedMatrix[,int] = transposedMatrix[,int+1]
    transposedMatrix[,int+1] = temp
  }
  
  #Perform Modified Gauss Jordan for Simplex Approach
  GaussMod(transposedMatrix)
}


shinyApp(ui= ui, server= server)