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
                       numericInput(inputId="realX", label = "Input a real number", value =1)),
      conditionalPanel(condition = "input.Method == 'Spline Interpolation'",
                       fileInput(inputId="csvfile2", label = "Attach CSV file",
                                 multiple = FALSE, accept = c("text/csv", "text/comma-separated-values, text/plain",
                                                              ".csv")),
                       numericInput(inputId="realX", label = "Input a real number", value=1)),
      conditionalPanel(condition = "input.Method == 'Simplex Method'",
                       #INPUT CODE HERE
      )
    ),
    mainPanel(conditionalPanel(condition = "input.Method == 'Polynomial Regression'",
                               h2("Polynomial Function: "),
                               textOutput("polynomial_function"),
                               h2("Estimate of f(x): "),
                               textOutput("estimate"),plotOutput("trend")),
              conditionalPanel(condition = "input.Method == 'Spline Interpolation'",
                               #SPLINE INTERPOLATION
              ),
              conditionalPanel(condition = "input.Method == 'Simplex Method'",
                               #SIMPLEX METHOD
              ),
    )
  )
)

server<-function(input, output){
  
  setwd('~/Desktop/CMSC150/Project')
  
  # source('Gauss.R')
  
  output$polynomial_function = renderText({
    mainFile <- input$csvfile
    PolynomialRegression(mainFile$datapath, input$order, input$realX, 1) })
  output$estimate = renderText({                    
    mainFile <- input$csvfile
    PolynomialRegression(mainFile$datapath, input$order, input$realX, 2) })
  #                    eval(parse(text= output$polynomial_function )) })
  output$trend = renderPlot({ 
    mainFile <- input$csvfile
    PolynomialRegression(mainFile$datapath, input$order, input$realX, 3) } )
}

#AugCoeffMatrix
AugCoeffMatrix <-function(mainList){
  #checking if unknowns are equal to the number of
  numOfEq = length(mainList) #number of equations
  unknowns = 0
  checked = 0 #variable for equality check
  loopCount = 0 #variable for equality check (loop)
  varCount = 0
  outerLoopCount = 0

  #Checker for number of unknowns
  for (a in mainList){
    loopCount = 0
    tempDeparse = deparse(a, width.cutoff = 400L); #deparsed function
    tempCoeffs = strsplit(tempDeparse[2],split = "+", fixed = T) #splitted the string that was deparsed to elements
    splitCoeffs = gsub(" ", "", tempCoeffs[[1]], fixed = T)
    splittedCoeffs = strsplit(splitCoeffs, split = "*", fixed = T) #Splitted Coefficient and Variables
    while ( loopCount < length(splittedCoeffs) - 1){
      loopCount = loopCount + 1
      if (as.numeric(substring(splittedCoeffs[[loopCount]][2],2)) > unknowns){
        unknowns = as.numeric(substring(splittedCoeffs[[loopCount]][2],2))
      }
    }
  }

  unknowns = unknowns + 1 #increment for x0 count

  if (unknowns <= numOfEq){
    #  print("Number of equations and unknowns are equal.")
  } else{
    #  print("Number of equations and unknowns are not equal.")
  }

  #create a list of all the unknown variables
  mainVariables = (c())
  a = 0
  while(a<unknowns){
    a = a + 1
    mainVariables[a] = paste("x", a-1, sep="")
  }

  #initialize mainData
  mainData = list(variables = mainVariables, augcoeffmatrix = matrix(0L, byrow = TRUE, nrow = numOfEq, ncol = unknowns+1,dimnames = list(c(1:numOfEq),c(mainVariables,"RHS"))))

  for (a in mainList){
    outerLoopCount = outerLoopCount + 1
    loopCount = 0
    tempDeparse = deparse(a, width.cutoff = 200L); #deparsed function
    tempCoeffs = strsplit(tempDeparse[2],split = "+", fixed = T) #splitted the string that was deparsed to elements
    splitCoeffs = gsub(" ", "", tempCoeffs[[1]], fixed = T)
    splittedCoeffs = strsplit(splitCoeffs, split = "*", fixed = T) #Splitted Coefficient and Variables
    #split all the elements and get their respective values via indexing, negate all the right hand side

    while (loopCount < length(splittedCoeffs)){
      loopCount = loopCount + 1
      for (var in mainVariables){ #loop the coefficients variable to all unknowns that was initialized.
        if(is.na(splittedCoeffs[[loopCount]][2])){ #if there is no variable, it will automatically be set to the right hand side
          mainData$augcoeffmatrix[outerLoopCount,unknowns+1] = (-1)*as.numeric(splittedCoeffs[[loopCount]][1])
        } else if(splittedCoeffs[[loopCount]][2] == var){ #if the variable is equal to the looped variable in the vector of unknowns, it will replace the 0s in the array
          mainData$augcoeffmatrix[outerLoopCount,as.numeric(substring(var,2,2))+1] = as.numeric(splittedCoeffs[[loopCount]][1])
        }
      }
    }
  }

  return (mainData)
}

#GaussJordan Function
GaussJordan<-function(newList){

  numOfEq = length(newList); #number of equations
  newMatrix = AugCoeffMatrix(newList) #Augmented Coefficient Matrix
  unknowns = length(newMatrix$variables) #number of unknown variables

  variables = c(newMatrix$variables) #creates a vector of the variables
  variableValues = rep(0, (unknowns)) #creates a vector of zeroes corresponding to variables

  #print(newMatrix$augcoeffmatrix)

  for(pivotIteration in 1:(unknowns)){ #loops through the rows of the system based on number of unknowns
    max = 0
    pivotRow = 0
    PivotElement = 0;

    for (i in pivotIteration:(numOfEq)){#loop for the number of equations
      #get the maximum value on the ith column
      if (max <= abs(newMatrix$augcoeffmatrix[i,pivotIteration])){
        max = abs(newMatrix$augcoeffmatrix[i,pivotIteration])
        pivotRow = i
      }
    }

    #swap first row and the new pivot row with max
    if (pivotRow != pivotIteration){
      temp = newMatrix$augcoeffmatrix[pivotIteration,];
      newMatrix$augcoeffmatrix[pivotIteration,] = newMatrix$augcoeffmatrix[pivotRow,]
      newMatrix$augcoeffmatrix[pivotRow,] = temp
    }

    #  cat("\n\nResult of pivot: \n")
    #  print(newMatrix$augcoeffmatrix)

    #Pivot element before normalization
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration] #pivot element in main diagonal

    if(PivotElement == 0){
      #   print("Division by zero error. No solution.")
      return(variableValues)
    }

    newMatrix$augcoeffmatrix[pivotIteration,] = newMatrix$augcoeffmatrix[pivotIteration,]/PivotElement

    #pivot element after normalization
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration]

    #  cat("\n\nResult of Normalization: \n")
    # print(newMatrix$augcoeffmatrix)

    for (upperTriangle in 1:(numOfEq)){ #looping all values to be eliminated to create an identity matrix
      if(upperTriangle != pivotIteration){ #ignore the values in the main diagonal and RHS
        VTBE = newMatrix$augcoeffmatrix[upperTriangle,pivotIteration]
        multiplier = VTBE/PivotElement

        #temporary storage of pivot row
        tempVector = c(newMatrix$augcoeffmatrix[pivotIteration,])
        tempVector = tempVector*multiplier #M x PivotRow

        #row - vector(MxPivotRow) loop
        for(element in 1:(unknowns+1)){
          newMatrix$augcoeffmatrix[upperTriangle, element] = newMatrix$augcoeffmatrix[upperTriangle, element] - (tempVector[element])
        }

        #    cat("\n\nCurrent pivot row:  ")
        #   cat(newMatrix$augcoeffmatrix[pivotIteration,])
        #   cat("\nPivot element: ", PivotElement,"\n")
        #   cat("Value to be eliminated:  ", VTBE, "\n")
        #   cat("Vector: \n ")
        #    print(tempVector)
        #    cat("\n")
        #    cat("\nResulting matrix:  \n")
        #    print(newMatrix$augcoeffmatrix)
        #    cat("\n")
      }
    }
  }

  for(i in 1:(unknowns)) #input results to a vector
    variableValues[i] = newMatrix$augcoeffmatrix[i,(unknowns+1)]

  #create new list for results
  results = list(variables = newMatrix$variables, augcoeffmatrix = newMatrix$augcoeffmatrix, values = variableValues)
  #print(variableValues)

  return(results$values)
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
      mainData$polynomialString = paste(mainData$polynomialString, qModel$coefficients[[degCount]])
    else{
      mainData$polynomialString = paste(mainData$polynomialString, qModel$coefficients[[degCount]], " * x ^ ", (degCount-1) ," + ", sep = "", collapse = NULL)
    }
  }
  
  #create the function through evaluation and parsing functions
  #print(mainData$polynomialString)
  mainData$polynomialFunction = eval(parse(text = mainData$polynomialString))
  #print(mainData$polynomialFunction)
  if(int == 1)
    return(mainData$polynomialString)
  else if(int == 2)
    return(mainData$polynomialFunction(realX))
  else
    return(plot_f)
}

shinyApp(ui= ui, server= server)