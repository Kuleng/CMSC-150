#source('~/Desktop/CMSC150/Project/Gauss.R')

source('~/Desktop/CMSC150/Week5/tumpalan_exer04.R')
options(max.print=999999)
MatrixGaussJordan<-function(newMatrix, RHS){

  unknowns = length(newMatrix[,1])
  numOfEq = length(newMatrix[,1])
  
  #Add the right hand side values to the matrix parameter
  newMatrix = cbind(newMatrix,RHS)
  
  #Testprint
  #print(newMatrix)
  
  #variables = c(newMatrix$variables) #creates a vector of the variables
  variableValues = rep(0, (unknowns)) #creates a vector of zeroes corresponding to variables
  
  #print(newMatrix)
  
  for(pivotIteration in 1:(unknowns)){ #loops through the rows of the system based on number of unknowns
    max = 0
    pivotRow = 0
    PivotElement = 0;
    
    for (i in pivotIteration:(numOfEq)){#loop for the number of equations
      #get the maximum value on the ith column
      if (max <= abs(newMatrix[i,pivotIteration])){
        max = abs(newMatrix[i,pivotIteration])
        pivotRow = i
      }
    }
    
    #swap first row and the new pivot row with max
    if (pivotRow != pivotIteration){
      temp = newMatrix[pivotIteration,];
      newMatrix[pivotIteration,] = newMatrix[pivotRow,]
      newMatrix[pivotRow,] = temp
    }
    
    #  cat("\n\nResult of pivot: \n")
    #  print(newMatrix)
    
    #Pivot element before normalization
    PivotElement = newMatrix[pivotIteration,pivotIteration] #pivot element in main diagonal
    
    if(PivotElement == 0){
      #   print("Division by zero error. No solution.")
      return(variableValues)
    }
    
    newMatrix[pivotIteration,] = newMatrix[pivotIteration,]/PivotElement
    
    #pivot element after normalization
    PivotElement = newMatrix[pivotIteration,pivotIteration]
    
    #  cat("\n\nResult of Normalization: \n")
    # print(newMatrix)
    
    for (upperTriangle in 1:(numOfEq)){ #looping all values to be eliminated to create an identity matrix
      if(upperTriangle != pivotIteration){ #ignore the values in the main diagonal and RHS
        VTBE = newMatrix[upperTriangle,pivotIteration]
        multiplier = VTBE/PivotElement
        
        #temporary storage of pivot row
        tempVector = c(newMatrix[pivotIteration,])
        tempVector = tempVector*multiplier #M x PivotRow
        
        #row - vector(MxPivotRow) loop
        for(element in 1:(unknowns+1)){
          newMatrix[upperTriangle, element] = newMatrix[upperTriangle, element] - (tempVector[element])
        }
        
        #    cat("\n\nCurrent pivot row:  ")
        #   cat(newMatrix[pivotIteration,])
        #   cat("\nPivot element: ", PivotElement,"\n")
        #   cat("Value to be eliminated:  ", VTBE, "\n")
        #   cat("Vector: \n ")
        #    print(tempVector)
        #    cat("\n")
        #    cat("\nResulting matrix:  \n")
        #    print(newMatrix)
        #    cat("\n")
      }
    }
  }
  
  for(i in 1:(unknowns)) #input results to a vector
    variableValues[i] = newMatrix[i,(unknowns+1)]
  
  #create new list for results
  #results = list(variables = newMatrix$variables, augcoeffmatrix = newMatrix, values = variableValues)
  #print(variableValues)
  
  return(variableValues)
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
        mainMatrix[rowInput,colInput] =(-1)
        colInput = colInput + 1
      }
    }
    RHSvalues = c(RHSvalues,0)
    colInput = colInput - 2
    rowInput = rowInput + 1
    xi_input = xi_input + 1
  }

  #print(RHSvalues)

  mainMatrix = cbind(mainMatrix,RHSvalues)
  print(mainMatrix)
  #unknownValues = MatrixGaussJordan(mainMatrix, RHSvalues)
  
  unknownValues = GaussJordan(mainMatrix)
  print(unknownValues)
  
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

  #interval = 0
  print(x)
  for(value in 1:(length(x))){
    if(realX > x[value] && realX < x[value+1]){
      #interval 1
      interval = value
    }
  }
  
  print(interval)
  
  #realX value in interval
  #print(functionList[[interval]](realX))
  #print(solve(mainMatrix, RHSvalues))
  
  if(returnType == 1){
    return(functionList)
  } else if(returnType == 2){
    print(functionList)
    print(interval)
    print(functionList[1])
    print(eval(parse(text=functionList[interval])))
    function_fx<-eval(parse(text = functionList[interval]))
    return(function_fx(realX))
  }
}

#print(QuadSpline('~/Desktop/CMSC150/Project/qsidata.csv', 6, 1))
print(QuadSpline('~/Desktop/CMSC150/Project/test.csv', 0, 2))
