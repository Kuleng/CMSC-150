source('~/Desktop/CMSC150/Week4/tumpalan_exer03_v2.R')

Gaussian<-function(mainMatrix){

  unknowns = length(mainMatrix[1,])-1 #number of unknown variables
  
  variableValues = rep(0, (unknowns)) #creates a vector of zeroes corresponding to variables
  
  for(pivotIteration in 1:(unknowns-1)){ #iterates all values below the main diagonal
    max = 0
    pivotRow = 0
    PivotElement = 0;

    for (i in pivotIteration:(unknowns)){#loop from the value below the pivot element to the number of equations.
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
    
    #pivot element in main diagonal of the matrix
    PivotElement = mainMatrix[pivotIteration,pivotIteration]

    #if division by zero, prompt    
    if(PivotElement == 0){
      print("Division by zero error. No solution.")
      return(variableValues)
    }
    
    
    for (upperTriangle in pivotIteration:(unknowns-1)){ #looping all values to be eliminated
      #Value to be eliminated
      VTBE = mainMatrix[upperTriangle+1,pivotIteration]
      #Multiplier
      multiplier = VTBE/PivotElement
      
      #Variable for the temporary pivot rows
      tempVector = c(mainMatrix[pivotIteration,])
      tempVector = tempVector*multiplier #Scalar multiplication to rows

      #row - vector(M x PivotRow) loop
      for(element in 1:(unknowns+1)){
        mainMatrix[upperTriangle+1, element] = mainMatrix[upperTriangle+1, element] - (tempVector[element])
      }
    }
  }

  for(var in (unknowns):1){     #for every unknown
    for(coeffs in (unknowns+1):1){ #for every column
      if(coeffs != var){            #if the coefficient is not in the main diagonal, add.
        if(coeffs != (unknowns+1)){ #if the value to be added is not in the right hand side, multiply to its coefficient before adding.
          variableValues[var] = variableValues[var] + (-1)*mainMatrix[var,coeffs]*variableValues[coeffs]
        }
        else{               #if the value to be added is in the right hand sige, add only
          variableValues[var] = variableValues[var] + mainMatrix[var,coeffs]
        }
      }
    }
    #divide both sides by the coefficient of the unknown variable
    variableValues[var] = variableValues[var]/mainMatrix[var,var]
  }
  
  #create new list for results
  results = list(augcoeffmatrix = mainMatrix, values = variableValues)
  
  return(results$values)
}

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
