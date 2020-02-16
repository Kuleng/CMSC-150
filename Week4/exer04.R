source('~/Desktop/CMSC150/Week4/tumpalan_exer03_v2.R')

E1 <- function(x0,x1,x2) 1 * x0 + -0 * x1 + 1 * x2 + -4;
E2 <- function(x0,x1,x2) 8 * x0 + 4 * x2 + -3 * x1 + -12; #no constant
E3 <- function(x0,x1,x2) 6 * x0 + -3 * x1 + 1 * x2 + -8; #missing 1 variable

newList = list(E1, E2, E3);
numOfEq = length(newList); #number of equations
newMatrix = AugCoeffMatrix(newList) #Augmented Coefficient Matrix

Gaussian<-function(newList){
  
  variables = c(newMatrix$variables) #creates a vector of the variables
  variableValues = rep(0, 3) #creates a vector of zeroes corresponding to variables
 
  print(newMatrix$augcoeffmatrix)
  
  cat("Forward Elimination: ")
  
  for(pivotIteration in 1:(length(newList))){
    max = 0
    pivotRow = 0
    PivotElement = 0;

    for (i in pivotIteration:numOfEq){#loop for the number of equations.
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
    
    cat("\n\nResult of pivot: \n")
    print(newMatrix$augcoeffmatrix)
    
    #pivot element in main diagonal of the matrix
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration]
    
    
    for (upperTriangle in pivotIteration:(length(newList))){ #looping all values to be eliminated
      
      #Value to be eliminated
      VTBE = newMatrix$augcoeffmatrix[upperTriangle+1,pivotIteration]
      #Multiplier
      multiplier = VTBE/PivotElement
      
      #Variable for the temporary pivot rows
      tempVector = c(newMatrix$augcoeffmatrix[pivotIteration,])
      tempVector = tempVector*multiplier #Scalar multiplication to rows

      
      #row - vector(M x PivotRow) loop
      for(element in 1:(numOfEq+1)){
        newMatrix$augcoeffmatrix[upperTriangle+1, element] = newMatrix$augcoeffmatrix[upperTriangle+1, element] - (tempVector[element])
      }
      
      cat("\n\nCurrent pivot row:  ")
      cat(newMatrix$augcoeffmatrix[pivotIteration,])
      cat("\n")  
      cat("Pivot element: ")
      cat(PivotElement)
      cat("\n")
      cat("Value to be eliminated:  ")
      cat(VTBE)
      cat("\n")
      cat("Vector:  ")
      cat(tempVector)
      cat("\n")
      
      cat("\nResulting matrix:  \n")
      print(newMatrix$augcoeffmatrix)
      cat("\n")
    }
  }
  
  cat("Backward substitution:  \n")
  
  for(var in (numOfEq):1){#for every unknown
    for(coeffs in (numOfEq+1):1){ #for every column
      if(coeffs != var){ #if the coefficient is not in the main diagonal, add.
        if(coeffs != (numOfEq+1)){ #if the value to be added is not in the right hand side, multiply to its coefficient.
          variableValues[var] = variableValues[var] + (-1)*newMatrix$augcoeffmatrix[var,coeffs]*variableValues[coeffs]
        }
        else{ #if the value to be added is in the right hand sige, add only
          variableValues[var] = variableValues[var] + newMatrix$augcoeffmatrix[var,coeffs]
        }
      }
    }
    
    for(a in 1:(numOfEq+1)){ #Main loop for display
      if(newMatrix$augcoeffmatrix[var,a] != 0){
        cat(newMatrix$augcoeffmatrix[var,a]) #Displays the expression if it doesn't have 0 coeff
        if(!is.na(variables[a])){
          cat(" * ") #if it has a variable, display asterisk (multiplication)
          cat(variables[a])
        }
        if (a == numOfEq) #if the variable is the last variable, display equality
          cat(" = ")
        else if(a < numOfEq) # else, display +
          cat(" + ")
      }
    }
    
    cat("\n")
    cat(newMatrix$augcoeffmatrix[var,var])
    cat(" * ")
    cat(variables[var])
    cat(" = ")
    cat(variableValues[var])
    cat("\n")
    #divide both sides by the coefficient of the unknown variable
    variableValues[var] = variableValues[var]/newMatrix$augcoeffmatrix[var,var]
    cat(variables[var])
    cat(" = ")
    cat(variableValues[var])
    cat("\n\n")
  }
  
  #output results
  print(variableValues)
}

GaussJordan<-function(newList){
  
  variables = c(newMatrix$variables) #creates a vector of the variables
  variableValues = rep(0, 3) #creates a vector of zeroes corresponding to variables
  
  print(newMatrix$augcoeffmatrix)
  
  for(pivotIteration in 1:(numOfEq)){ #loops through the rows of the system
    max = 0
    pivotRow = 0
    PivotElement = 0;
  
    for (i in pivotIteration:numOfEq){#loop for the number of equations
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
    
    cat("\n\nResult of pivot: \n")
    print(newMatrix$augcoeffmatrix)
    
    #Pivot element before normalization
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration] #pivot element in main diagonal
    
    newMatrix$augcoeffmatrix[pivotIteration,] = newMatrix$augcoeffmatrix[pivotIteration,]/PivotElement
    
    #pivot element after normalization
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration]
    
    cat("\n\nResult of Normalization: \n")
    print(newMatrix$augcoeffmatrix)
    
    for (upperTriangle in 1:(numOfEq)){ #looping all values to be eliminated
      if(upperTriangle != pivotIteration){ #ignore the values in the main diagonal and RHS
        VTBE = newMatrix$augcoeffmatrix[upperTriangle,pivotIteration]
        multiplier = VTBE/PivotElement
    
        #temporary storage of pivot row
        tempVector = c(newMatrix$augcoeffmatrix[pivotIteration,])
        tempVector = tempVector*multiplier #M x PivotRow
  
        #row - vector(MxPivotRow) loop
        for(element in 1:(numOfEq+1)){
          newMatrix$augcoeffmatrix[upperTriangle, element] = newMatrix$augcoeffmatrix[upperTriangle, element] - (tempVector[element])
        }
        
        cat("\n\nCurrent pivot row:  ")
        cat(newMatrix$augcoeffmatrix[pivotIteration,])
        cat("\n")  
        cat("Pivot element: ")
        cat(PivotElement)
        cat("\n")
        cat("Value to be eliminated:  ")
        cat(VTBE)
        cat("\n")
        cat("Vector:  ")
        cat(tempVector)
        cat("\n")
        
        cat("\nResulting matrix:  \n")
        print(newMatrix$augcoeffmatrix)
        cat("\n")
      }
    }
  }
  
  for(i in 1:(numOfEq)) #input results to a vector
    variableValues[i] = newMatrix$augcoeffmatrix[i,(numOfEq+1)]
  
  #output results
  print(variableValues)
}

#function call
GaussJordan(newMatrix)

Gaussian(newMatrix)