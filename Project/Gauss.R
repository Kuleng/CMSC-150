source('~/Desktop/CMSC150/Week4/tumpalan_exer03_v2.R')

# #Handout test case
# E1 <- function(x0,x1,x2) 1 * x0 + 0 * x1 + 1 * x2 + -4;
# E2 <- function(x0,x1,x2) 8 * x0 + -3 * x1 + 4 * x2 + -12;
# E3 <- function(x0,x1,x2) 6 * x0 + -3 * x1 + 1 * x2 + -8;
# 
# #Problem 1 - Exer
# D1 <- function(x0,x1,x2) 1 * x0 + 1 * x1 + 1 * x2 + -100;
# D2 <- function(x0,x1,x2) 0.5 * x0 + -0.334 * x2 + 0;
# D3 <- function(x0,x1,x2) -0.5 * x0 + 0.2 * x1 + 0;
# D4 <- function(x0,x1,x2) -.20 * x1 + 0.334 * x2 + 0;
# 
# #Problem 2 - Exer
# C1 <- function(x0,x1) 0.06 * x0 + 0.10 * x1 + -3520
# C2 <- function(x0,x1) 1 * x0 + -2 * x1
# 
# #Problem 3 - Exer
# F1 <- function(x0,x1,x2,x3,x4,x5) .7 * x0 + .2 * x1 + .1 * x2 + .1 * x3 + .05 * x4 + 0 * x5 + -1000
# F2 <- function(x0,x1,x2,x3,x4,x5) .1 * x0 + .6 * x1 + 0 * x2 + .1 * x3 + .05 * x4 + 0.1 * x5 + -700
# F3 <- function(x0,x1,x2,x3,x4,x5) .05 * x0 + .1 * x1 + .75 * x2 + .1 * x3 + 0 * x4 + 0.05 * x5 + -1300
# F4 <- function(x0,x1,x2,x3,x4,x5) .05 * x0 + 0 * x1 + .05 * x2 + .6 * x3 + .05 * x4 + 0.05 * x5 + -900
# F5 <- function(x0,x1,x2,x3,x4,x5) 0 * x0 + .05 * x1 + .05 * x2 + .1 * x3 + .75 * x4 + 0.2 * x5 + -1500
# F6 <- function(x0,x1,x2,x3,x4,x5) .1 * x0 + .05 * x1 + .05 * x2 + 0 * x3 + .1 * x4 + 0.6 * x5 + -1200
# 
# #newList = list(E1, E2, E3); #HandOut test case
# #newList = list(D1,D2,D3,D4); #problem 1
# #newList = list(C1,C2) #problem 2
# newList = list(F1,F2,F3,F4,F5,F6); #problem 3

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

#function call
# 
# GaussJordan(newList)
# 
# cat("\n\n\n\n\n")
# 
# Gaussian(newList)