source('~/Desktop/CMSC150/Week4/tumpalan_exer03_v2.R')

#Handout test case
E1 <- function(x0,x1,x2) 1 * x0 + 0 * x1 + 1 * x2 + -4;
E2 <- function(x0,x1,x2) 8 * x0 + -3 * x1 + 4 * x2 + -12;
E3 <- function(x0,x1,x2) 6 * x0 + -3 * x1 + 1 * x2 + -8;

#Problem 1 - Exer
D1 <- function(x0,x1,x2) 1 * x0 + 1 * x1 + 1 * x2 + -100;
D2 <- function(x0,x1,x2) 0.5 * x0 + -0.334 * x2 + 0;
D3 <- function(x0,x1,x2) -0.5 * x0 + 0.2 * x1 + 0;
D4 <- function(x0,x1,x2) -.20 * x1 + 0.334 * x2 + 0;

#Problem 2 - Exer
C1 <- function(x0,x1) 0.06 * x0 + 0.10 * x1 + -3520
C2 <- function(x0,x1) 1 * x0 + -2 * x1

#Problem 3 - Exer

F1 <- function(x0,x1,x2,x3,x4,x5) .7 * x0 + .2 * x1 + .1 * x2 + .1 * x3 + .05 * x4 + 0 * x5 + -1000;
F2 <- function(x0,x1,x2,x3,x4,x5) .1 * x0 + .6 * x1 + 0 * x2 + .1 * x3 + .05 * x4 + .1 * x5 + -700;
F3 <- function(x0,x1,x2,x3,x4,x5) .05 * x0 + .10 * x1 + .75 * x2 + .1 * x3 + 0 * x4 + .05 * x5 + -1300;
F4 <- function(x0,x1,x2,x3,x4,x5) .05 * x0 + 0 * x1 + .05 * x2 + .6 * x3 + .05 * x4 + .05 * x5 + -900;
F5 <- function(x0,x1,x2,x3,x4,x5) 0 * x0 + .05 * x1 + .05 * x2 + .1 * x3 + .75 * x4 + .2 * x5 + -1500;
F6 <- function(x0,x1,x2,x3,x4,x5) .1 * x0 + .05 * x1 + .05 * x2 + 0 * x3 + .1 * x4 + .6 * x5 + -1200;

#Problem 1 - Assignment
Q1 <-function(x0,x1,x2) 1 * x0 + -.0190 * x1 + -.0403 * x2 + -5213.04
Q2 <-function(x0,x1,x2) -.0331 * x0 + .9911 * x1 + -.1935 * x2 + -105182.07
Q3 <-function(x0,x1,x2) -.0003 * x0 + -.0004 * x1 + .9994 * x2 + -124023.7


#Probset Lab #2

Z1 <-function(x0,x1,x2) 1 * x0 + 1 * x1 + 1 * x2 + -25000;
Z2 <-function(x0,x1,x2) .06 * x0 + .07 * x1 + .08 * x2 + -1620;
Z3 <-function(x0,x1,x2) 1 * x1 + -1 * x2 + -6000;

#Z1 <-function(x0,x1,x2) 2 * x0 + 3 * x1 + 1 * x2 + -9;
#Z2 <-function(x0,x1,x2) 3 * x0 + 1 * x1 + 3 * x2 + -14;
#Z3 <-function(x0,x1,x2) 8 * x0 + 5 * x1 + 7 * x2 + -32;

#Z1 <-function(x0,x1) 1 * x0 + 1 * x1 + -1000;
#Z2 <-function(x0,x1) .1 * x0 + .7 * x1 + -160;


#newList = list(E1, E2, E3); #HandOutZA test case
#newList = list(D1,D2,D3, D4); #problem 1
#newList = list(C1,C2) #problem 2
#newList = list(F1, F2, F3, F4, F5, F6); #PROBLEM 3
#newList = list(Q1,Q2,Q3)
newList = list(Z1,Z2,Z3);
matrix = AugCoeffMatrix(newList)

Gaussian<-function(matrix){

  numOfEq = length(newList); #number of equations
  newMatrix = AugCoeffMatrix(newList) #Augmented Coefficient Matrix
  unknowns = length(newMatrix$variables) #number of unknown variables
  
  variables = c(newMatrix$variables) #creates a vector of the variables
  variableValues = rep(0, (unknowns)) #creates a vector of zeroes corresponding to variables
  
  for(pivotIteration in 1:(numOfEq-1)){ #iterates all values below the main diagonal
    max = 0
    pivotRow = 0
    PivotElement = 0;

    for (i in pivotIteration:(numOfEq)){#loop from the value below the pivot element to the number of equations.
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
    
    #pivot element in main diagonal of the matrix
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration]

    #if division by zero, prompt    
    if(PivotElement == 0){
      print("Division by zero error. No solution.")
      return(variableValues)
    }
    
    
    for (upperTriangle in pivotIteration:(numOfEq-1)){ #looping all values to be eliminated
      #Value to be eliminated
      VTBE = newMatrix$augcoeffmatrix[upperTriangle+1,pivotIteration]
      #Multiplier
      multiplier = VTBE/PivotElement
      
      #Variable for the temporary pivot rows
      tempVector = c(newMatrix$augcoeffmatrix[pivotIteration,])
      tempVector = tempVector*multiplier #Scalar multiplication to rows

      #row - vector(M x PivotRow) loop
      for(element in 1:(unknowns+1)){
        newMatrix$augcoeffmatrix[upperTriangle+1, element] = newMatrix$augcoeffmatrix[upperTriangle+1, element] - (tempVector[element])
      }
    }
  }

  for(var in (unknowns):1){     #for every unknown
    for(coeffs in (unknowns+1):1){ #for every column
      if(coeffs != var){            #if the coefficient is not in the main diagonal, add.
        if(coeffs != (unknowns+1)){ #if the value to be added is not in the right hand side, multiply to its coefficient before adding.
          variableValues[var] = variableValues[var] + (-1)*newMatrix$augcoeffmatrix[var,coeffs]*variableValues[coeffs]
        }
        else{               #if the value to be added is in the right hand sige, add only
          variableValues[var] = variableValues[var] + newMatrix$augcoeffmatrix[var,coeffs]
        }
      }
    }

    #divide both sides by the coefficient of the unknown variable
    variableValues[var] = variableValues[var]/newMatrix$augcoeffmatrix[var,var]
  }
  
  #create new list for results
  results = list(variables = newMatrix$variables, augcoeffmatrix = newMatrix$augcoeffmatrix, values = variableValues)
  print(variableValues)
  return(results)
}

GaussJordan<-function(matrix){

  numOfEq = length(newList); #number of equations
  newMatrix = AugCoeffMatrix(newList) #Augmented Coefficient Matrix
  unknowns = length(newMatrix$variables) #number of unknown variables
  
  variables = c(newMatrix$variables) #creates a vector of the variables
  variableValues = rep(0, (unknowns)) #creates a vector of zeroes corresponding to variables

  for(pivotIteration in 1:(numOfEq)){ #loops through the rows of the system
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
 
    #Pivot element before normalization
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration] #pivot element in main diagonal
    
    if(PivotElement == 0){
      print("Division by zero error. No solution.")
      return(variableValues)
    }
    
    newMatrix$augcoeffmatrix[pivotIteration,] = newMatrix$augcoeffmatrix[pivotIteration,]/PivotElement
    
    #pivot element after normalization
    PivotElement = newMatrix$augcoeffmatrix[pivotIteration,pivotIteration]
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
      }
    }
  }
  
  for(i in 1:(unknowns)) #input results to a vector
    variableValues[i] = newMatrix$augcoeffmatrix[i,(unknowns+1)]
  
  #create new list for results
  results = list(variables = newMatrix$variables, augcoeffmatrix = newMatrix$augcoeffmatrix, values = variableValues)
  print(variableValues)
  return(results)
}

#function call
#GaussJordan(matrix)
cat("\n\n\n\n\n")
Gaussian(matrix)