E1 <- function(x0,x1,x2) 0.3*x0+-0.2*x1+10*x2+9.8;
E2 <- function(x0,x1,x2) 3*x0+-0.2*x2+0.1*x1+-7.85;
E3 <- function(x0,x2,x3) 0.4*x0+7*x3+-0.3*x2+19.3;
E4 <- function(x0,x1,x3) 0.2*x0+8*x1+-0.4*x3+20.3;

mainList = list(E1,E2,E3,E4)

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
    tempDeparse = deparse(a, width.cutoff = 200L); #deparsed function
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
  
  if (unknowns == numOfEq){
    print("Number of equations and unknowns are equal.")
  } else{
    print("Number of equations and uknowns are not equal.")
    return (0)
  }
  
  #create a list of all the unknown variables
  mainVariables = (c())
  a = 0
  while(a<unknowns){
    a = a + 1
    mainVariables[a] = paste("x", a-1, sep="")
  }
  #initialize mainData
  mainData = list(variables = mainVariables, augcoeffmatrix = matrix(0L, byrow = TRUE, nrow = numOfEq, ncol = numOfEq+1,dimnames = list(c(1:numOfEq),c(mainVariables,"RHS"))))
  
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
      if(loopCount == (numOfEq+1) && is.na(as.numeric(substring(splittedCoeffs[[loopCount]][2],2))) == FALSE){ #if last col or right hand side
        mainData$augcoeffmatrix[outerLoopCount,loopCount] = (-1)*as.numeric(substring(splittedCoeffs[[loopCount]][1],1))
      } #if the coefficient and variable before the right hand side is not existent, it catches the NA and produces 0 in the matrix.
      else if(is.na(as.numeric(substring(splittedCoeffs[[loopCount]][2],2)))){
        mainData$augcoeffmatrix[outerLoopCount,loopCount] = 0
      } #given xn of the function element, if n == length-1 of matrix index (means that it has correct position
      else if(as.numeric(substring(splittedCoeffs[[loopCount]][2],2)) == loopCount-1){ #while x^n, if n == col.
        mainData$augcoeffmatrix[outerLoopCount,loopCount] = as.numeric(substring(splittedCoeffs[[loopCount]][1],1))
      } else{
        mainData$augcoeffmatrix[outerLoopCount,loopCount] = 0 #assuming x^n, if n != to the coefficient's variable n
      }
    }
    #the coefficient matrix is now built but those coefficients with misplaced variables
    #are going to be initialized with 0 value.
 
    #re-loop all the elements of the function and automatically place the elements in its original position.
    loopCount = 0
    while ( loopCount < numOfEq+1){
      loopCount = loopCount + 1
      if(loopCount == numOfEq+1 && as.numeric(substring(splittedCoeffs[[length(splittedCoeffs)]][1],1)) != 0){
        mainData$augcoeffmatrix[outerLoopCount,numOfEq+1] = (-1)*as.numeric(substring(splittedCoeffs[[length(splittedCoeffs)]][1],1))
      } 
      else if (loopCount <= length(splittedCoeffs)){
        mainData$augcoeffmatrix[outerLoopCount,as.numeric(substring(splittedCoeffs[[loopCount]][2],2))+1] = as.numeric(substring(splittedCoeffs[[loopCount]][1],1))
      }
    }
  }
  
  return (mainData)
}