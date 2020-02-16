E1 <- function(x0,x1,x2) 0.3*x0+-0.2*x1+10*x2+-71.4;
E2 <- function(x0,x1,x2) 3*x0+-0.2*x2+-0.1*x1+-7.85;
E3 <- function(x0,x1,x2) 0.1*x0+7*x1+-0.3*x2+19.3;

mainList = list(E1,E2,E3)

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
  mainData = list(variables = mainVariables, augcoeffmatrix = matrix(c(1:(numOfEq*(numOfEq+1))), byrow = TRUE, nrow = numOfEq, ncol = numOfEq+1,dimnames = list(c(1:numOfEq),c(mainVariables,"RHS"))))
  print(mainData)
  
  for (a in mainList){
    outerLoopCount = outerLoopCount + 1
    loopCount = 0
    tempDeparse = deparse(a, width.cutoff = 200L); #deparsed function
    tempCoeffs = strsplit(tempDeparse[2],split = "+", fixed = T) #splitted the string that was deparsed to elements
    splitCoeffs = gsub(" ", "", tempCoeffs[[1]], fixed = T)
    splittedCoeffs = strsplit(splitCoeffs, split = "*", fixed = T) #Splitted Coefficient and Variables
    while ( loopCount < length(splittedCoeffs)){
      loopCount = loopCount + 1
      if(loopCount == length(splittedCoeffs)){
        mainData$augcoeffmatrix[outerLoopCount,loopCount] = (-1)*as.numeric(substring(splittedCoeffs[[loopCount]][1],1))
      } else{
        mainData$augcoeffmatrix[outerLoopCount,loopCount] = as.numeric(substring(splittedCoeffs[[loopCount]][1],1))
      }
    }
  }
  print(mainData)
  
  return (mainData)
}