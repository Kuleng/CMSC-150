E1 <- function(x0,x1,x2) 0.3*x0+-0.2*x1+10*x2+-71.4;
E2 <- function(x0,x1,x2) 3*x0+-0.2*x2+-0.1*x1+-7.85;
E3 <- function(x0,x1,x2) 0.1*x0+7*x1+-0.3*x2+19.3;

mainList = list(E1,E2,E3)

AugCoeffMatrix <-function(mainList){
  #checking if unknowns are equal to the number of 
  numOfEq = length(mainList)
  checked = 0
  loopCount = 0
  listVars = list(c())
  mainVarCount = 1

  for (a in mainList){
    tempDeparse = deparse(a, width.cutoff = 200L); #deparsed function
    tempCoeffs = strsplit(tempDeparse[2],split = "+", fixed = T) #splitted the string that was deparsed to elements
    if (numOfEq == (length(tempCoeffs[[1]]) -1 )){
      checked = checked + 1
    }
    loopCount = loopCount + 1
  }
  
  if (checked == numOfEq){
    print("Number of equations and unknowns are equal.")
  } else{
    print("Number of equations and uknowns are not equal.")
  }

  for (a in mainList){
    tempDeparse = deparse(a, width.cutoff = 200L); #deparsed function
    tempCoeffs = strsplit(tempDeparse[2],split = "+", fixed = T) #splitted the string that was deparsed to elements
    varCount = 1 #counter for the loop
    while (varCount < length(tempCoeffs[[1]])){
      tempSplittedCoeffs = strsplit(tempCoeffs[[1]][varCount], split = "*", fixed = T) #splits the elements to coefficient and variable
      listVars[mainVarCount] = gsub(" ", "", tempSplittedCoeffs[[1]][2], fixed = T) #the variable will be stored in the variable list
      varCount = varCount + 1
      mainVarCount = mainVarCount + 1
    }
  }
  
  mainAug = matrix(nrow = numOfEq, ncol = length(listVars), dimnames = list(c(1:numOfEq),c(listVars)))

  for (a in mainList){
    tempDeparse = deparse(a, width.cutoff = 200L); #deparsed function
    tempCoeffs = strsplit(tempDeparse[2],split = "+", fixed = T) #splitted the string that was deparsed to elements
    varCount = 1 #counter for the loop
    while (varCount < length(tempCoeffs[[1]])){
      tempSplittedCoeffs = strsplit(tempCoeffs[[1]][varCount], split = "*", fixed = T) #splits the elements to coefficient and variable
      mainAug[varCount][varCount] = gsub(" ", "", tempSplittedCoeffs[[1]][1], fixed = T) 
      varCount = varCount + 1
      mainVarCount = mainVarCount + 1
    }
  }
  print(listVars)
  print(mainAug)
  return (checked)
}