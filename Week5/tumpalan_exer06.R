source('~/Desktop/CMSC150/Week5/tumpalan_exer04.R')

#variables
#x = c(0,4,8,12,16,20)
#y = c(67,84,98,125,149,185)

#Exer 6 word problem
#x = c(26.67,93.33,148.89, 315.56)
#y = c(1.35, 0.085, 0.012, 0.00075)

#Exer 7 word problem
#x = c(0,8,16,24,32,40)
#y = c(14.621,11.843,9.870,8.418, 7.305, 6.413)

#Exer8 word problem
#x = c(0.10377, 0.11144, 0.1254)
#y = c(6.4147, 6.5453, 6.7664)

#csv integration
mydata = read.csv("dataframe.csv", sep = ",")

mainList = list(mydata[,1],mydata[,2])


PolynomialRegression<-function(degree, list){

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
  plot(mainList[[1]] , mainList[[2]], pch = 20, col = "red", main = "Function plot", xlab = "X", ylab = "Y")
  
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
  mainData$polynomialFunction = eval(parse(text = mainData$polynomialString))
  
  print(mainData)
  
  print(mainData$polynomialFunction(0.108))
}

PolynomialRegression(1,mainList)