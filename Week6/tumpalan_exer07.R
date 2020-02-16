#x = c(3,5,10,12,15)
#y = c(1.5, 6.5, 0.5, 9.9, 4.3)

#x = c(0,8,16,24,32,40)
#y = c(14.621,11.843,9.870,8.418, 7.305, 6.413)
x = c(0.10377, 0.11144, 0.1254)
y = c(6.4147, 6.5453, 6.7664)


dataList = list(x,y)

NDD<-function(dataList){
  
  #create the augmented coefficient matrix
  augcoeffmatrix = matrix(0L, nrow = length(dataList[[1]]), ncol = (length(dataList[[2]])+1))
  
  #create a labeled list with labeled vector of coefficients, polynomialString, and polynomialFunction
  mainData = list(coefficients = c(), polynomial_string = "function(x) ", polynomial_function = "")

  #initialize the first 2 column of the matrix with the given data points  
  for(data in 1:length(dataList[[1]])){
    augcoeffmatrix[data,1] = dataList[[1]][data]
    augcoeffmatrix[data,2] = dataList[[2]][data]
  }
  
  #create a variable for the number of row loops to reduce every iteration
  rowReduction = 1
  #fill the matrix according to its formula according to the pre-computed and given data points
  for(col in 3:(length(dataList[[1]])+1)){
    for(row in 1:(length(dataList[[1]])-rowReduction)){
      augcoeffmatrix[row,col] = ((augcoeffmatrix[row+1,col-1] - augcoeffmatrix[row,col-1])/(augcoeffmatrix[row+rowReduction,1]-augcoeffmatrix[row,1]))
    }
    rowReduction = rowReduction + 1
  }

  #fill the coefficients vector from the first row of the created matrix  
  for(i in 2:(length(dataList[[2]])+1)){
    mainData$coefficients[i-1] = augcoeffmatrix[1,i]
  }

    
  #generate the polynomial string based on the coefficient solved from the matrix and the data X values
  for(FDD in 1:(length(dataList[[2]]))){
    if(FDD == 1) #if first loop, generate b0, which is the first Y
      mainData$polynomial_string = paste(mainData$polynomial_string, dataList[[2]][FDD], " + ", sep = "", collapse = NULL)
    else{
      #paste operations and expressions
      mainData$polynomial_string = paste(mainData$polynomial_string, mainData$coefficients[[FDD]], sep = "", collapse = NULL)

      for(orderCount in 1:(FDD-1)){
        mainData$polynomial_string = paste(mainData$polynomial_string," * (x - ", dataList[[1]][orderCount] ,")", sep = "", collapse = NULL)
      }
      
      #condition for "+" paste break
      if(FDD < length(dataList[[2]]))
        mainData$polynomial_string = paste(mainData$polynomial_string," + ", sep = "", collapse = NULL)
    }
  }
  
  #create a polynomial function based from the created polynomial string  
  mainData$polynomial_function = eval(parse(text = mainData$polynomial_string))
  print(augcoeffmatrix)
  print(mainData)
  print(mainData$polynomial_function(0.103))
  
  eq = mainData$polynomial_function
  curve(eq,from = 1, to=50, xlab = "x", ylab = "y")
}

NDD(dataList)