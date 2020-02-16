#data points given
#x = c(8,9,11,12)
#y = c(0.9031, 0.9542, 1.0414, 1.0792)

#Exer8 word problem
x = c(0.10377, 0.11144, 0.1254)
y = c(6.4147, 6.5453, 6.7664)

list = list(x,y)

Lagrange<-function(list){
  #data X will be [[1]] while Y will be [[2]]
  polynomial_string = "function(x) "
  
  #numerator loop
  for(term in 1:(length(list[[1]]))){
    
    polynomial_string = paste(polynomial_string, list[[2]][term], " * ",  sep="", collapse = NULL)
    polynomial_string = paste(polynomial_string, "((", sep ="", collapse = NULL)
    termAdded = 1
    for(numeratorTerm in 1:(length(list[[1]]))){
      if(numeratorTerm != term  && termAdded <= length(list[[1]])-1){
        if(termAdded == (length(list[[1]]) - 1)){
          polynomial_string = paste(polynomial_string, "(x - ",list[[1]][numeratorTerm],")", sep ="", collapse = NULL)
          termAdded = termAdded + 1
        }
        else{
          polynomial_string = paste(polynomial_string, "(x - ",list[[1]][numeratorTerm],")*", sep ="", collapse = NULL)
          termAdded = termAdded + 1
        }
      }
    }
    
    polynomial_string = paste(polynomial_string, ")/", sep = "", collapse = NULL)

    
    polynomial_string = paste(polynomial_string, "(",sep = "", collapse = NULL)
    
    termAdded = 1
    for(denominatorTerm in 1:(length(list[[1]]))){
      if(denominatorTerm != term && termAdded <= length(list[[1]]) - 1){
        if(termAdded == (length(list[[1]]) - 1)){
          polynomial_string = paste(polynomial_string, "(" , list[[1]][term], " - ",list[[1]][denominatorTerm],")", sep ="", collapse = NULL)
          termAdded = termAdded + 1
        }else{
          polynomial_string = paste(polynomial_string, "(", list[[1]][term], " - ",list[[1]][denominatorTerm],")*", sep ="", collapse = NULL)
          termAdded = termAdded + 1
        }
      }
    }
    polynomial_string = paste(polynomial_string, "))",sep = "", collapse = NULL)
    
    if(term != length(list[[1]]))
      polynomial_string = paste(polynomial_string, " + ",sep = "", collapse = NULL)
    
  }
  
  #evaluate the polynomial string
  polynomial_function = eval(parse(text = polynomial_string))
  
  #insert to a labeled list
  mainData = list(polynomial_string = polynomial_string, polynomial_function = polynomial_function)
  
  #return labeled list
  return(mainData)
}

Neville<-function(x , list){
  y_label = c("xi", "|x-xi|")
  
  for(Pi_iteration in 1: length(list[[1]])){
    y_label = c(y_label, paste("Pi",Pi_iteration, sep="", collapse = NULL))
  }

  newMatrix = matrix(0L , nrow = length(list[[1]]), ncol = (length(y_label)), byrow = T)
  neville_matrix = matrix(0L , nrow = length(list[[1]]), ncol = (length(y_label)), byrow = T)
  colnames(neville_matrix) = y_label
  
  for(y in 1:(length(list[[2]]))){
    newMatrix[y, 1] = list[[1]][y]
    newMatrix[y, 2] = abs(x - list[[1]][y])
    newMatrix[y, 3] = list[[2]][y]
  }
  
  #perform swapping
  sortedDistance = c(order(newMatrix[,2])) #create a vector of sorted distance of xi to x
  for(i in 1:(length(list[[2]]))){
    neville_matrix[i,] = newMatrix[sortedDistance[i],]    #insert the matrix data of new matrix to the original matrix
  }
  
  
  #perform value of Pi using the NIP formula
  for(col in 1:(length(y_label) - 3)){ #for every column
    for(row in 1:(length(list[[1]]) - col)){ #for every row
      neville_matrix[row, col+3] = ((x - neville_matrix[row, 1])*neville_matrix[row+1, ((col+3)-1)] + (neville_matrix[row+col, 1] - x)*neville_matrix[row, (col+3)-1])/(neville_matrix[row+col, 1] - neville_matrix[row,1])
    }    
  }

  #convert to numeric
  y_val = as.numeric(neville_matrix[1,length(list[[2]])+2])
  
  #insert to a labeled list
  mainData = list(table = neville_matrix, y_value = y_val)

  #return labeled list
  return(mainData)
}

print(Lagrange(list))
print(Neville(.108, list))