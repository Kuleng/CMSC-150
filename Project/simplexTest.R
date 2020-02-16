source('~/Desktop/CMSC150/Project/GaussMod.R')

vector1 = c(310,260,280,180,80,200,160,220)
matrix1 = matrix(c(100,10,8,6,5,4, 100, 6,5,4,3,6,100,3,4,5,5,9), nrow = 3, byrow = TRUE)

simplexSetup<-function(demand_supply_vector, costMatrix, rType){
  
  #Sets up the initial matrix based from the initial positions of the constraints and equations
  mainMatrix = matrix(0L, nrow = 9, ncol = 16)
  
  #sets up the supply constraints
  mainMatrix[1,] = c(-1,-1,-1,-1,-1, 0,0,0,0,0,0,0,0,0,0, demand_supply_vector[1]*(1))
  mainMatrix[2,] = c(0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0, demand_supply_vector[2]*(1))
  mainMatrix[3,] = c(0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1, demand_supply_vector[3]*(1))
  
  
  #Sets up the demand constraints to the matrix
  col = 4
  i = 1
  demandInput = 4
  while(col < 9){
      while(i <= 15){
        mainMatrix[col, i] = 1
        i = i + 5
      }
    mainMatrix[col, 16] = demand_supply_vector[demandInput]*(-1)
    demandInput = demandInput + 1
    col = col + 1
    i = col - 3
  }
  
  #mainMatrix[9,] = c(5,6,7,8,9,6,7,8,9,10,3,5,7,11,13, 0)
  
  #Gets the initial cost values per shipping area to the matrix
  newMatrix <- matrix(costMatrix[,-1], nrow = 3, byrow = FALSE)
  costVector <- c()
  
  input = 1
  for(i in 1:(length(newMatrix[,1]))){
    for(j in 1:(length(newMatrix[1,]))){
      costVector[input] = newMatrix[i, j]
      input = input + 1
    }
  }
  
  #print(costVector)
  mainMatrix[9,] = c(costVector, 0L)
  
  print(mainMatrix)
  
  #Tranposes the initial matrix
  transposedMatrix = t(mainMatrix)
  
  #creates a matrix of slack variables for the identification of unknowns
  unknownsMatrix = matrix(0L, nrow = 16, ncol = 16)
  for(i in 1:(length(unknownsMatrix[,1]))){
    for(j in 1:(length(unknownsMatrix[1,]))){
      if(i == j)
        if(i > 3)
          unknownsMatrix[i,j] = 1
        else
          unknownsMatrix[i,j] = 1
    }
  }
  
  #binds the matrix of slack variables to the initial matrix
  transposedMatrix = cbind(transposedMatrix, unknownsMatrix)
  
  #Sort the merged matrix to adjust RHS to the rightmost side of the matrix
  for(int in 9:(length(transposedMatrix[1,])-1)){
    temp = transposedMatrix[,int]
    transposedMatrix[,int] = transposedMatrix[,int+1]
    transposedMatrix[,int+1] = temp
  }
  
  #Perform Modified Gauss Jordan for Simplex Approach
  mainSolution = GaussMod(transposedMatrix, 1)
  matrixSolution = matrix(0L, nrow = length(mainSolution), ncol = length(mainSolution[[1]]) +1)
  
  colnames(matrixSolution) <- 1:(length(mainSolution[[1]]) +1)
  for(row in 1:(length(mainSolution))){
    matrixSolution[row,1] = paste("Iteration ", row, sep="")
    for(solution in 1:(length(mainSolution[[1]]))){
      matrixSolution[row,solution+1] = mainSolution[[row]][solution]
    }
  }

  mainTableaus = GaussMod(transposedMatrix, 2)
  
  #print(mainTableaus[[length(mainTableaus)]][16,25])
  if(rType ==1){
    print(mainSolution)
    print(matrixSolution)
  }
  else
    print(mainTableaus)
  
}

#print(matrix1)
simplexSetup(vector1,matrix1, 1)
simplexSetup(vector1,matrix1, 2)
#print(GaussMod(transposedMatrix))

