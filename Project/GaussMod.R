min_matrix = matrix(0L, nrow = 3, ncol = 6);
min_matrix[1,] = c(2,1,1,0,0,3)
min_matrix[2,] = c(1,2,0,1,0,9)
min_matrix[3,] = c(-8,-8,0,0,1,0)

matrix2 = matrix(0L, nrow =4, ncol = 4)
# matrix2[1,] = c(400,300,200,20000)
# matrix2[2,] = c(300,400,500,25000)
# matrix2[3,] = c(-25000,-27000,-30000, 0)
matrix2[1,] = c(1,0,-1,2)
matrix2[2,] = c(1,1,2,10)
matrix2[3,] = c(1,2,2, 8)
matrix2[4,] = c(-6,-8,-4,0)

GaussMod<-function(mainMatrix, rType){
  
  tableausAdded = 1
  tableauList = list()
  solutionList = list()
  
  highestNeg = 0
  lowestPos = 999999
  isNeg = 0
  numOfIterations = 0
  #check the number of negative numbers in the objective function
  for(num in mainMatrix[length(mainMatrix[,1]),]){
    if(num < 0)
      isNeg = isNeg + 1
  }
  
  while(isNeg > 0){
    isNeg <- 0
    numOfIterations =   numOfIterations + 1
    
    #find maximum pivot column
    for(PC in 1:(length(mainMatrix[1,])-1)){
      if( mainMatrix[length(mainMatrix[,1]), PC] < highestNeg){
        highestNeg = mainMatrix[length(mainMatrix[,1]), PC]
        maxPC = PC
      }
    }

    #find minimum pivot row
    for(PR in 1:(length(mainMatrix[,1])-1)){
      if(mainMatrix[PR, maxPC] != 0 && mainMatrix[PR, length(mainMatrix[1,])]/mainMatrix[PR, maxPC] > 0){
        if( mainMatrix[PR, length(mainMatrix[1,])]/mainMatrix[PR, maxPC] < lowestPos ){
          lowestPos = mainMatrix[PR, length(mainMatrix[1,])]/mainMatrix[PR, maxPC]
          minPR = PR
        }
      }
    }

    #normalize
    mainMatrix[minPR,] =mainMatrix[minPR,]/mainMatrix[minPR, maxPC];

    #choose pivot element    
    PE = mainMatrix[minPR, maxPC]
    
    cat("Pivot element :")
    cat(PE)
    cat("\n")
    
    cat("Result of normalization:\n")
    print(mainMatrix)
    cat("\n")
    
    
    for(VTBE in 1:(length(mainMatrix[,PC]))){
      if(VTBE != minPR && mainMatrix[VTBE,maxPC] != 0){
    
        #get multiplier by value to be eliminated over PE
        multiplier = mainMatrix[VTBE, maxPC]/PE
        cat("Multiplier: ")
        cat(multiplier)
        cat("PE:")
        cat(PE)
        cat("\n")
        
        vectorToSubtract = multiplier*mainMatrix[minPR,]

        mainMatrix[VTBE,] = mainMatrix[VTBE,] - vectorToSubtract
      }
    }
    cat("\n")
    cat("-------")
    cat("\n")
    highestNeg = 0
    lowestPos = 9999
    for(num in mainMatrix[length(mainMatrix[,1]),]){
      if(num < 0)
        isNeg = isNeg + 1
    }
    
    cat("isNeg")
    cat(isNeg)
    cat("\n")
    tableauList[[paste("Iteration ", numOfIterations, sep="")]] = mainMatrix
    newVector = c()
    for(i in 9:23){
      newVector = c(newVector, tableauList[[tableausAdded]][length(mainMatrix[,1]), i])
    }
    
    solutionList[[paste("Iteration ", numOfIterations, sep="")]] = newVector
    tableausAdded = tableausAdded + 1
  }

  if(rType == 1)
    return(solutionList)
  else if(rType == 2)
    return(tableauList)
}


# print(GaussMod(matrix2))