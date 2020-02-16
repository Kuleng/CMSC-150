source('~/Desktop/CMSC150/Project/GaussMod.R')

mainMatrix = matrix(0L, nrow = 8, ncol = 16)
mainMatrix[1,] = c(1,1,1,1,1, 0,0,0,0,0,0,0,0,0,0, -100)
mainMatrix[2,] = c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0, -100)
mainMatrix[3,] = c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1, -100)
col = 4
i = 1
while(col < 9){
  while(i <= 15){
    mainMatrix[col, i] = 1
    i = i + 5
  }
  mainMatrix[col, 16] = -20
  col = col + 1
  i = col - 3
}

zeroConst = matrix(0L, nrow = 16, ncol = 16)
for(i in 1:15){
  for(j in 1:16){
    if(i == j)
      zeroConst[i,j] = 1
  }
}

mainMatrix = rbind(mainMatrix, zeroConst)
mainMatrix[24,] = c(5,6,7,8,9,6,7,8,9,10,3,5,7,11,13, 1)

print(mainMatrix)


transposedMatrix = t(mainMatrix)
unknownsMatrix = matrix(0L, nrow = 16, ncol = 15)
for(i in 1:(length(unknownsMatrix[,1]))){
  for(j in 1:(length(unknownsMatrix[1,]))){
    if(i == j)
      unknownsMatrix[i,j] = 1
  }
}

transposedMatrix = cbind(transposedMatrix, unknownsMatrix)

for(int in 24:(length(transposedMatrix[1,])-1)){
  temp = transposedMatrix[,int]
  transposedMatrix[,int] = transposedMatrix[,int+1]
  transposedMatrix[,int+1] = temp
}

#print(transposedMatrix)

print(GaussMod(transposedMatrix))

