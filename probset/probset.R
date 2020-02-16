source('~/Desktop/CMSC150/Week4/tumpalan_exer03_v2.R')

x = c(0,6.5,6.5,-1,4.1,4.1,7.2,7.2,7.2,0.5,9.86)
y = c(0,2.15,-2.15,0,1.63,-1.63,3.26,0,-3.26,0,0)
z = c(0,0,0,5.54,5.54,5.54,5.54,5.54,5.54,7.83,6.54)


beamConnections = list(c(1,4),c(1,5),c(1,6),c(2,5),c(2,7),c(2,8),c(3,8),c(3,9),c(3,6),c(4,5),c(4,6),c(5,6),c(5,7),c(7,8),c(5,8),c(6,8),c(8,9),c(6,9),c(4,10),c(5,10),c(7,10),c(8,10),c(9,10),c(6,10),c(10,11),c(7,11),c(9,11))

#beamConnections = list(c(1,4),c(1,5),c(1,6),c(2,5),c(2,7),c(2,8),c(3,8),c(3,9),c(3,6),c(4,5),c(4,6),c(4,10),c(5,6),c(5,7),c(5,8),c(5,10), c(6,8), c(6,9),c(6,10),c(7,8), c(7,10),c(7,11),c(8,9),c(8,10),c(9,10),c(9,11),c(10,11))

distance<-function(j,k){
  x = c(0,6.5,6.5,-1,4.1,4.1,7.2,7.2,7.2,0.5,9.86)
  y = c(0,2.15,-2.15,0,1.63,-1.63,3.26,0,-3.26,0,0)
  z = c(0,0,0,5.54,5.54,5.54,5.54,5.54,5.54,7.83,6.54)
  
  value = ((x[k] - x[j])^2) + ((y[k] - y[j])^2) +((z[k] - z[j])^2)
  return (sqrt(value))
}


angleCosine<-function(j,k,component1){
  #cos(k;j;x) = x(k) - x(j) / dist(k,j)
  cosjkx = (component1[k] - component1[j])/distance(j,k)
  
  return (cosjkx)
}

angleCosine(11,10,x)

fillMatrix<-function(a,b,c,connections){
  
  mainData = list(mainMatrix = matrix(0L, byrow = TRUE, ncol = 28, nrow = 27, dimnames = list(c(1:27), c(paste("f",c(1:27), sep="", collapse = NULL),"Load"))))
  
  joint = 1
  equationInput = 1
  for(row in 1:33){
    if (row %% 3 == 1){
      direction = a;
      directionNum = 1
    } else if (row %% 3 == 2){
      direction = b;
      directionNum = 2
    } else{
      direction = c;
      directionNum = 3
    }
    
    for(col in 1:27){
      if((joint == 1 && row%%3 == 2) || (joint == 1 && row%%3 == 0) || (joint == 2 && row%%3 == 0) || (joint == 2 && row%%3 == 1) || (joint == 3 && row%%3 == 0) || (joint == 3 && row%%3 == 1)){
     #   cat("Joint: ")
    #    cat(joint)
    #    cat("\n Direction: ")
    #    cat(direction)
    #    cat("\n")
    #    print("Equation excluded.")
        equationInput = equationInput - 1
        break;
      }
      else if(connections[[col]][1] == joint || connections[[col]][2] == joint){
         mainData$mainMatrix[equationInput,col] = angleCosine(connections[[col]][1], connections[[col]][2], direction)
         cat("J:")
         cat(connections[[col]][1])
         cat("\nK:")
         cat(connections[[col]][2])
         cat("\nCos: ")
         print(mainData$mainMatrix[equationInput,col])
         if(joint == 11){
           if(directionNum == 1){
             mainData$mainMatrix[equationInput,28] = 0
           } else if(directionNum == 2){
             mainData$mainMatrix[equationInput,28] = -1
           } else{
             mainData$mainMatrix[equationInput,28] = 10
           }
         }
      }
    }
    
    equationInput = equationInput + 1
    direction = direction + 1
    
   # for(RHS in 1:27){
    #  for(RHS2 in 1:27){
    #    mainData$mainMatrix[RHS,28] = mainData$mainMatrix[RHS,28] + mainData$mainMatrix[RHS,RHS2]   
    #    round(mainData$mainMatrix[RHS,28], digits = 2)
    #  }
    #}
  }
  
  
  return (mainData$mainMatrix)
  
}

finalData = fillMatrix(x,y,z,beamConnections)

unknowns = Gaussian(finalData)

#print(finalData)
#print(unknowns)