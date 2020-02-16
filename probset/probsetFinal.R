source('~/Desktop/CMSC150/Week5/tumpalan_exer04.R') #Source file for Gaussian and Gauss-Jordan

#Data sets for the x, y, and z components of the joints.

x = c(0,6.5,6.5,-1,4.1,4.1,7.2,7.2,7.2,0.5,9.86)
y = c(0,2.15,-2.15,0,1.63,-1.63,3.26,0,-3.26,0,0)
z = c(0,0,0,5.54,5.54,5.54,5.54,5.54,5.54,7.83,6.54)

#list of vectors for the beams connecting 2 different joints
beamConnections = list(c(1,4),c(1,5),c(1,6),c(2,5),c(2,7),c(2,8),c(3,8),c(3,9),c(3,6),c(4,5),c(4,6),c(5,6),c(5,7),c(7,8),c(5,8),c(6,8),c(8,9),c(6,9),c(4,10),c(5,10),c(7,10),c(8,10),c(9,10),c(6,10),c(10,11),c(7,11),c(9,11))

distance<-function(j,k){
  
  #Data sets for the x, y, and z components of the joints.

  x = c(0,6.5,6.5,-1,4.1,4.1,7.2,7.2,7.2,0.5,9.86)
  y = c(0,2.15,-2.15,0,1.63,-1.63,3.26,0,-3.26,0,0)
  z = c(0,0,0,5.54,5.54,5.54,5.54,5.54,5.54,7.83,6.54)

  
  value = (x[k] - x[j])^2 + (y[k] - y[j])^2 + (z[k] - z[j])^2 #Distance formula.
  
  return (sqrt(value)) #return square root of the value to complete the distance formula
}


angleCosine<-function(j,k,component){
  
  #cos(k;j;x) = x(k) - x(j) / dist(k,j)
  cosjkx = (component[k] - component[j])/distance(j,k)
  
  return (cosjkx)
}

solve<-function(a,b,c,connections){
  
  #main list for the data
  mainData = list(mainMatrix = matrix(0L, byrow = TRUE, ncol = 28, nrow = 27, dimnames = list(c(1:27), c(paste("f",c(1:27), sep="", collapse = NULL),"Load"))), unknowns = rep(0, 27) )
  
  #variables for looping
  joint = 1
  rowInput = 1

  #for the 33 equations that the 11 joints would produce
  for(row in 1:33){
    
    if(row%%3 == 1){      #loops for the direction. Loops to X Y and Z points even joint exceptions will be made.
      direction = a;
      directionNum = 1;  #variable for direction identification, 1 == x
    } else if(row %%3 == 2){
      direction = b;
      directionNum = 2; #variable for direction identification, 2 == y
    } else{
      direction = c;
      directionNum = 3;   #variable for direction identification, 3 == z
    }
    
    if(joint == 11){ #right hand side initialization of the given joint 11 forces
      if(directionNum == 1){
        mainData$mainMatrix[rowInput, 28] = 0
      } else if(directionNum == 2){
        mainData$mainMatrix[rowInput, 28] = -1 #multiplied to -1 since RHS is negative load
      } else if(directionNum == 3){
        mainData$mainMatrix[rowInput, 28] = 10 #multiplied to -1 since RHS is negative load
      }
    }
    
    for(element in 1:27){ #loop through the connections
      #exceptions for the 6 joint equations to be dropped
      if((joint == 1 && directionNum == 2)||(joint == 1 && directionNum == 3)||(joint == 2 && directionNum == 1)||(joint == 2 && directionNum == 3) || (joint == 3 && directionNum == 1) || (joint == 3 && directionNum == 3)){
        rowInput = rowInput - 1 #decrement row input variable to avoid skipping of rows on coefficient insertion on the matrix
        break;
      }
      
      #loop through the connections and input all cos(j,k,component) values to the matrix
      if(connections[[element]][1] == joint){
        mainData$mainMatrix[rowInput,element] = angleCosine(joint, connections[[element]][2], direction)
      }else if(connections[[element]][2] == joint){
        mainData$mainMatrix[rowInput,element] = angleCosine(joint, connections[[element]][1], direction)
      }
  
    }
    
    rowInput = rowInput + 1
    
    if(directionNum == 3){
      joint = joint + 1 #if the Z component of a joint is done, proceed to the next joint
    }
  }

  #since cos(j;k;x) = -cos(k;j;x), all RHS will be equal to 0
  
  #store unknowns' values to the main Data
  mainData$unknowns = Gaussian(mainData$mainMatrix)

  return (mainData)
}

finalData = solve(x,y,z,beamConnections)
print(finalData)