#TUMPALAN, JOHN KARL B.
#2018 - 02385
# WEEK 9 - ROOTS OF EQUATION

SecantMethod<-function(f,x0, x1, macheps, max, verbose){
  
  #SET INITIAL VALUES
  iteration = 1;
  ea = 100;
  current_x0 = x0;
  current_x1 = x1;
  f_x0 = f(current_x0);
  f_x1 = f(current_x1);
  current_x2 = ((current_x1)-(f_x1)*((current_x1-current_x0)/(f_x1-f_x0))); #x2
  f_x2 = f(current_x2);
  
  #CREATE A MATRIX WITH PRE-SET COLUMN NAMES BASED FROM THE INITIAL VALUES
  mainMatrix = matrix(c(iteration, current_x0, current_x1, f_x0, f_x1, current_x2, f_x2, ea), ncol = 8)
  colnames(mainMatrix) = c("Iteration", "x0", "x1", "f(x0)", "f(x1)", "x2", "f(x)", "ea")
  
  while(iteration < max && ea > macheps){
    
    #CHANGE VALUES FOR THE NEXT ITERATION
    current_x0 = current_x1;
    current_x1 = current_x2;
    
    #UPDATE VARIABLE VALUES USING THE NEW VALUES OF X0 AND X1
    f_x0 = f(current_x0);
    f_x1 = f(current_x1);
    current_x2 = ((current_x1)-(f_x1)*((current_x1-current_x0)/(f_x1-f_x0))); #x2
    f_x2 = f(current_x2);
    ea = abs((current_x2 - current_x1)/(current_x2))*100

    #LOOP UPDATE
    iteration = iteration + 1
    
    #MATRIX UPDATE
    mainMatrix = rbind(mainMatrix,c(iteration, current_x0, current_x1, f_x0, f_x1, current_x2, f_x2, ea))
    
  }
  if(verbose == TRUE)
    print(mainMatrix)
  
  #CREATE A LABELED LIST
  mainData = list('f' = f, given_x0 = x0, given_x1 = x1, x = current_x2, iterations = iteration, ea = ea)
  
  return(mainData)
}

MullerMethod<-function(f,x0, x1, x2, macheps, max, verbose){
  
  #SET INITIAL VALUES
  iteration = 1;
  
  current_x0 = x0
  current_x1 = x1
  current_x2 = x2
  current_h0 = current_x1-current_x0
  current_h1 = current_x2-current_x1
  current_d0 = (f(current_x1)-f(current_x0))/current_h0
  current_d1 = (f(current_x2)-f(current_x1))/current_h1
  current_A = (current_d1-current_d0)/(current_h1+current_h0)
  current_B = (current_d1+current_A*current_h1)
  current_C = (f(current_x2))
  current_denominatorA = current_B+sqrt((current_B)^2-(4*current_A*current_C))
  current_denominatorB = current_B-sqrt((current_B)^2-(4*current_A*current_C))
  
  #CONDITIONALS FOR DENOMINATOR IN THE FORMULA OF X3
  if(current_denominatorA > current_denominatorB)
    current_x3 = current_x2-((2*current_C)/current_denominatorA)
  else
    current_x3 = current_x2-((2*current_C)/current_denominatorB)
  
  #RELATIVE ERROR
  ea = abs(((current_x3-current_x2)/current_x3)*100)
  
  #CREATE A MATRIX BASED FROM THE INITIAL VALUES AND CREATE A PRE-SET COLUMN NAMES
  mainMatrix = matrix(c(iteration, current_x0,current_x1,current_x2, f(current_x0), f(current_x1), f(current_x2), current_A, current_B, current_C, current_x3, f(current_x3), ea), ncol = 13)
  colnames(mainMatrix) = c("i", "x0", "x1", "x2", "f(x0)", "f(x1)", "f(x2)", "A","B","C", "x3", "f(x3)", "ea")
  
  while(iteration < max && ea > macheps){
    
    #CHANGE VALUES FOR NEXT ITERATION
    current_x0 = current_x1
    current_x1 = current_x2
    current_x2 = current_x3
    
    #RE-COMPUTE VALUES BASED FROM NEW X0, X1, AND X2
    current_h0 = current_x1-current_x0
    current_h1 = current_x2-current_x1
    current_d0 = (f(current_x1)-f(current_x0))/current_h0
    current_d1 = (f(current_x2)-f(current_x1))/current_h1
    current_A = (current_d1-current_d0)/(current_h1+current_h0)
    current_B = (current_d1+current_A*current_h1)
    current_C = (f(current_x2))
    current_denominatorA = current_B+sqrt((current_B)^2-(4*current_A*current_C))
    current_denominatorB = current_B-sqrt((current_B)^2-(4*current_A*current_C))
    if(current_denominatorA > current_denominatorB)
      current_x3 = current_x2-((2*current_C)/current_denominatorA)
    else
      current_x3 = current_x2-((2*current_C)/current_denominatorB)
    
    #ERROR
    ea = abs(((current_x3-current_x2)/current_x3)*100)

    #UPDATE
    iteration = iteration + 1
    
    #MATRIX UPDATE
    mainMatrix = rbind(mainMatrix, c(iteration, current_x0,current_x1,current_x2, f(current_x0), f(current_x1), f(current_x2), current_A, current_B, current_C, current_x3, f(current_x3), ea))
    }
  
  #CREATE A LABELED LIST
  mainData = list('f' = f, given_x0 = x0, given_x1 = x1, given_x2 = x2, x3 = current_x3, iterations = iteration, ea = ea)
  
  #PRINT IF NEEDED
  if(verbose == TRUE){
    print(mainMatrix)
  }
  return(mainData)
}

#print(SecantMethod(function(x) sin(x) + cos(1+x^2) - 1, 1.0, 3.0, 1e-9, 100000, TRUE))

print(SecantMethod(function(x) 3.1415*x^3-9*3.1415*x^2+28.6487*3.1415, 0, 6, 1e-5, 5, TRUE))
#print(SecantMethod(function(x) -26 + 85*x - 91*x^2 - 44*x^3 - 8*x^4 + x^5, 10, 15, 1e-5, 1000, TRUE))
#print(MullerMethod(function(x) -26 + 85*x - 91*x^2 - 44*x^3 - 8*x^4 + x^5, 8, 10, 15, 1e-5, 100, TRUE))