functionNum = 1
E2 <- function(x0,x1,x2) 8 * x0 + -3 * x1 + 4 * x2 + -12;
E3 <- function(x0,x1,x2) 6 * x0 + -3 * x1 + 1 * x2 + -8;

mainList = list()

for(i in 1:2){
  if(i == 1)
    mainList[[1]] <- E2
  else
    mainList[[2]] <- E3
}

print(mainList)