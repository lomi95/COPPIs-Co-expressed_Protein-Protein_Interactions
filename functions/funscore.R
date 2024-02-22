funscore <- function(x,y){
  R <- x/y
  if (R<=0.5){
    score <- log((1-R^2)/sqrt(R^3)*x^(5/2))
  } else {
    R1 <- 0.5
    score <- log((1-R1^2)/sqrt(R1^3)*(y/2)^(5/2)) + 
      5/(2*x)*(x - y/2)
  }
  return(score)
}