fit_distr_edge <- function(N, filt_ed){
  
  x <- 1:length(N)
  y <- sort(N,decreasing = T)
  
  my_func <- function(x,b) {
    return((max(N)-filt_ed)*exp(-b * x) + filt_ed)
  }
  start_values <- c(b=0.001)
  # Esegui il fitting dei parametri utilizzando nls()
  fit <- nls(y ~ my_func(x,b), start = start_values)
  
  # Estrai i parametri stimati
  parameters <- coef(fit)
  
  x.new <- seq(1,length(N))
  y.new <- (max(N)-filt_ed) *exp(-parameters* x.new) + filt_ed
  
  
  err <- y.new - y
  
  func_err <- function(x,b) {
    return(abs(min(err))*exp(-b * x))
  }
  
  y.err <- abs(err[which.min(err):length(err)])
  x.err  <- 1:length(y.err)
  
  start_values.err  <- c(b=0.001)
  fit.err <- nls(y.err  ~ func_err(x.err,b),   start = start_values.err)
  
  
  
  x.final <- seq(filt_ed,length(N), 
                 by = (length(N)-filt_ed)/10000)
  y.fit <- my_func(x.final, coef(fit))
  y.err <- func_err(x.final, coef(fit.err))
  start_err <- which.min(err)*length(x.final)/length(N)
  y.final <- y.fit + y.err 
  
  return(y.final)

}
