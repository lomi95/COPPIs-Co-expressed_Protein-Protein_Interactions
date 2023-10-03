path_sig <- function(tb_k,
                     g1,g2,
                     tb_meancor,
                     significance,
                     p.adj,
                     sd_groups){
  g2_0 <- which(tb_meancor[,g2+1]==0)
  tb_meancor[g2_0,] <- 1000
  if (tb_k[g1,g2] > 1){
    k_x <- 1
  } else {
    k_x <- tb_k[g1,g2]
  }
  
  ratiog1g2 <- cbind(tb_meancor[,1],tb_meancor[,g1+1]/tb_meancor[,g2+1])
  p.v <- apply(ratiog1g2,1, function(x){
    return(pnorm(x[2],mean = k_x, sqrt((sd_groups[g1]^2+sd_groups[g2]^2)/x[1])))
  })
  
  p.v.adj <- p.adjust(p.v,p.adj)
  p.sign <- p.v.adj[intersect(which(p.v.adj < significance),which(ratiog1g2[,2] < k_x))]
  
  if (length(p.sign)>0){
    
    scores <- matrix(nrow = length(p.sign), ncol = 3)
    n <- 1
    for (i in names(p.sign)){
      
      # trovare la distanza minima tra il punto e la curva 
      
      
      y0 <- ratiog1g2[i,2]
      x0 <- ratiog1g2[i,1]
      d0 <- abs(k_x - 2*sqrt((sd_groups[g1]^2+sd_groups[g2]^2)/x0) - y0)
      
      if (x0 == 1){
        d.i <- d0
      } else {
        d.i <- d0 + 1
        x.i <- x0
        while (x.i != 1){
          x.i <- x.i - 1
          y.i <- abs(k_x - 2*sqrt((sd_groups[g1]^2+sd_groups[g2]^2)/x.i) - y0)
          d.i <- sqrt((x.i-x0)^2+y.i^2)
          if (d.i > d0) {
            break
          } else {
            d0 <- d.i
          }
        }
      }
      if (x0 < d.i){
        d0 <- x0
      }
      
      scores[n,1] <- names(p.sign)[n]
      scores[n,2] <- p.sign[n]
      scores[n,3] <- d0* (k_x-ratiog1g2[i,2])/k_x
      
      n <- n+1
    }
    colnames(scores) <- c("Description","p.value","Score")
    scores <- as.data.frame(scores)
    scores[,3] <- as.numeric(scores[,3])
    scores[,2] <- as.numeric(scores[,2])
    return(scores)
  }
}
