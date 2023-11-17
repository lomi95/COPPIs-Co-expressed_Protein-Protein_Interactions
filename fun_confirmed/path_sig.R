path_sig <- function(tb_k,
                     g1,g2,
                     tb_meancor,
                     significance,
                     p.adj,
                     sd_groups,
                     asym.dist){
  g2_0 <- which(tb_meancor[,g2+1]==0)
  tb_meancor[g2_0,g1+1] <- .Machine$double.xmin
  if (tb_k[g1,g2] > 1){
    k_x <- 1
  } else {
    k_x <- tb_k[g1,g2]
  }
  
  if (asym.dist){
    sd.x <- (k_x/abs(qnorm(significance/2)))^2
  } else {
    sd.x <- sd_groups[g1]^2+sd_groups[g2]^2
  }
  ratiog1g2 <- cbind(tb_meancor[,1],tb_meancor[,g1+1]/tb_meancor[,g2+1])
  p.v <- apply(ratiog1g2,1, function(x){
    return(pnorm(x[2],mean = k_x, sqrt(sd.x/x[1])))
  })
  #p.v[which(ratiog1g2[,2]==0)] <- 0
  p.v.adj <- p.adjust(p.v,p.adj)
  p.sign <- p.v.adj[intersect(which(p.v.adj < significance),which(ratiog1g2[,2] < k_x))]
  
  if (length(p.sign)>0){
    
    scores <- matrix(nrow = length(p.sign), ncol = 3)
    n <- 1
    for (i in names(p.sign)){
      # trovare la distanza minima tra il punto e la curva 
      
      y0 <- ratiog1g2[i,2]
      x0 <- ratiog1g2[i,1]
      d0 <- k_x - abs(qnorm(significance/2))*sqrt((sd.x)/x0) - y0
      
      # se è significativo ma la media-qnorm(sig)*sd è negativa
      if ((d0+y0) < 0){
        print("eh?")
        scores[n,1] <- names(p.sign)[n]
        scores[n,2] <- p.sign[n]
        scores[n,3] <- 0
      } else {
        if (x0 == 1){
          d.i <- d0
        } else {
          d.i <- d0 + 1
          x.i <- x0
          while (x.i != 1){
            x.i <- x.i - 1
            y.i <- abs(k_x - abs(qnorm(significance/2))*sqrt((sd.x)/x.i) - y0)
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
      }
      n <- n+1
    }
    colnames(scores) <- c("Description","p.value","Score")
    scores <- as.data.frame(scores)
    scores[,3] <- as.numeric(scores[,3])
    scores[,2] <- as.numeric(scores[,2])
    return(scores)
  }
}
