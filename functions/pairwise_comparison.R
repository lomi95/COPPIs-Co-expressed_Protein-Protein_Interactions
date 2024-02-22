pairwise_comparison <- function(tab_k,
                                g1,g2,
                                sd_groups,
                                tab_ppicor,
                                significance.coppi,
                                correction.coppi){
  sd.x <- sd_groups[g1]^2+sd_groups[g2]^2
  if (tab_k[g1,g2] > 1){
    k_x <- 1
    sd.x <- sd.x - (tab_k[g1,g2] - 1)
  } else {
    k_x <- tab_k[g1,g2]
  }
  
  
  ratiog1g2 <- cbind(tab_ppicor[,"N.edges"],tab_ppicor[,g1+1]/tab_ppicor[,g2+1])
  p.v <- apply(ratiog1g2,1, function(x){
    return(pnorm(x[2],mean = k_x, sqrt(sd.x/x[1])))
  })
  
  if (is.null(correction.coppi)){
    p.v.adj <- p.v
  } else {
    p.v.adj <- p.adjust(p.v,correction.coppi)
  }
  p.sign <- p.v.adj[p.v.adj < significance.coppi & ratiog1g2[,2] < k_x]
  
  if (length(p.sign)>0){
    n.ppi <- tab_ppicor[names(p.sign),"N.ppi"]
    n.edg <- tab_ppicor[names(p.sign),"N.edges"]
    z.scr <- (k_x - ratiog1g2[names(p.sign),2])/ratiog1g2[names(p.sign),1]
    dci   <- log((n.ppi^2-n.edg^2)*n.edg/(sqrt(n.ppi)))
    
    
    scores1 <- data.frame(description = names(p.sign),
                          N.edges = n.edg,
                          N.PPI   = n.ppi,
                          p.adj   = p.sign,
                          score   = mapply(funscore,n.edg,n.ppi)*sqrt(z.scr))
    
    
  } else {
    scores1 <- NULL
  }
  
  
  
  
  
  sd.x <- sd_groups[g1]^2+sd_groups[g2]^2
  if (tab_k[g2,g1] > 1){
    k_x <- 1
    sd.x <- sd.x - (tab_k[g2,g1] - 1)
  } else {
    k_x <- tab_k[g2,g1]
  }
  ratiog2g1 <- cbind(tab_ppicor[,"N.edges"],tab_ppicor[,g2+1]/tab_ppicor[,g1+1])
  p.v <- apply(ratiog2g1,1, function(x){
    return(pnorm(x[2],mean = k_x, sqrt(sd.x/x[1])))
  })
  
  if (is.null(correction.coppi)){
    p.v.adj <- p.v
  } else {
    p.v.adj <- p.adjust(p.v,correction.coppi)
  }
  p.sign <- p.v.adj[p.v.adj < significance.coppi & ratiog2g1[,2] < k_x]
  
  if (length(p.sign)>0){
    n.ppi <- tab_ppicor[names(p.sign),"N.ppi"]
    n.edg <- tab_ppicor[names(p.sign),"N.edges"]
    z.scr <- (k_x - ratiog2g1[names(p.sign),2])/ratiog2g1[names(p.sign),1]
    dci   <- log((n.ppi^2-n.edg^2)*n.edg/(sqrt(n.ppi)))
    scores2 <- data.frame(description = names(p.sign),
                         N.edges = n.edg,
                         N.PPI   = n.ppi,
                         p.adj   = p.sign,
                         score   = -mapply(funscore,n.edg,n.ppi)*sqrt(z.scr))
    
    
  } else {
    scores2 <- NULL
  }
  
  scores <- rbind(scores1,scores2)
  return(scores)
}
