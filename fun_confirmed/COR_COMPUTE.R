COR_COMPUTE <- function(dataset,
                        test_type,
                        significance,
                        p.adj,
                        compute_weights){
  
  # matrice delle correlazioni
  Rcorr <- rcorr(as.matrix(dataset),
                 type = test_type)
  
  
  
  # flatten to read easily
  flatten <- flattenCorrMatrix(Rcorr$r,Rcorr$P) 
  
  flatten$cor[which(is.na(flatten$p))] <- 0
  flatten$p[which(is.na(flatten$p))] <- 0
  
  flatten.adj <- flatten[order(flatten$p),]
  flatten.adj$p.adj <- p.adjust(flatten.adj$p,p.adj)
  
  flatten_norm <- flatten.adj
  ind_sig <- which(flatten_norm$p.adj < significance)
  
  if (compute_weights){
    ### normalizzazione per togliere gli 1 dal pvalue, per rendere possibile il logaritmo
    # L'obiettivo è trasformare l'intervallo del pvalue da (Significativià ; 1 ] a (0 ; 1)
    
    # Il nuovo massimo del pvalue sarà la media tra 1 e il secondo massimo
    if (max(flatten_norm$p.adj)==1){
      new_max <- mean(c(1,max(flatten_norm$p.adj[-which(flatten_norm$p.adj==1)])))
    } else {
      new_max <- max(flatten_norm$p.adj)
    }
    
    # calcolo parametri normalizzazione: (pvalue + epsilon ) * k = p_norm1
    # ( Significance + eps ) * k = significance
    # ( max(pvalue)  + eps ) * k = new_max
    
    # risolvendo questo sistema si ottiene
    eps <- ( new_max - max(flatten_norm$p.adj) ) * significance / ( significance - new_max )
    k   <- ( significance - new_max ) / ( significance - max(flatten.adj$p.adj) )
    
    # applichiamo la normalizzazione solo ai pvalue non significativi
    ind_notsig <- which(flatten_norm$p.adj > significance)
    norm1 <- (flatten_norm$p.adj[ind_notsig] + eps) * k
    
    # adesso applichiamo -log10 per "invertire" i pvalue bassi in pesi alti e viceversa 
    lg    <- -log10(norm1)
    
    ## normalizziamo lg tra 0 e 1, calcolo parametri ( lg + eps_log) * k_log = p_norm2
    # ( MAX_lg  + eps_log ) * k_log = 1
    # ( min(lg) + eps_log ) * k_log = min(lg)
    MAX_lg <- -log10(significance)
    
    # risolvendo il sistema si ottiene 
    eps_log <- ( min(lg) - MAX_lg * min(lg) ) / ( min(lg) - 1 )
    k_log   <- ( min(lg) - 1 ) / ( min(lg) - MAX_lg)
    norm2   <- (lg + eps_log) * k_log
    
    flatten_norm$p_norm              <- flatten_norm$p.adj
    flatten_norm$p_norm[ind_notsig]  <- norm2
    
    flatten_norm$p_norm2             <- flatten_norm$p.adj
    flatten_norm$p_norm2[ind_notsig] <- norm2^2
    
    flatten_norm$p_norm3             <- flatten_norm$p.adj
    flatten_norm$p_norm3[ind_notsig] <- norm2^3
    
    
    #moltiplichiamo la correlazione per il pnorm3
    flatten_norm$weights <- flatten_norm$cor
    flatten_norm$weights[ind_notsig] <- flatten_norm$cor[ind_notsig]*norm2^3
    
  } else {
    flatten_norm$weights <- flatten_norm$cor
    flatten_norm$weights[-ind_sig] <- 0
  }
  
  flatten_norm$cor_features <- paste(flatten_norm$row, 
                                     "and", 
                                     flatten_norm$column)
  flatten_norm$cor_features_1 <- paste(flatten_norm$column, 
                                       "and",
                                       flatten_norm$row)
  Sign <- flatten_norm[ind_sig,] #storing significative ones
  Sign <- Sign[order(rownames(Sign)),]
  
  results <- list(ALL = flatten_norm,
                  SIGN = Sign, 
                  Rcor = Rcorr)
  return(results)
}

