cor_fun <- function(list_condizioni,
                    test_type,
                    significance,
                    p.adj,
                    compute_weights,
                    Rcorr = NULL,
                    threshold.n.corr = NULL,
                    flatten = NULL
                    ){
  # Correlazione di spearman e calcolo dei pesi per i grafi
  names_cond <- names(list_condizioni)
  names(names_cond) <- names_cond
  cor_cond  <- lapply(names_cond, function(x){
    COR_COMPUTE(list_condizioni[[x]],
                test_type = test_type,
                significance = significance,
                p.adj = p.adj,
                compute_weights = compute_weights,
                Rcorr = Rcorr[[x]],
                thr_n = threshold.n.corr,
                flatten = flatten[[x]])
  })
  
  # ordinare le interazioni in tutti i gruppi per ordine alfabetico, in modo che 
  # l'ordine sia lo stesso per tutti i gruppi
  cor_cond2 <- lapply(cor_cond,function(x){
    x$ALL <- x$ALL[order(x$ALL$cor_features),]
    return(x)
  })
  
  graph_cor_notsign.cond <- lapply(cor_cond2, function(x){graph_from_correlation(x$ALL ,x$ALL$weights)})
  graph_cor_sign.cond    <- lapply(cor_cond2, function(x){graph_from_correlation(x$SIGN,x$SIGN$weights)})
  
  graph_cor.cond <- append(graph_cor_notsign.cond,graph_cor_sign.cond)
  names(graph_cor.cond)[1:(length(graph_cor.cond)/2)] <- 
    str_c("Notsign_",names(graph_cor.cond)[1:(length(graph_cor.cond)/2)])
  names(graph_cor.cond)[(length(graph_cor.cond)/2+1):length(graph_cor.cond)] <- 
    str_c("Sign_",names(graph_cor.cond)[(length(graph_cor.cond)/2+1):length(graph_cor.cond)])
  
  return(list(cor_cond   = cor_cond2,
              cor_graphs = graph_cor.cond,
              Rcor = lapply(cor_cond, function(x){x$Rcor})))
}

