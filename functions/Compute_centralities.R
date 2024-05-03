Compute_centralities <- function(gCOR.groups,
                                 interactome,
                                 signifCorr,
                                 correctionCorr,
                                 compute_weights,
                                 thr_degree,
                                 thr_corr){
  names_of_groups <- names(gCOR.groups)
  names(names_of_groups) <- names_of_groups
  g.interactome <- graph_from_edgelist(as.matrix(interactome[,3:4]), 
                                       directed = F)
  gcor.w <- lapply(gCOR.groups, function(gcor.el){
    
    if (is.null(correctionCorr)){
      gcor.el$p.adj <- gcor.el$p
    } else {
      gcor.el$p.adj <- p.adjust(gcor.el$p, correctionCorr)
    }
    
    ind_sig <- gcor.el$p.adj <= signifCorr
    
    if (compute_weights){
      ### normalizzazione per togliere gli 1 dal pvalue, per rendere possibile il logaritmo
      # L'obiettivo è trasformare l'intervallo del pvalue da (Significativià ; 1 ] a (0 ; 1)
      
      # Il nuovo massimo del pvalue sarà la media tra 1 e il secondo massimo
      if (max(gcor.el$p.adj)==1){
        new_max <- mean(c(1,max(gcor.el$p.adj[-which(gcor.el$p.adj==1)])))
      } else {
        new_max <- max(gcor.el$p.adj)
      }
      
      # calcolo parametri normalizzazione: (pvalue + epsilon ) * k = p_norm1
      # ( signifCorr + eps ) * k = signifCorr
      # ( max(pvalue)  + eps ) * k = new_max
      
      # risolvendo questo sistema si ottiene
      eps <- ( new_max - max(gcor.el$p.adj) ) * signifCorr / ( signifCorr - new_max )
      k   <- ( signifCorr - new_max ) / ( signifCorr - max(gcor.el$p.adj) )
      
      # applichiamo la normalizzazione solo ai pvalue non significativi
      ind_notsig <- which(gcor.el$p.adj > signifCorr)
      norm1 <- (gcor.el$p.adj[ind_notsig] + eps) * k
      
      # adesso applichiamo -log10 per "invertire" i pvalue bassi in pesi alti e viceversa 
      lg    <- -log10(norm1)
      
      ## normalizziamo lg tra 0 e 1, calcolo parametri ( lg + eps_log) * k_log = p_norm2
      # ( MAX_lg  + eps_log ) * k_log = 1
      # ( min(lg) + eps_log ) * k_log = min(lg)
      MAX_lg <- -log10(signifCorr)
      
      # risolvendo il sistema si ottiene 
      eps_log <- ( min(lg) - MAX_lg * min(lg) ) / ( min(lg) - 1 )
      k_log   <- ( min(lg) - 1 ) / ( min(lg) - MAX_lg)
      norm2   <- (lg + eps_log) * k_log
      
      #moltiplichiamo la correlazione per il pnorm3
      gcor.el$weights <- gcor.el$cor
      gcor.el$weights[ind_notsig] <- gcor.el$cor[ind_notsig]*norm2^3
      
    } else {
      gcor.el$weights <- gcor.el$cor
      gcor.el$weights[!ind_sig] <- 0
    }
    
    Sign <- gcor.el[ind_sig,] #storing significative ones
    Sign <- Sign[order(Sign$cor_features),]
    
    return(list(All_correlation  = gcor.el,
                Sign_correlation = Sign))
  })
  not_sign_graphs <- lapply(gcor.w,function(x){
    G <- graph_from_edgelist(as.matrix(x$All_correlation)[,1:2],directed = F)
    edge.attributes(G)$weights <- x$All_correlation$weights
    return(G)
  })
  sign_graphs <- lapply(gcor.w,function(x){
    G <- graph_from_edgelist(as.matrix(x$Sign_correlation)[,1:2],directed = F)
    edge.attributes(G)$weights <- x$Sign_correlation$weights
    return(G)
  })
  
  hubs.sign.percentile  <- lapply(sign_graphs,get.hubs, thr_degree, thr_corr)
  hubs.sign2.percentile <- hubs.diff(hubs.sign.percentile, names_of_groups)
  
  
  hubs.sign.av_degree  <- lapply(sign_graphs,function(x){
    D_x <- igraph::degree(x)
    h <- cbind(names(which(D_x > mean(D_x))),
               D_x[which(D_x > mean(D_x))])
    colnames(h) <- c("nodes","degree")
    return(h)
  })
  hubs.sign2.av_degree <- hubs.diff(hubs.sign.av_degree, names_of_groups)
  
  
  ppi_cor.groups <- lapply(not_sign_graphs,function(x){
    G.x <- intersection(x, g.interactome, keep.all.vertices = F)
    G.y <- delete_edges(G.x, which(edge.attributes(G.x)$weights==0))
    edge.attributes(G.y)$weights[which(edge.attributes(G.y)$weights==1)] <- 
      0.5 + unique(sort(edge.attributes(G.y)$weights))[2]/2 ## al posto di 1 sostituisco la via di mezzo tra 1 e il secondo massimo
    edge.attributes(G.y)$weights[which(edge.attributes(G.y)$weights==-1)] <- 
      -0.5 - unique(sort(edge.attributes(G.y)$weights))[2]/2 ## al posto di 1 sostituisco la via di mezzo tra 1 e il secondo massimo
    return(G.y)
  })
  bet.ppi <- sapply(ppi_cor.groups,function(x){
    betweenness(x,
                directed = F,
                weights = 1-abs(edge.attributes(x)$weights),
                normalized = T)
  })
  
  bottneck <- apply(bet.ppi,2,function(x){
    x[which(x >= mean(x))] <- x[which(x > mean(x))]
    x[which(x < mean(x))]  <- NA
    return(x)
  })
  bottneck <- bottneck[which(rowSums(bottneck,na.rm = T)> 0),]
  
  bottneck2 <- bott.diff(bottneck,names_of_groups)
  
  ed.bet.ppi <- sapply(ppi_cor.groups,function(x){
    ed.bet <- edge_betweenness(x,
                               directed = F,
                               weights = 1 - abs(edge.attributes(x)$weights))
    names(ed.bet) <- apply(as_edgelist(x),1,str_c, collapse = " ")
    return(ed.bet)
  })
  
  ed.bottneck.list <- sapply(ed.bet.ppi,function(x){
    return(x[which(x >= mean(x))])
  })
  
  n.row <- names(ed.bottneck.list[[1]])
  for (i in 1:length(ed.bottneck.list)){
    n.row <- base::union(n.row,names(ed.bottneck.list[[i]]))
  }
  
  ed.bottneck <- matrix(NA, nrow = length(n.row), ncol = length(ed.bottneck.list),
                        dimnames = list(n.row,names(ed.bottneck.list)))
  for (i in colnames(ed.bottneck)){
    ed.bottneck[names(ed.bottneck.list[[i]]),i] <- ed.bottneck.list[[i]]
  }
  ed.bottneck2 <- bott.diff(ed.bottneck,names_of_groups)
  
  return(list(PPI_weighted = ppi_cor.groups,
              Hubs.percentile = list(Hubs = hubs.sign.percentile,
                                     Hubs_Diff = hubs.sign2.percentile),
              Hubs.av_degree = list(Hubs = hubs.sign.av_degree,
                                    Hubs_Diff = hubs.sign2.av_degree),
              Bet.ppi = bet.ppi,
              Bottleneck.av_betw = list(Bottleneck = bottneck,
                                        Bottleneck_Diff = bottneck2),
              Edge.bet.ppi = ed.bet.ppi,
              Edge.Bottleneck.av_edbetw = list(Edge.Bottleneck = ed.bottneck,
                                               Edge.Bottleneck_Diff = ed.bottneck2)))
}
