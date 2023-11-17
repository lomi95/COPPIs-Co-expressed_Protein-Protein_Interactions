Compute_centralities <- function(cor_groups,
                                 g.interactome,
                                 names_of_groups,
                                 thr_degree,
                                 thr_corr,
                                 ignore_case){
  not_sign_graphs <- cor_groups$cor_graphs[1:(length(cor_groups$cor_graphs)/2)]
  sign_graphs     <- cor_groups$cor_graphs[((length(cor_groups$cor_graphs)/2+1):
                                              length(cor_groups$cor_graphs))]
  names(sign_graphs)     <- str_remove(names(sign_graphs), "Sign_")
  names(not_sign_graphs) <- str_remove(names(not_sign_graphs), "Notsign_")
  
  hubs.sign.percentile  <- lapply(sign_graphs,get.hubs, thr_degree, thr_corr)
  hubs.sign2.percentile <- hubs.diff(hubs.sign.percentile, names_of_groups,ignore_case)
  
  
  hubs.sign.av_degree  <- lapply(sign_graphs,function(x){
    D_x <- igraph::degree(x)
    h <- cbind(names(which(D_x > mean(D_x))),
               D_x[which(D_x > mean(D_x))])
    colnames(h) <- c("nodes","degree")
    return(h)
  })
  hubs.sign2.av_degree <- hubs.diff(hubs.sign.av_degree, names_of_groups,ignore_case)
  
  
  ppi_cor.groups <- lapply(not_sign_graphs,function(x){
    G.x <- intersection(x, g.interactome, keep.all.vertices = F)
    G.y <- delete.edges(G.x, which(edge.attributes(G.x)$weights==0))
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
  
  bottneck2 <- bott.diff(bottneck,names_of_groups,ignore_case)
  
  ed.bet.ppi <- sapply(ppi_cor.groups,function(x){
    ed.bet <- edge.betweenness(x,
                               directed = F,
                               weights = 1 - abs(edge.attributes(x)$weights))
    names(ed.bet) <- apply(as_edgelist(x),1,str_c, collapse = " ")
    return(ed.bet)
  })
  
  ed.bottneck.list <- sapply(ed.bet.ppi,function(x){
    return(x[which(x >= mean(x))])
  })
  
  n.row <- names(ed.bottneck.list$dada)
  for (i in 1:length(ed.bottneck.list)){
    n.row <- base::union(n.row,names(ed.bottneck.list[[i]]))
  }
  
  ed.bottneck <- matrix(NA, nrow = length(n.row), ncol = length(ed.bottneck.list),
                        dimnames = list(n.row,names(ed.bottneck.list)))
  for (i in colnames(ed.bottneck)){
    ed.bottneck[names(ed.bottneck.list[[i]]),i] <- ed.bottneck.list[[i]]
  }
  ed.bottneck2 <- bott.diff(ed.bottneck,names_of_groups,ignore_case)
  
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
