get.hubs <- function(G, thr_degree = 0.75, thr_corr = NULL){
  
  if (sum(abs(edge.attributes(G)$weights) < thr_corr)){
    G <- delete.edges(G, which(abs(edge.attributes(G)$weights) < thr_corr))
  }
  G <- delete.vertices(G, which(degree(G)==0))
  
  if (sum(degree(G)<quantile(degree(G),thr_degree))){
    G.1 <- delete.vertices(G, which(degree(G)<quantile(degree(G),thr_degree)))
  } else {
    message("no nodes were filtered")
    G.1  <- G
  }
  n.deg <- cbind(names(V(G.1)),degree(G.1))
  colnames(n.deg) <- c("nodes","degree")
  return(n.deg[order(n.deg[,1]),])
}
