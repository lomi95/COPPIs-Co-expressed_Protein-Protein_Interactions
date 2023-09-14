gpp_fun <- function(W,S,filtering = T,
                    node.comp_thr = node.comp_thr, 
                    quantile.node_thr = quantile.node_thr,
                    quantile.edge_thr = quantile.edge_thr){
  #weights
  #scores
  W1 <- as.matrix(W) %*% (t(as.matrix(W)))
  W2 <- (W1/rowMaxs(W1) + t(W1/colMaxs(W1)))/2
  
  
  
  graph_pathpath <- graph_from_adjacency_matrix(W2,mode="undirected",
                                                weighted = T, diag = T)
  vertex.attributes(graph_pathpath)$scores <- as.numeric(S[,3])
  if (filtering){
    graph_pathpath <- filter_network(graph_pathpath,
                                     node.comp_thr = node.comp_thr, 
                                     quantile.node_thr = quantile.node_thr,
                                     quantile.edge_thr = quantile.edge_thr)
  }
  
  
  el_gpp <- cbind(as_edgelist(graph_pathpath),
                  edge.attributes(graph_pathpath)$weight)
  
  # aggiungere gli scores e i pvalue
  el_gpp <- cbind(el_gpp,S[match(el_gpp[,1],S[,1]),2:3],S[match(el_gpp[,2],S[,1]),2:3])
  colnames(el_gpp) <- c("Node_path1","Node_path2",
                        "Similarity",
                        "p.value1","score1",
                        "p.value2","score2")
  return(list(gpp = graph_pathpath,
              el_gpp = el_gpp))
}


