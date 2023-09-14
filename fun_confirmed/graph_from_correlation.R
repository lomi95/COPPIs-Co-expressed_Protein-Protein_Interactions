graph_from_correlation <- function(cor_group,weights){
  str_group   <- str_split(cor_group$cor_features,pattern = " and ", simplify = T)
  graph_group <- graph_from_edgelist(cbind(str_group[,1],str_group[,2]), directed = F)
  edge.attributes(graph_group)$weights <- weights
  return(graph_group)
}
