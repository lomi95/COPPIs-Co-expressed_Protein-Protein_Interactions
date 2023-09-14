graph_ideal_fun <- function(col_i, tab_proteinpathways){
  ind_p            <- rownames(tab_proteinpathways)[which(tab_proteinpathways[,col_i]==1)]
  m_temp           <- matrix(1,ncol=length(ind_p),nrow = length(ind_p))
  rownames(m_temp) <- ind_p
  colnames(m_temp) <- ind_p
  g_i              <- graph_from_adjacency_matrix(m_temp, diag=F,mode = "undirected")
  return(g_i)
}
