intersect_graph_list <- function(graph_lists,g1){
  gx <- lapply(X = graph_lists, FUN = intersection, g1, keep.all.vertices = F)
  gx <- lapply(gx, function(x) delete.vertices(x,which(degree(x)==0)))
  names(gx) <- names(graph_lists)
  return(gx)
}
