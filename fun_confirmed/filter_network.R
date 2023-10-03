filter_network <- function(netw,node.comp_thr = 25, 
                           quantile.node_thr = 0.15,
                           quantile.edge_thr = 0.33){
  netw.rmnodes <- delete.vertices(netw, vertex.attributes(netw)$scores <=
                                   quantile(vertex.attributes(netw)$scores,quantile.node_thr))
  
  netw.rmedges <- delete.edges(netw.rmnodes, which(edge.attributes(netw.rmnodes)$weight <=
                                                   quantile(edge.attributes(netw.rmnodes)$weight,quantile.edge_thr)))
  c.netw <- components(netw.rmedges)
  
  while(sum(c.netw$csize>node.comp_thr)>0){
    
    n.overcomponent <- which(c.netw$csize>node.comp_thr)
    nodes.comp.i <- which(c.netw$membership==n.overcomponent[1])
    score.component.i <- vertex.attributes(netw.rmedges,nodes.comp.i)$scores
    netw.rmedges <- delete.vertices(netw.rmedges,nodes.comp.i[score.component.i <=
                                                              quantile(score.component.i,
                                                                       quantile.node_thr)])
    if (sum(components(netw.rmedges)$csize>node.comp_thr)==0){
      break
    }
    netw.rmedges <- delete.edges(netw.rmedges, which(edge.attributes(netw.rmedges)$weight <=
                                                     quantile(edge.attributes(netw.rmedges)$weight,
                                                              quantile.edge_thr)))
    
    c.netw <- components(netw.rmedges)
  }
  return(netw.rmedges)
}
