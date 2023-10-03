tm_components <- function(SUID_list = getNetworkList(getSUIDs = T)){
  lapply(SUID_list, function(netw) {
    message(getNetworkName(netw)," - ", 
            getCollectionName(getCollectionSuid(netw)))
    layoutNetwork(network = netw, layout.name = "attributes-layout nodeAttribute=tag")
    node.tags <- getTableColumns(network = netw,columns = "tag")
    
    for (i in names(table(node.tags$tag))){
      xy <- getNodePosition(network = netw,rownames(node.tags)[which(node.tags$tag==i)])
      if (is.null(nrow(xy))){
        pos <- xy
        pos[2] <- xy[2] - 100
      } else {
        pos <- colMeans(xy)
        pos[2] <- colMins(as.matrix(xy))[2] - 100
      }
      addAnnotationText(network = netw,text = i,x.pos = pos[1]-nchar(i)*7,y.pos = pos[2],
                        fontSize = 30)
    }
  })
}
