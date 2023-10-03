setCytoStyle <- function(SUID_list = getNetworkList(getSUIDs = T),
                         style.name = "score",
                         palette.nodes = c('#CCFFFF','#FC0000')){
  all_Suid.collection <- unique(sapply(SUID_list,getCollectionSuid))
  default_param <- sapply(getVisualPropertyNames(),getVisualPropertyDefault,"default")
  nodeLabels <- mapVisualProperty('node label','id','p')
  nodeFills <- mapVisualProperty(visual.prop = 'node fill color',table.column = 'scores',mapping.type = 'c',
                                   table.column.values = as.numeric(c(0,1)),
                                   visual.prop.values = palette.nodes)
  createVisualStyle(style.name, default_param, list(nodeLabels,nodeFills))
  lockNodeDimensions(FALSE, style.name)
  lapply(SUID_list,function(x){
    message(getNetworkName(x), " - ", getCollectionName(getCollectionSuid(x)))
    setVisualStyle(style.name = style.name,network = x)
  })
}
