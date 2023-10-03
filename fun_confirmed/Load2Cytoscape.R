Load2Cytoscape <- function(file_path, filter_similarity = 0.6){
  all_sheets <- excel_sheets(file_path)
  names(all_sheets) <- all_sheets
  max_score <- 0
  min_score <- Inf
  for (i in all_sheets[-1]){
    message("__________________________________________\n",file_path," - ",i)
    xl.i <- read.xlsx(file_path,sheet = i)
    if (xl.i[1,1]==""){
      next
    } else {
      if (ncol(xl.i)>3){
        nodes <- unique(data.frame(id = xl.i[,1], 
                                   scores = as.numeric(xl.i[,5]),
                                   cluster = as.factor(xl.i[,8]),
                                   tag = as.factor(xl.i[,10])),
                        stringsAsFactors=FALSE)
      } else {
        nodes <- unique(data.frame(id = xl.i[,1],
                                   stringsAsFactors=FALSE))
      }
      max_score <- max(max_score,max(nodes$scores))
      min_score <- min(min_score,min(nodes$scores))
      edges <- data.frame(source=xl.i[,1],
                          target=xl.i[,2],
                          weight=as.numeric(xl.i[,3]),
                          stringsAsFactors=FALSE)
      createNetworkFromDataFrames(nodes,edges[edges$weight > filter_similarity,], 
                                  title= i, collection=str_remove(file_path,".xlsx"))
      loadTableData(read.xlsx(file_path,all_sheets[1]), data.key.column = "name_pathways")
    }
  }
  return(c(min_score,max_score))
}
