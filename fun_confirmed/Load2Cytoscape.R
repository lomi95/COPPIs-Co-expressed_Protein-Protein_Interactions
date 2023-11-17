Load2Cytoscape <- function(file_path, filter_similarity = 0.6, 
                           filter_score = NULL,
                           filter_first_x = NULL){
  all_sheets <- getSheetNames(file_path)
  names(all_sheets) <- all_sheets
  for (i in all_sheets[-1]){
    message("__________________________________________\n",file_path," - ",i)
    xl.i <- read.xlsx(file_path,sheet = i)
    if (ncol(xl.i)>3){
      if (sum(as.numeric(xl.i[,5]) >= filter_score)){
        xl.i <- xl.i[as.numeric(xl.i[,5]) >= filter_score,]
      }
      if (sum(filter_first_x)){
        filt.rows <- min(c(length(unique(xl.i[,1])),filter_first_x))
        min_score <- sort(as.numeric(xl.i[which(xl.i[,1] == xl.i[,2]),5]),decreasing = T)[filt.rows]
        
        xl.i <- xl.i[which(as.numeric(xl.i[,5]) >= min_score),]
        xl.i <- xl.i[which(as.numeric(xl.i[,7]) >= min_score),]
        
      }
    }
    
    
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
      edges <- data.frame(source=xl.i[,1],
                          target=xl.i[,2],
                          weight=as.numeric(xl.i[,3]),
                          stringsAsFactors=FALSE)
      createNetworkFromDataFrames(nodes,edges[which(edges$weight > filter_similarity),], 
                                  title= i, collection=str_remove(file_path,".xlsx"))
      loadTableData(read.xlsx(file_path,all_sheets[1]), data.key.column = "name_pathways")
    }
  }
}
