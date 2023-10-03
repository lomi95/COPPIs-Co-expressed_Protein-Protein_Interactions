gpp_fun <- function(W,S){
  #weights
  #scores
  W1 <- as.matrix(W) %*% (t(as.matrix(W)))
  W2 <- (W1/rowMaxs(W1) + t(W1/colMaxs(W1)))/2
  
  
  
  graph_pathpath <- graph_from_adjacency_matrix(W2,mode="undirected",
                                                weighted = T, diag = T)
  vertex.attributes(graph_pathpath)$scores <- as.numeric(S[,3])
  #provo louvain con differenti resolutions
  list_c <- lapply(seq(0.05,2, by = 0.05), function(x){
    return(cluster_louvain(graph_pathpath,  
                           edge.attributes(graph_pathpath)$weigth,x))
  })
  
  # numero cluster
  list_k <- sapply(list_c,function(x){
    # Numero di cluster ottenuti dall'algoritmo di Louvain
    k <- length(unique(membership(x)))
  })
  D <- dist(W2)
  
  stats_clusters <- lapply(list_c,function(x){
    # Numero di cluster ottenuti dall'algoritmo di Louvain
    k <- length(unique(membership(x)))
    if (k>1 & k<length(graph_pathpath)){
      # Calcola gli indici di valutazione
      return(cluster.stats(D, membership(x)))
      
    } else return(0)
    
  })
  
  #minimum average dissimilarity between two cluster / 
  #maximum average within cluster dissimilarity
  dunn2 <- sapply(stats_clusters, function(x) {
    if (is.numeric(x)){ return(0)} else x$dunn2
  })
  
  k.min <- length(graph_pathpath)/40 #40 = numero medio di pathways raggruppati in un cluster
  if (length(intersect(order(dunn2,decreasing = T),which(list_k>=k.min)))) {
    clusters.opt <- list_c[[intersect(order(dunn2,decreasing = T),which(list_k>=k.min))[1]]]
  } else {
    clusters.opt <- list_c[[which.max(dunn2)]]
  }
  
  names_cluster <- list()
  clusters.opt$tag <- vector(length = length(clusters.opt$membership))
  for (i in names(table(clusters.opt$membership))){
    names_vertex <- names(V(graph_pathpath))[which(clusters.opt$membership==i)]
    docs <- Corpus(VectorSource(names_vertex))
    docs <- tm_map(docs, toSpace, ":")
    docs <- tm_map(docs, toSpace, "'")
    docs <- tm_map(docs, toSpace, "'")
    docs <- tm_map(docs, toSpace, "-")
    docs <- tm_map(docs, removePunctuation)
    docs <- tm_map(docs, content_transformer(tolower))
    docs <- tm_map(docs, removeWords, stopwords())
    docs <- tm_map(docs, stripWhitespace)
    docs <- tm_map(docs, stemDocument)
    docs <- tm_map(docs, removeWords, c("process","regul","negat","posit","pathway"))
    dtm <- DocumentTermMatrix(docs)
    
    freq <- colSums(as.matrix(dtm))
    ord  <- sort(freq,decreasing = T)
    
    if (length(docs)<15){
      freq.proc <- names(ord[1:2])
    } else if (sum(ord>(length(docs)/4))){
      freq.proc <- c(names(which(ord>(length(docs)/4))),
                     setdiff(names(which(ord>(length(docs)/5))),
                             names(which(ord>(length(docs)/4)))))
    } else {
      freq.proc <- names(head(ord,round(length(ord)/10)))
    }
    
    
    all_words <- tolower(unlist(str_split(str_c(
      str_replace(str_remove(names_vertex,"/"),"-"," "),
      collapse = " ")," ")))
    freq.proc.stem.compl <- sapply(freq.proc,stemCompletion,all_words)
    
    names_cluster[[as.numeric(i)]] <- freq.proc.stem.compl
    names(names_cluster)[as.numeric(i)] <- stemCompletion(names(ord)[1],all_words)
    
    clusters.opt$tag[which(clusters.opt$membership==i)] <- str_c(freq.proc.stem.compl, collapse = ", ")
    clusters.opt$membership[which(clusters.opt$membership==i)] <- stemCompletion(names(ord)[1],all_words)
  }
  graph_pathpath <- set_vertex_attr(graph = graph_pathpath,
                                    name = "tag",
                                    value = clusters.opt$tag)
  graph_pathpath <- set_vertex_attr(graph = graph_pathpath,
                                    name = "cluster",
                                    value = clusters.opt$membership)
  
  el_gpp <- cbind(as_edgelist(graph_pathpath),
                  edge.attributes(graph_pathpath)$weight)
  
  # aggiungere gli scores e i pvalue
  el_gpp <- cbind(el_gpp,S[match(el_gpp[,1],S[,1]),2:3],S[match(el_gpp[,2],S[,1]),2:3],
                  get.vertex.attribute(graph_pathpath,"cluster",el_gpp[,1]),
                  get.vertex.attribute(graph_pathpath,"cluster",el_gpp[,2]),
                  get.vertex.attribute(graph_pathpath,"tag",el_gpp[,1]),
                  get.vertex.attribute(graph_pathpath,"tag",el_gpp[,2]))
  
  colnames(el_gpp) <- c("Node_path1","Node_path2",
                        "Similarity",
                        "p.value1","score1",
                        "p.value2","score2",
                        "cluster1","cluster2",
                        "tag1","tag2")
  return(list(gpp = graph_pathpath,
              el_gpp = el_gpp))
}
