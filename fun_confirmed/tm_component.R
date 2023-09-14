tm_components <- function(netw, min_comp = 4,
                          freq.thr = 2.5,
                          freq.thr1 = 4.5){
  print(netw)
  g.tm <- createIgraphFromNetwork(netw)
  c.g <- components(g.tm)
  groups <- which(c.g$csize>=min_comp)
  
  if (length(groups)>0){
    for (i in groups){
      m.g <- names(which(c.g$membership==i))
      corpus <- Corpus(VectorSource(m.g))
      corpus <- tm_map(corpus, toSpace, ":")
      corpus <- tm_map(corpus, toSpace, "'")
      corpus <- tm_map(corpus, toSpace, "/")
      corpus <- tm_map(corpus, toSpace, "-")
      corpus <- tm_map(corpus, removePunctuation)
      corpus <- tm_map(corpus, content_transformer(tolower))
      corpus <- tm_map(corpus, removeWords, stopwords())
      corpus <- tm_map(corpus, stripWhitespace)
      corpus <- corpus %>% tm_map(stemDocument)
      dtm <- DocumentTermMatrix(corpus)
      freq.terms <- colSums(as.matrix(dtm))
      freq.terms.1 <- names(which(freq.terms> length(freq.terms)/freq.thr))
      freq.terms.2 <- names(which(freq.terms> length(freq.terms)/freq.thr1))
      
      if (length(freq.terms.2)){
        
        pos_nodes <- apply(getNodePosition(selectNodes(m.g, 'id',
                                                       network = netw,
                                                       preserve.current.selection = F)$nodes,network = netw),2,as.numeric)
        pos_mean <- colMeans(as.matrix(pos_nodes))
        
        if (length(freq.terms.1)){
          addAnnotationText(network = netw,
                            text = str_flatten_comma(freq.terms.1), 
                            x.pos = pos_mean[1]-nchar(str_flatten_comma(freq.terms.1))/0.5, 
                            y.pos = pos_mean[2]-11,
                            fontSize = 20)
          addAnnotationText(network = netw,
                            text = str_flatten_comma(setdiff(freq.terms.2,freq.terms.1)), 
                            x.pos = pos_mean[1]-nchar(str_flatten_comma(setdiff(freq.terms.2,freq.terms.1))),
                            y.pos = pos_mean[2]+11,
                            fontSize = 15)
        } else {
          addAnnotationText(network = netw,
                            text = str_flatten_comma(setdiff(freq.terms.2,freq.terms.1)), 
                            x.pos = pos_mean[1]-nchar(str_flatten_comma(setdiff(freq.terms.2,freq.terms.1))),
                            y.pos = pos_mean[2],
                            fontSize = 17)
        }
      }
    }
  }
}

