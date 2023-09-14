pipe_for_term <- function(enr_sample,
                          cor_groups1,
                          sd_groups,
                          graph_cor.groups,
                          significance,
                          p.adj,
                          interac,
                          min_edge,
                          max_edge,
                          filter_percentage){
  
  start_time <- Sys.time()
  
  proteins <- unique(unlist(enr_sample$inputGenes))
  pathways <- enr_sample$description
  
  tab_pathways_protein <- matrix(ncol=length(proteins),nrow = length(pathways),data = 0,
                                 dimnames = list(pathways,proteins))
  
  for (i in 1:nrow(enr_sample)){
    tab_pathways_protein[i,match(unlist(enr_sample[i,]$inputGenes),proteins)] <- 1
  }
  
  tab_pathways_protein <- aggregate(tab_pathways_protein[],
                                    list(rownames(tab_pathways_protein)), 
                                    FUN = max, na.rm = TRUE)
  rownames(tab_pathways_protein) <- tab_pathways_protein$Group.1
  tab_pathways_protein <- as.matrix(tab_pathways_protein[,-1])
  #creating path-path network
  g.interac <- graph_from_edgelist(as.matrix(interac[,3:4]),directed = F)
  
  
  ## intersezione dell'interattoma con le interazioni INTRA-pathway
  
  # step 1: creare un grafo per ogni pathway, immaginando che tutte le proteine 
  # interagiscano tra loro
  
  # step 2: intersecare con il grafo dell'interattoma, considerando solo le interazioni
  # reali
  
  # step 3: unire tutti questi grafi per avere il grafo completo intersecato con
  # l'interattoma
  tab_protein_pathways <- t(tab_pathways_protein)
  #step 1
  Graph_pathways_ideal <- lapply(1:ncol(tab_protein_pathways),graph_ideal_fun, tab_protein_pathways)  
  names(Graph_pathways_ideal) <- colnames(tab_protein_pathways)
  
  #step 2
  Graph_pathways.1 <- lapply(Graph_pathways_ideal,intersection, ... = g.interac, keep.all.vertices = F)
  names(Graph_pathways.1) <- colnames(tab_protein_pathways)
  N_edges.1 <- unlist(lapply(Graph_pathways.1, function(x){length(E(x))}))
  Graph_pathways <- Graph_pathways.1[intersect(which(N_edges.1>=min_edge),
                                               which(N_edges.1<=max_edge))]
  
  #step 3
  Graph.all_pathways <- graph.union(Graph_pathways)
  
  # intersezione delle interazioni con le correlazioni
  graph_inter.groups <- lapply(graph_cor.groups, intersection, g.interac, keep.all.vertices = F)
  
  graph_allpathways.groups <- lapply(graph_inter.groups, intersection, Graph.all_pathways, keep.all.vertices = F)
  
  Graph_pathways.groups <- list()
  for (i in 1:length(graph_allpathways.groups)){
    Graph_pathways.groups[[i]] <- intersect_graph_list(Graph_pathways,graph_allpathways.groups[[i]])
  }
  names(Graph_pathways.groups) <- names(graph_allpathways.groups)
  
  
  N.g <- length(Graph_pathways.groups)/2
  N_edges_sign <- lapply(Graph_pathways.groups[(N.g+1):(2*N.g)],function(x) {
    sapply(x,function(y) length(E(y)))
  })
  
  
  # calcoliamo il rapporto tra la media dei pesi dei gruppi (k)
  tab_k <- matrix(nrow=length(cor_groups1),ncol=length(cor_groups1))
  for (i in 1:nrow(tab_k)){
    for (j in 1:ncol(tab_k)){
      tab_k[i,j] <- mean(abs(cor_groups1[[i]]$ALL$weights))/
        mean(abs(cor_groups1[[j]]$ALL$weights))
    }
  }
  colnames(tab_k) <- names(cor_groups1)
  rownames(tab_k) <- names(cor_groups1)
  
  
  mean_abs_cor <- function(Graph_pathways.groups,g1){
    asd <- unlist(lapply(Graph_pathways.groups[[g1]],function(x){mean(abs(edge.attributes(x)$weights))}))
  }
  
  tab_meancor <- matrix(nrow=length(Graph_pathways.groups[[1]]), ncol=length(cor_groups1)+1)
  for (i in 1:length(cor_groups1)){
    tab_meancor[,1] <- unlist(lapply(Graph_pathways.groups[[i]],function(x){length(E(x))}))
    tab_meancor[,i+1] <- mean_abs_cor(Graph_pathways.groups,i)
  }
  rownames(tab_meancor) <- names(Graph_pathways.groups[[1]])
  colnames(tab_meancor) <- c("N_edges",names(cor_groups1))
  tab_meancor <- tab_meancor[order(tab_meancor[,1]),]
  
  tab_meancor.sig <- tab_meancor
  for (i in 1:length(N_edges_sign)){
    tab_meancor.sig[names(which(N_edges_sign[[i]]==0)),i+1] <- 0
  }
  ## filter by edges e tolgo i path non significativi in tutti i gruppi
  if (sum(rowSums(tab_meancor.sig[,2:ncol(tab_meancor.sig)])==0)){
    tab_meancor.sig.1 <- tab_meancor.sig[-which(rowSums(tab_meancor.sig[,2:ncol(tab_meancor.sig)])==0),]
  } else {
    tab_meancor.sig.1 <- tab_meancor.sig
  }
  
  tab_meancor.fe <- tab_meancor.sig.1[intersect(which(tab_meancor.sig.1[,1]>=min_edge),
                                                which(tab_meancor.sig.1[,1]<=max_edge)),]
  
  
  # funzione che prende i significativi e assegna lo score
  path_sig <- function(tb_k,
                       g1,g2,
                       tb_meancor,
                       p.adj){
    g2_0 <- which(tb_meancor[,g2+1]==0)
    tb_meancor[g2_0,] <- 1000
    if (tb_k[g1,g2] > 1){
      k_x <- 1
    } else {
      k_x <- tb_k[g1,g2]
    }
    
    ratiog1g2 <- cbind(tb_meancor[,1],tb_meancor[,g1+1]/tb_meancor[,g2+1])
    p.v <- apply(ratiog1g2,1, function(x){
      return(pnorm(x[2],mean = k_x, sqrt((sd_groups[g1]^2+sd_groups[g2]^2)/x[1])))
    })
    
    p.v.adj <- p.adjust(p.v,p.adj)
    p.sign <- p.v.adj[intersect(which(p.v.adj < significance),which(ratiog1g2[,2] < k_x))]
    
    if (length(p.sign)>0){
      
      scores <- matrix(nrow = length(p.sign), ncol = 3)
      n <- 1
      for (i in names(p.sign)){
        
        # trovare la distanza minima tra il punto e la curva 
        
        
        y0 <- ratiog1g2[i,2]
        x0 <- ratiog1g2[i,1]
        d0 <- abs(k_x - 2*sqrt((sd_groups[g1]^2+sd_groups[g2]^2)/x0) - y0)
        
        if (x0 == 1){
          d.i <- d0
        } else {
          d.i <- d0 + 1
          x.i <- x0
          while (x.i != 1){
            x.i <- x.i - 1
            y.i <- abs(k_x - 2*sqrt((sd_groups[g1]^2+sd_groups[g2]^2)/x.i) - y0)
            d.i <- sqrt((x.i-x0)^2+y.i^2)
            if (d.i > d0) {
              break
            } else {
              d0 <- d.i
            }
          }
        }
        if (x0 < d.i){
          d0 <- x0
        }
        
        scores[n,1] <- names(p.sign)[n]
        scores[n,2] <- p.sign[n]
        scores[n,3] <- d0* (k_x-ratiog1g2[i,2])/k_x
        
        n <- n+1
      }
      rownames(scores) <- p.sign
      scores <- as.data.frame(scores)
      scores[,3] <- as.numeric(scores[,3])
      scores[,2] <- as.numeric(scores[,2])
      return(scores)
    }
  }
  
  
  sig_path <- list()
  
  n <-1
  for (i in 1:length(cor_groups1)){
    for (j in 1:length(cor_groups1)){
      if(i == j){
        next
      }
      
      ps <- path_sig(tb_k = tab_k,
                     g1 = i,g2 = j,
                     tb_meancor = tab_meancor.fe,
                     p.adj = p.adj)
      if (!is.null(ps)){
        sig_path[[n]] <- ps
        
      } else {
        sig_path[[n]] <- as.data.frame(matrix(nrow=0,ncol = 3))
      }
      
      names(sig_path)[n] <- str_c(names(cor_groups1)[i]," diviso ",names(cor_groups1)[j])
      n <- n+1
    }
  }
  
  tab_X <- final_table(sig_path,Graph_pathways.groups,tab_meancor.fe,N.g,names(cor_groups1))
  sig_path.p <- sig_path
  
  for (i in 1:N.g){
    ind_s <- grep(str_c("diviso ",names(cor_groups1)[i]), names(sig_path.p), ignore.case = T)
    sig_path.p[ind_s] <- lapply(sig_path.p[ind_s],function(x){
      filt <- x[which(tab_X$tab_s_perc[x[,1],i]>filter_percentage),]
    })
  }
  
  
  time_employed <- round(Sys.time()- start_time,digits = 1)
  print(str_c(enr_sample$category[1],"finished in", time_employed,sep = " "))
  
  return(list(sig_path = sig_path,
              sig_path.p = sig_path.p,
              tab_meancor = tab_meancor.fe,
              tab_pathprot = tab_pathways_protein,
              Graph_everypath.group = Graph_pathways.groups,
              tab_summary =tab_X$tab_s))
  
}
