# ho una matrice P in cui le righe (N) rappresentano i processi, le colonne (M) rappresentano le proteine e le singole celle hanno 1 se la proteina è nel pathway, 0 altrimenti
# 
# Ho un'altra matrice C (M X M) che contiene le correlazioni tra le proteine
# 
# Prendiamo ad esempio due pathway casuali ni e nj:
# pij sono le proteine comuni tra i due pathway
# pii sono le proteine solo di ni
# pjj sono le proteine solo di nj
# cii sono le correlazioni delle proteine pii
# cjj sono le correlazioni delle proteine pjj
# cij sono le correlazioni tra pii e pjj
# NE_ii è la grandezza di cii
# NE_jj è la grandezza di cjj
# NE_ij è la grandezza di cij
# Voglio una matrice  W che abbia dimensione (N X N), in cui sia le righe, 
# sia le colonne rappresentano i pathway. le celle di W hanno come valore:
# pij/(pii+pjj+pij) + somma del valore assoluto di cij/(NE_ii + NE_jj +NE_ij)
# c'è un modo per rappresentare W in forma matriciale

# path i e j

fun_Wcor <- function(W,S,cor_groups){
  W.cor1 <- lapply(cor_groups$Rcor,function(x){
    st.time <- Sys.time()
    signif_matrix <- ifelse(x$P < 0.05, x$r, 0)
    W1 <- matrix(nrow=nrow(W),ncol = nrow(W),
                dimnames = list(rownames(W),rownames(W)))
    for (i in 1:(nrow(W1)-1)){
      p.i  <- W[i,]
      
      for (j in (i+1):ncol(W1)){
        p.j  <- W[j,]
        p.ij <- as.numeric(p.i & p.j)
        p.ii <- p.i - p.ij
        p.jj <- p.j - p.ij
        
        
        c.ii <- signif_matrix[names(which(p.ii>0)),names(which(p.ii>0))]
        c.jj <- signif_matrix[names(which(p.jj>0)),names(which(p.jj>0))]
        
        asd <- union(names(which(p.ii>0)), names(which(p.jj>0)))
        c.ij <- signif_matrix[asd,asd]
        
        W1[i,j] <- sum(p.ij)/sum(p.ii+p.jj+p.ij) + sum(abs(c.ij),na.rm = T)/(length(c.ii)+length(c.ij)+length(c.jj))
        
      }
    }
    diag(W1) <- 1
    return(W1)
  })
  
  graph_pathpath <- graph_from_adjacency_matrix(W.cor1,mode="undirected",
                                                weighted = T, diag = T)
  el_gpp <- cbind(as_edgelist(graph_pathpath),
                  edge.attributes(graph_pathpath)$weight)
  
  # aggiungere gli scores
  el_gpp <- cbind(el_gpp,S[match(el_gpp[,1],S[,2]),1],S[match(el_gpp[,2],S[,2]),1])
  
  return(list(gpp = graph_pathpath,
              el_gpp = el_gpp))
  
}




