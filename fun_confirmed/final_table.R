final_table <- function(enr,Gpath.group,tb_mncor,n_groups,names_of_groups){
  
  sp.all <- sort(unique(unlist(lapply(enr,function(x){
    if ( !is.null(dim(x)) )  {
      x[,1]}}))))
  
  if (length(sp.all)){
    # le colonne sono nomi, nodi, edges, media pesi
    tab_s <- matrix(0,nrow=length(sp.all),ncol=3+n_groups*3)
    colnames(tab_s) <- 1:ncol(tab_s)
    # nome path
    rownames(tab_s) <- sp.all
    tab_s[,1] <- sp.all
    
    colnames(tab_s)[1:3] <- c("name_pathways","n_nodes","n_edges")
    
    # n_nodi
    n.1 <- unlist(lapply(Gpath.group[[1]][sp.all],length))
    tab_s[sp.all,2] <- n.1
    
    # n_edge
    e.1 <- unlist(lapply(Gpath.group[[1]][sp.all],function(x){
      length(E(x))
    }))
    tab_s[sp.all,3] <- e.1
    
    n <- 3
    # per ogni gruppo
    for (i in 1:n_groups){
      # nodi
      tab_s[sp.all,n+i]   <- unlist(lapply(Gpath.group[[i+n_groups]][sp.all],length))
      # edges
      tab_s[sp.all,n+n_groups+i] <- unlist(lapply(Gpath.group[[i+n_groups]][sp.all],function(x){
        length(E(x))
      }))
      
      # mean cor
      tab_s[sp.all,n+2*n_groups+i]  <- tb_mncor[sp.all,1+i]
      
      colnames(tab_s)[c(n+i,n+n_groups+i,n+2*n_groups+i)] <- 
        str_c(c("nodes","edges","mean_cor"),colnames(tb_mncor)[i+1],sep = "_")
      
    }
    
    tab_s <- as.data.frame(tab_s)
    tab_s[,2:ncol(tab_s)] <- apply(tab_s[,2:ncol(tab_s)],2,as.numeric)
    
    tab_s_perc <- matrix(nrow=nrow(tab_s),ncol = n_groups,
                         dimnames = list(rownames(tab_s),str_c("perc_edges.",names_of_groups)))
    for (i in 1:n_groups){
      tab_s_perc[,i] <- tab_s[,i+n_groups+3]/tab_s[,3]
    }
    tab_s <- cbind(tab_s,tab_s_perc)
    return(list(tab_s = tab_s,
                tab_s_perc = tab_s_perc))
  } else {
    return(list(tab_s = NULL,
                tab_s_perc = NULL))
  }
  
}
