filtering_percentage_edge <- function(sig.all, threshold){
  new.sig <- lapply(sig.all,function(x){
    n_groups <- (ncol(x$tab_summary)-8)/4
    N.edges.thr <- apply(x$tab_summary[,(3+n_groups*3+1):(3+n_groups*4)],2,function(y){
      return(y>threshold)
    })
    for (i in 1:ncol(N.edges.thr)){
      ind_sigpath <- grep(str_c("diviso ",colnames(x$tab_meancor)[i+1]),names(x$sig_path.p))
      for (j in ind_sigpath){
        x$sig_path.p[[j]] <- x$sig_path.p[[j]][na.omit(match(names(which(N.edges.thr[,i])),x$sig_path.p[[j]]$Description)),]
      }
    }
    return(x)
  })
  return(new.sig)
}
