filter.edge.pipe <- function(coppi.res, 
                             threshold, 
                             min_edges = 1, 
                             max_edges = Inf,
                             ignore_case = T){
  sig.all <- coppi.res$SIG.all
  names_of_groups <- colnames(sig.all[[1]]$tab_meancor)[-1]
  
  # filtro numero edges minimo e massimo considerato
  sig.all <- lapply(sig.all,function(x){
    path_filt <- x$tab_summary$name_pathways[intersect(which(x$tab_summary$n_edges > min_edges),
                                                       which(x$tab_summary$n_edges < max_edges))]
    x$sig_path.p <- lapply(x$sig_path.p, function(y){
      return(y[na.omit(match(path_filt,y$Description)),])
    })
    return(x)
  })
  
  # filtro percentuale significativi su numero edges dei pathways 
  sig.new <- filtering_percentage_edge(sig.all, threshold)
  
  start_time <- Sys.time() 
  message("\nCreating pathways-pathways graphs")
  gpp.all <- lapply(sig.new, function(x){
    if (!is.null(x$sig_path.p)){
      gpp.pair <- lapply(x$sig_path.p, function(y){
        if (length(y[,2])>0){
          gpp.pairwise <- gpp_fun(as.data.frame(x$tab_pathprot)[y[,1],],y)
        } else {
          gpp.pairwise <- list(gpp = make_empty_graph(directed = F),
                               el_gpp = matrix(nrow=0,ncol=7))
        }
      })
      if (length(names_of_groups)>2){
        gpp.complete <- gpp.compare(gpp.pair,names_of_groups,ignore_case)
        return(list(pairwise = gpp.pair,
                    complete = gpp.complete))
      } else {
        return(list(pairwise = gpp.pair))
      }
    }
  })
  time_employed <- format_time_difference(difftime(Sys.time(),start_time,units = "s"))
  message(str_c("Pathways-pathways graphs created in", time_employed,sep = " "))
  
  coppi.res$SIG.all <- sig.new
  coppi.res$gpp.all <- gpp.all
  return(coppi.res)
}
