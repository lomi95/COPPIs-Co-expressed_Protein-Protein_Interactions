gpp.compare <- function(gpp.i,names_groups, ign_case){
  all_compared    <- list()
  all_compared.el <- list()
  
  n <- 1
  for (i in names_groups){
    ind_gpp <- grep(str_c("diviso ",i),names(gpp.i),ignore.case = ign_case)
    l_gpp <- lapply(gpp.i,function(x) x$gpp)
    inters <- l_gpp[[ind_gpp[1]]]
    for (j in 2:length(ind_gpp)){
      inters <- intersection(inters, l_gpp[[ind_gpp[j]]],keep.all.vertices = F)
    }
    all_compared[[n]] <- inters
    if (length(inters)>0){
      all_compared.el[[n]] <- cbind(as_edgelist(inters),edge.attributes(inters)[[1]])
    } else {
      all_compared.el[[n]] <- matrix(nrow = 0, ncol=3)
    }
    names(all_compared)[n] <- str_c("all diviso ", i)
    names(all_compared.el)[n] <- str_c("all diviso ", i)
    
    n <- n+1
    ind_gpp <- grep(str_c(i," diviso"),names(gpp.i),ignore.case = ign_case)
    l_gpp <- lapply(gpp.i,function(x) x$gpp)
    inters <- l_gpp[[ind_gpp[1]]]
    for (j in 2:length(ind_gpp)){
      inters <- intersection(inters, l_gpp[[ind_gpp[j]]], keep.all.vertices = F)
    }
    all_compared[[n]] <- inters
    names(all_compared)[n] <- str_c(i, " diviso all")
    if (length(inters)>0){
      all_compared.el[[n]] <- cbind(as_edgelist(inters),edge.attributes(inters)[[1]])
    } else {
      all_compared.el[[n]] <- matrix(nrow = 0, ncol=3)
    }
    names(all_compared.el)[n] <- str_c(i," diviso all")
    n <- n+1
  }
  return(list(all = all_compared,
              all_el = all_compared.el))
}
