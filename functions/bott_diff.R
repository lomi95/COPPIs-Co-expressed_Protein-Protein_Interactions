bott.diff <- function(bott.sign, names_of_groups){
  pairwise.diff <- list()
  n <- 1
  for (i in 1:length(names_of_groups)){
    for (j in 1:length(names_of_groups)){
      if (i == j){
        next
      }
      I <- grep(names_of_groups[i],colnames(bott.sign))
      J <- grep(names_of_groups[j],colnames(bott.sign))
      pairwise.diff[[n]] <- setdiff(names(which(bott.sign[,I]>0)),
                                    names(which(bott.sign[,J]>0)))
      names(pairwise.diff)[n] <- str_c(names_of_groups[i]," UP vs ",names_of_groups[j])
      n <- n+1
    }
  }
  
  
  complete.diff <- list()
  n <- 1
  for (i in names_of_groups){
    if (sum(grepl(str_c(i, " UP vs"), names(pairwise.diff)))){
      list_x <- pairwise.diff[grep(str_c(i, " UP vs"), names(pairwise.diff))]
      complete.diff[[n]] <- list_x[[1]]
      for (j in 1:length(list_x)){
        complete.diff[[n]] <- intersect(complete.diff[[n]], list_x[[j]])
      }
      names(complete.diff)[n] <- str_c(i, " UP vs all")
      n <- n+1
    }
    if (sum(grepl(str_c(" UP vs ",i), names(pairwise.diff)))){
      list_x <- pairwise.diff[grep(str_c(" UP vs ",i), names(pairwise.diff))]
      complete.diff[[n]] <- list_x[[1]]
      for (j in 1:length(list_x)){
        complete.diff[[n]] <- intersect(complete.diff[[n]], list_x[[j]])
      }
      names(complete.diff)[n] <- str_c(i, " DOWN vs all")
      n <- n+1
    }
  }
  
  return(list(Pairwise = pairwise.diff,
              Complete = complete.diff))
}
