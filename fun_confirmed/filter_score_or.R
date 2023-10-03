filter_score_or <- function(interactome,
                            score_types, 
                            score_values){
  
  if (length(score_types)!=length(score_values)){
    print("length of score_types")
    print(length(score_types))
    print("------------")
    print("length of score_values")
    print(length(score_values))
  }
  
  conditions <- vector()
  
  for (i in 1:length(score_types)){
    
    if (!length(conditions)){
      conditions <- which(interactome[,score_types[i]]>score_values[i])
    }
    conditions <- unique(base::union(conditions,which(interactome[,score_types[i]]>score_values[i])))
  }
  return(interactome[conditions,])
}
