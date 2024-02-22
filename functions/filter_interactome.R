filter_interactome <- function(interactome, score_type, scores){
  
  score.selected <- interactome[,score_type]
  score.reached <- matrix(nrow=nrow(score.selected),
                          ncol=ncol(score.selected))
  for (i in 1:length(scores)){
    score.reached[,i] <- score.selected[,i] > scores[i]
  }
  
  return(interactome[rowSums(score.reached)>0,])
}
