save_excel.centralities <- function(coppi.obj, name_analysis){
  
  Centr.wb <- createWorkbook()
  
  
  addWorksheet(Centr.wb,"Hubs_percentile")
  H <- coppi.obj$Centralities$Hubs.percentile$Hubs
  All.hubs <- unique(unlist(sapply(H, function(x) x[,1])))
  tab.hubs <- matrix(NA,nrow = length(All.hubs),ncol = length(H),
                     dimnames = list(All.hubs, names(H)))
  
  for (i in 1:length(H)){
    tab.hubs[H[[i]][,1],i] <- as.numeric(H[[i]][,2])
  }
  writeDataTable(Centr.wb,"Hubs_percentile",as.data.frame(tab.hubs),
                 rowNames = T)
  
  
  addWorksheet(Centr.wb,"Hubs_percentile_comparison")
  n.row <- max(sapply(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Pairwise,length)) 
  tab.hubs.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Pairwise)+
                                  length(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Complete))
  colnames(tab.hubs.comparison) <- c(names(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Pairwise),
                                     names(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Complete))
  for (i in names(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Pairwise)){
    tab.hubs.comparison[1:length(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Pairwise[[i]]),i] <- 
      coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Pairwise[[i]]
  }
  for (i in names(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Complete)){
    tab.hubs.comparison[1:length(coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Complete[[i]]),i] <- 
      coppi.obj$Centralities$Hubs.percentile$Hubs_Diff$Complete[[i]]
  }
  writeDataTable(Centr.wb,"Hubs_percentile_comparison",as.data.frame(tab.hubs.comparison))
  
  
  
  
  addWorksheet(Centr.wb,"Hubs_av_degree")
  H <- coppi.obj$Centralities$Hubs.av_degree$Hubs
  All.hubs <- unique(unlist(sapply(H, function(x) x[,1])))
  tab.hubs <- matrix(NA,nrow = length(All.hubs),ncol = length(H),
                     dimnames = list(All.hubs, names(H)))
  
  for (i in 1:length(H)){
    tab.hubs[H[[i]][,1],i] <- as.numeric(H[[i]][,2])
  }
  writeDataTable(Centr.wb,"Hubs_av_degree",as.data.frame(tab.hubs),
                 rowNames = T)
  
  
  addWorksheet(Centr.wb,"Hubs_av_degree_comparison")
  n.row <- max(sapply(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Pairwise,length)) 
  tab.hubs.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Pairwise)+
                                  length(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Complete))
  colnames(tab.hubs.comparison) <- c(names(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Pairwise),
                                     names(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Complete))
  for (i in names(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Pairwise)){
    tab.hubs.comparison[1:length(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Pairwise[[i]]),i] <- 
      coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Pairwise[[i]]
  }
  for (i in names(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Complete)){
    tab.hubs.comparison[1:length(coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Complete[[i]]),i] <- 
      coppi.obj$Centralities$Hubs.av_degree$Hubs_Diff$Complete[[i]]
  }
  writeDataTable(Centr.wb,"Hubs_av_degree_comparison",as.data.frame(tab.hubs.comparison))
  
  
  
  addWorksheet(Centr.wb,"Betweenness")
  writeDataTable(Centr.wb,"Betweenness", as.data.frame(coppi.obj$Centralities$Bet.ppi),
                 rowNames = T)
  
  addWorksheet(Centr.wb,"Bottleneck_av_betw")
  writeDataTable(Centr.wb,"Bottleneck_av_betw",as.data.frame(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck),
                 rowNames = T)
  
  addWorksheet(Centr.wb,"Bottleneck_av_betw_comparison")
  n.row <- max(sapply(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise,length)) 
  tab.bott.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise)+
                                  length(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete))
  colnames(tab.bott.comparison) <- c(names(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise),
                                     names(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete))
  for (i in names(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise)){
    tab.bott.comparison[1:length(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise[[i]]),i] <- 
      coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise[[i]]
  }
  for (i in names(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete)){
    tab.bott.comparison[1:length(coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete[[i]]),i] <- 
      coppi.obj$Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete[[i]]
  }
  
  writeDataTable(Centr.wb,"Bottleneck_av_betw_comparison",as.data.frame(tab.bott.comparison))
  
  addWorksheet(Centr.wb,"Edge betweenness")
  Ed.bt <- coppi$Centralities$Edge.bet.ppi
  name_edges <- unique(unlist(lapply(Ed.bt, names)))
  tab.ed.bet <- matrix(NA,nrow = length(name_edges),
                       ncol = length(Ed.bt),
                       dimnames = list(name_edges,names(Ed.bt)))
  
  for (i in 1:length(Ed.bt)){
    tab.ed.bet[names(Ed.bt[[i]]),i] <- as.numeric(Ed.bt[[i]])
  }
  writeDataTable(Centr.wb,"Edge betweenness", as.data.frame(tab.ed.bet),
                 rowNames = T)
  
  addWorksheet(Centr.wb,"Edge Bottleneck_av_betw")
  writeDataTable(Centr.wb,"Edge Bottleneck_av_betw",as.data.frame(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck),
                 rowNames = T)
  
  
  addWorksheet(Centr.wb,"EdBottneck_avbetw_comparison")
  n.row <- max(sapply(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise,length)) 
  tab.bott.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise)+
                                  length(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete))
  colnames(tab.bott.comparison) <- c(names(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise),
                                     names(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete))
  for (i in names(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise)){
    tab.bott.comparison[1:length(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise[[i]]),i] <- 
      coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise[[i]]
  }
  for (i in names(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete)){
    tab.bott.comparison[1:length(coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete[[i]]),i] <- 
      coppi.obj$Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete[[i]]
  }
  
  writeDataTable(Centr.wb,"EdBottneck_avbetw_comparison",as.data.frame(tab.bott.comparison))
  saveWorkbook(Centr.wb, 
               file = str_c("Centralities_",name_analysis,".xlsx"),
               overwrite = T)
}
