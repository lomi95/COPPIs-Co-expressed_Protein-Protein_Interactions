save_excel.centralities <- function(Centralities, name_analysis){
  
  Centr.wb <- createWorkbook()
  
  
  addWorksheet(Centr.wb,"Hubs_percentile")
  H <- Centralities$Hubs.percentile$Hubs
  All.hubs <- unique(unlist(sapply(H, function(x) x[,1])))
  tab.hubs <- matrix(NA,nrow = length(All.hubs),ncol = length(H),
                     dimnames = list(All.hubs, names(H)))
  
  for (i in 1:length(H)){
    tab.hubs[H[[i]][,1],i] <- as.numeric(H[[i]][,2])
  }
  writeDataTable(Centr.wb,"Hubs_percentile",as.data.frame(tab.hubs),
                 rowNames = T)
  
  
  addWorksheet(Centr.wb,"Hubs_percentile_comparison")
  n.row <- max(sapply(Centralities$Hubs.percentile$Hubs_Diff$Pairwise,length)) 
  tab.hubs.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(Centralities$Hubs.percentile$Hubs_Diff$Pairwise)+
                                  length(Centralities$Hubs.percentile$Hubs_Diff$Complete))
  colnames(tab.hubs.comparison) <- c(names(Centralities$Hubs.percentile$Hubs_Diff$Pairwise),
                                     names(Centralities$Hubs.percentile$Hubs_Diff$Complete))
  for (i in names(Centralities$Hubs.percentile$Hubs_Diff$Pairwise)){
    if (length(Centralities$Hubs.percentile$Hubs_Diff$Pairwise[[i]])){
      tab.hubs.comparison[1:length(Centralities$Hubs.percentile$Hubs_Diff$Pairwise[[i]]),i] <- 
        Centralities$Hubs.percentile$Hubs_Diff$Pairwise[[i]]
    }
  }
  for (i in names(Centralities$Hubs.percentile$Hubs_Diff$Complete)){
    if (length(Centralities$Hubs.percentile$Hubs_Diff$Complete[[i]])){
      tab.hubs.comparison[1:length(Centralities$Hubs.percentile$Hubs_Diff$Complete[[i]]),i] <- 
        Centralities$Hubs.percentile$Hubs_Diff$Complete[[i]]
    }
  }
  writeDataTable(Centr.wb,"Hubs_percentile_comparison",as.data.frame(tab.hubs.comparison))
  
  
  
  
  addWorksheet(Centr.wb,"Hubs_av_degree")
  H <- Centralities$Hubs.av_degree$Hubs
  All.hubs <- unique(unlist(sapply(H, function(x) x[,1])))
  tab.hubs <- matrix(NA,nrow = length(All.hubs),ncol = length(H),
                     dimnames = list(All.hubs, names(H)))
  
  for (i in 1:length(H)){
    if (length(as.numeric(H[[i]][,2]))){
      tab.hubs[H[[i]][,1],i] <- as.numeric(H[[i]][,2])
    }
  }
  writeDataTable(Centr.wb,"Hubs_av_degree",as.data.frame(tab.hubs),
                 rowNames = T)
  
  
  addWorksheet(Centr.wb,"Hubs_av_degree_comparison")
  n.row <- max(sapply(Centralities$Hubs.av_degree$Hubs_Diff$Pairwise,length)) 
  tab.hubs.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(Centralities$Hubs.av_degree$Hubs_Diff$Pairwise)+
                                  length(Centralities$Hubs.av_degree$Hubs_Diff$Complete))
  colnames(tab.hubs.comparison) <- c(names(Centralities$Hubs.av_degree$Hubs_Diff$Pairwise),
                                     names(Centralities$Hubs.av_degree$Hubs_Diff$Complete))
  for (i in names(Centralities$Hubs.av_degree$Hubs_Diff$Pairwise)){
    if (length(Centralities$Hubs.av_degree$Hubs_Diff$Pairwise[[i]])){
      tab.hubs.comparison[1:length(Centralities$Hubs.av_degree$Hubs_Diff$Pairwise[[i]]),i] <- 
        Centralities$Hubs.av_degree$Hubs_Diff$Pairwise[[i]]
    }
  }
  for (i in names(Centralities$Hubs.av_degree$Hubs_Diff$Complete)){
    if (length(Centralities$Hubs.av_degree$Hubs_Diff$Complete[[i]])){
      tab.hubs.comparison[1:length(Centralities$Hubs.av_degree$Hubs_Diff$Complete[[i]]),i] <- 
        Centralities$Hubs.av_degree$Hubs_Diff$Complete[[i]]
    }
  }
  writeDataTable(Centr.wb,"Hubs_av_degree_comparison",as.data.frame(tab.hubs.comparison))
  
  
  
  addWorksheet(Centr.wb,"Betweenness")
  writeDataTable(Centr.wb,"Betweenness", as.data.frame(Centralities$Bet.ppi),
                 rowNames = T)
  
  addWorksheet(Centr.wb,"Bottleneck_av_betw")
  writeDataTable(Centr.wb,"Bottleneck_av_betw",as.data.frame(Centralities$Bottleneck.av_betw$Bottleneck),
                 rowNames = T)
  
  addWorksheet(Centr.wb,"Bottleneck_av_betw_comparison")
  n.row <- max(sapply(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise,length)) 
  tab.bott.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise)+
                                  length(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete))
  colnames(tab.bott.comparison) <- c(names(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise),
                                     names(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete))
  for (i in names(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise)){
    if (length(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise[[i]])){
      tab.bott.comparison[1:length(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise[[i]]),i] <- 
        Centralities$Bottleneck.av_betw$Bottleneck_Diff$Pairwise[[i]]
    }
  }
  for (i in names(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete)){
    if (length(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete[[i]])){
      tab.bott.comparison[1:length(Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete[[i]]),i] <- 
        Centralities$Bottleneck.av_betw$Bottleneck_Diff$Complete[[i]]
    }
  }
  
  writeDataTable(Centr.wb,"Bottleneck_av_betw_comparison",as.data.frame(tab.bott.comparison))
  
  addWorksheet(Centr.wb,"Edge betweenness")
  Ed.bt <- Centralities$Edge.bet.ppi
  name_edges <- unique(unlist(lapply(Ed.bt, names)))
  tab.ed.bet <- matrix(NA,nrow = length(name_edges),
                       ncol = length(Ed.bt),
                       dimnames = list(name_edges,names(Ed.bt)))
  
  for (i in 1:length(Ed.bt)){
    if (length(as.numeric(Ed.bt[[i]]))){
      tab.ed.bet[names(Ed.bt[[i]]),i] <- as.numeric(Ed.bt[[i]])
    }
  }
  writeDataTable(Centr.wb,"Edge betweenness", as.data.frame(tab.ed.bet),
                 rowNames = T)
  
  addWorksheet(Centr.wb,"Edge Bottleneck_av_betw")
  writeDataTable(Centr.wb,"Edge Bottleneck_av_betw",as.data.frame(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck),
                 rowNames = T)
  
  
  addWorksheet(Centr.wb,"EdBottneck_avbetw_comparison")
  n.row <- max(sapply(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise,length)) 
  tab.bott.comparison <- matrix(NA,nrow = n.row, 
                                ncol = length(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise)+
                                  length(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete))
  colnames(tab.bott.comparison) <- c(names(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise),
                                     names(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete))
  for (i in names(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise)){
    if (length(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise[[i]])){
      tab.bott.comparison[1:length(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise[[i]]),i] <- 
        Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Pairwise[[i]]
    }
  }
  for (i in names(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete)){
    if (length(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete[[i]])){
      tab.bott.comparison[1:length(Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete[[i]]),i] <- 
        Centralities$Edge.Bottleneck.av_edbetw$Edge.Bottleneck_Diff$Complete[[i]]
    }
  }
  
  writeDataTable(Centr.wb,"EdBottneck_avbetw_comparison",as.data.frame(tab.bott.comparison))
  saveWorkbook(Centr.wb, 
               file = str_c("Centralities_",name_analysis,".xlsx"),
               overwrite = T)
}
