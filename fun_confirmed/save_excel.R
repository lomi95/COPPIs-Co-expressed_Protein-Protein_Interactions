save_excel <- function(all_categories, name_analysis){
  # creare un file per ogni enrichment
  for (i in names(all_categories$SIG.all)){
    if (nrow(as.data.frame(all_categories$SIG.all[[i]]$tab_summary))){
      
      
      summary_pathways <- createWorkbook()
      # nell'enrichment deve esserci 
      # uno sheet con la tabella con tutti i pathways risultati significativi 
      name_sheet <- "node_table_pathways"
      addWorksheet(summary_pathways, name_sheet)
      writeDataTable(summary_pathways, name_sheet, as.data.frame(all_categories$SIG.all[[i]]$tab_summary))
      
      # uno sheet per ogni confronto gpp, con score di similarità, score cor
      for (j in 1:length(all_categories$gpp.all[[i]]$pairwise)){
        name_sheet <- names(all_categories$gpp.all[[i]]$pairwise)[j]
        addWorksheet(summary_pathways, name_sheet)
        writeDataTable(summary_pathways, name_sheet, 
                       as.data.frame(all_categories$gpp.all[[i]]$pairwise[[j]]$el_gpp))
        
      }
      if (length(all_categories$gpp.all[[i]])>1){
        for (j in 1:length(all_categories$gpp.all[[i]]$complete$all_el)){
          name_sheet <- names(all_categories$gpp.all[[i]]$complete$all_el)[j]
          addWorksheet(summary_pathways, name_sheet)
          writeDataTable(summary_pathways, name_sheet, 
                         as.data.frame(all_categories$gpp.all[[i]]$complete$all_el[[j]]))
          
        }
      }
      
      
      saveWorkbook(summary_pathways, str_c(name_analysis,"_",i,".xlsx"), overwrite = TRUE)
    }
  }
}
