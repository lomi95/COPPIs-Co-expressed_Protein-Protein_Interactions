save_results <- function(category, name_analysis){
  name_categ <- category$table_summary$category[1]
  if (!is.null(name_categ)){
    category$table_summary$category <- NULL
    category$table_summary$inputGenes <- NULL
    summary.terms <- createWorkbook()
    # nell'enrichment deve esserci 
    # uno sheet con la tabella con tutti i pathways risultati significativi 
    name.sheet <- "node_table.terms"
    
    addWorksheet(summary.terms, name.sheet)
    writeDataTable(summary.terms, name.sheet, 
                   category$table_summary)
    
    for (i in 1:length(category$graph.similarity)){
      if (!is.null(dim(category$graph.similarity[[i]]$el_graph.term))){
        name.sheet <- names(category$graph.similarity)[[i]]
        addWorksheet(summary.terms, name.sheet)
        writeDataTable(summary.terms, name.sheet, 
                       category$graph.similarity[[i]]$el_graph.term)
      }
    }
    
    saveWorkbook(summary.terms, 
                 str_c(name_analysis,"_",
                       name_categ,".xlsx"), 
                 overwrite = TRUE) 
  }
}
