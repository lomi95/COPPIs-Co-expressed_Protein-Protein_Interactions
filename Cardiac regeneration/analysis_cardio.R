source("loading_function.R")
name_analysis <- "Cardiac_regeneration"
threshold_edge_cyto <- 0.6
tax_id <- 10090
names_of_groups <- c("mock","mir199","mir590")
name_genes <- "Gene"
Cardio <- as.data.frame(read_excel("Files/DataMatrixCmyo_2.xlsx"))
# keeping the highest value between duplicated genes
Cardio <- aggregate(Cardio[], list(Cardio[,name_genes]), FUN = max, na.rm = TRUE)
rownames(Cardio) <- Cardio[,1]


setwd("Cardiac regeneration")
load(str_c(name_analysis,".coppi.RDa"))
coppi <- all_analysis(Cardio, tax_ID = tax_id,
                      name_genes = name_genes,
                      names_of_groups =  names_of_groups,
                      cor_groups = coppi$cor.all,
                      sig.all = coppi$SIG.all,
                      interactome.01 = coppi$Interactome)
save(coppi, file = str_c(name_analysis,".coppi.RDa"))
save_excel(coppi,name_analysis = name_analysis)
all_excels <- list.files()[grep(".xlsx",list.files())]
names(all_excels) <- all_excels

# open cytoscape first
cytoscapePing()
closeSession(save.before.closing = F)
mm <- lapply(all_excels,Load2Cytoscape, threshold_edge_cyto)

# mappare score a colore
setCytoStyle()
saveSession(str_c(name_analysis,"_",threshold_edge_cyto,"filt.cys"))
# clusterizzare secondo louvain e taggare cluster
tm_components()
saveSession(str_c(name_analysis,"_tm.cys"))


setwd("..")

