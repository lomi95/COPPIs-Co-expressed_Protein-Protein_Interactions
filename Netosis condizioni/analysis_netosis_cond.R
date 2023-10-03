source("loading_function.R")
name_analysis <- "Netosis_cond"
threshold_edge_cyto <- 0.6
tax_id <- 9606
names_of_groups =  c("dada","dctr","hd")

Norm_IceR <- as.data.frame(read_excel("Files/DADA_IceR_imputata_normalizzata_manipolata.xlsx"))
# removing gene name = Nan
Norm_IceR <- Norm_IceR[-which(Norm_IceR$`gene name`==NaN),] ## non funziona con is.nan
# keeping the highest value between duplicated genes
Norm_IceR <- aggregate(Norm_IceR[], list(Norm_IceR$`gene name`), FUN = max, na.rm = TRUE)

rownames(Norm_IceR) <- Norm_IceR[,1]
setwd("Netosis condizioni/")
#load(str_c(name_analysis,".coppi.RDa"))
coppi <- all_analysis(Norm_IceR, 
                      name_genes = "gene name",
                      names_of_groups =  names_of_groups,
                      tax_ID = tax_id)
save(coppi, file = str_c(name_analysis,".coppi.RDa"))
save_excel(coppi,name_analysis = name_analysis)
all_excels <- list.files()[grep(".xlsx",list.files())]
names(all_excels) <- all_excels

# open cytoscape first
cytoscapePing()
openSession()
mm <- lapply(all_excels,Load2Cytoscape, threshold_edge_cyto)

# mappare score a colore
setCytoStyle()
saveSession(str_c(name_analysis,"_",threshold_edge_cyto,"filt.cys"))
# clusterizzare secondo louvain e taggare cluster
tm_components()
saveSession(str_c(name_analysis,"_tm.cys"))


setwd("..")