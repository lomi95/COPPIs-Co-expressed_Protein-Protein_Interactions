source("loading_function.R")


Cardio <- as.data.frame(read_excel("Files/DataMatrixCmyo_2.xlsx"))
# keeping the highest value between duplicated genes
Cardio <- aggregate(Cardio[], list(Cardio$Gene), FUN = max, na.rm = TRUE)
rownames(Cardio) <- Cardio[,1]
Cardio_t <- t(Cardio)

all_categ.cardio <- all_analysis(Cardio_t,
                                 names_of_groups =  c("Mock","miR590","miR199"),
                                 tax_ID = 10090)



save_excel(all_categ.cardio)

all_excels <- list.files()[grep(".xlsx",list.files())]
names(all_excels) <- all_excels

cytoscapePing()
mm <- lapply(all_excels,Load2Cytoscape, 0.6)
