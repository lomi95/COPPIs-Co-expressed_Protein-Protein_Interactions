source("loading_function.R")


Cardio <- as.data.frame(read_excel("Files/DataMatrixCmyo_2.xlsx"))
# keeping the highest value between duplicated genes
Cardio <- aggregate(Cardio[], list(Cardio$Gene), FUN = max, na.rm = TRUE)
rownames(Cardio) <- Cardio[,1]
Cardio_t <- t(Cardio)

load("all_categ.cardio.RDa")
all_categ.cardio <- all_analysis(Cardio_t,
                                 names_of_groups =  c("Mock","miR590","miR199"),
                                 tax_ID = 10090, 
                                 cor_groups = all_categ.cardio$cor.all,
                                 node.comp_thr = 50)
#save(all_categ.cardio,file = "all_categ.cardio.RDa")



save_excel(all_categ.cardio)

all_excels <- list.files()[grep(".xlsx",list.files())][c(5,11,12,15)]
names(all_excels) <- all_excels

cytoscapePing()
mm <- lapply(all_excels,Load2Cytoscape, 0.6)
