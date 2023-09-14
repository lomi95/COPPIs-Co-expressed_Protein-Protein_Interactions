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

SUID_list <- getNetworkList(getSUIDs = T)

#lapply(SUID_list,tm_components)

all_Suid.collection <- unique(sapply(SUID_list,getCollectionSuid))
default_param <- sapply(getVisualPropertyNames(),getVisualPropertyDefault,"default")
for (i in 1:length(all_Suid.collection)){
  print(all_Suid.collection[i])
  style.name  <- getCollectionName(all_Suid.collection[i])
  n_i <- getCollectionNetworks(all_Suid.collection[i])[1]
  nodeLabels <- mapVisualProperty('node label','id','p')
  ind_coll <- grep(style.name,names(mm))
  nodeFills <- mapVisualProperty(visual.prop = 'node fill color',table.column = 'scores',mapping.type = 'c',
                                 table.column.values = as.numeric(c(mm[[ind_coll]][1],mean(mm[[ind_coll]]),mm[[ind_coll]][2])),
                                 visual.prop.values = c('#FFF202','#FF7D00','#FC0000'))
  createVisualStyle(style.name, default_param, list(nodeLabels,nodeFills))
  lockNodeDimensions(FALSE, style.name)
  # lapply(getCollectionNetworks(all_Suid.collection[i]),function(x){
  #   setVisualStyle(style.name = style.name,network = x)})
}
