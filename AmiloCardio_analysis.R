rm(list=ls())
source("loading_function.R")
load("interactome_annotations_hs.RDa")
setwd("Cardiac Amiloidosis/")
AmiloCardio <- read.xlsx("amilocardio_mean_replicate.xlsx",rowNames = T)

countzeros <- apply(AmiloCardio,1, function(x) sum(x==0,na.rm = T))
countNAs   <- apply(AmiloCardio,1, function(x) sum(is.na(x)))

missingValues <- countNAs+countzeros

AmiloCardio.full <- AmiloCardio[missingValues==0,]

categories = c("Component","Process","RCTM","WikiPathways", "CORUM")
CoPPI_results <- pipeline(dataset = AmiloCardio.full,
                          names_of_groups = c("Control","AL","TTR"),
                          genes_id = AmiloCardio.full$Gene.Name,
                          summary_interactome = ppi.int, 
                          categories = categories,
                          name_analysis = "AmiloCardio_nonCorretto", 
                          correctionCorr = NULL)
