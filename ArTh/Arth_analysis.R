rm(list=ls())
source("loading_function.R")
load("interactome_annotations_arth.RDa")
setwd("ArTh/")
name_analysis <- "ArTh"

data <- read.xlsx("Arabidopsis.xlsx")
rownames(data) <- data$`Gene.names.(primary.)`

CoPPI_results <- pipeline(dataset = data,
                          genes_id = data$`Gene.names.(primary.)`,
                          names_of_groups =  c("wt_con","wt_lyn","gun_con","gun_lyn"),
                          summary_interactome = ppi.int_arth,
                          name_analysis = "Arth",
                          categories = c("Component","Function","Process","RCTM","WikiPathways"))
