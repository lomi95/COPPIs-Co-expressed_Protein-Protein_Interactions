source("loading_function.R")
load("interactome_annotations_hs.RDa")
setwd("Netosis/")

Norm_IceR <- read.xlsx("DADA_IceR_imputata_normalizzata_manipolata.xlsx")
Norm_IceR <- Norm_IceR[-which(Norm_IceR$gene.name==NaN),] ## non funziona con is.nan
Norm_IceR <- aggregate(Norm_IceR[], list(Norm_IceR$gene.name), FUN = max, na.rm = TRUE)
rownames(Norm_IceR) <- Norm_IceR[,1]
summary_interactome <- ppi.int

dataset <- Norm_IceR[,3:ncol(Norm_IceR)]
genes_id <- Norm_IceR$Group.1
corr.test <- "spearman"
rownames(dataset) <- genes_id

prot.annotation <- lapply(summary_interactome$annotation.interactome,function(x){
  compute_annotation(genes_id,x)
})

dataset.t <- t(dataset)
## correlazione




Corr <- rcorr(dataset.t,type = corr.test)
flattCorr <- flattenCorrMatrix(Corr$r,Corr$P,Corr$n)
flattCorr$cor_features <- paste(flattCorr$row,"and",flattCorr$column)

graph.Corr <- graph_from_edgelist(as.matrix(flattCorr[,1:2]),directed = F)

edge.attributes(graph.Corr)$p <- flattCorr$p

### graph PPI
graph.interactome <- graph_from_edgelist(as.matrix(ppi.int$interactome[,3:4]),directed = F)

### intersection graph correlation - PPI
Corr.PPI <- intersection(graph.Corr, graph.interactome, keep.all.vertices = F)
### difference corr- ppi
Corr.notPPI <- difference(graph.Corr, Corr.PPI)

N.corrPPI <- length(E(Corr.PPI))
N.corrnotPPI <- length(E(Corr.notPPI))
sign.PPI <- sum(edge.attributes(Corr.PPI)$p < 0.05)
sign.notPPI <- sum(edge.attributes(Corr.notPPI)$p < 0.05)

contingency.table <- cbind(c(sign.PPI,N.corrPPI),
                           c(sign.notPPI,N.corrnotPPI))


fish.test <- fisher.test(contingency.table)
fish.test




# vediamo se le correlazioni all'interno di component e compartment sono 
# in media maggiori di una distribuzione random
graph.Corr.1 <- graph_from_edgelist(as.matrix(flattCorr[,1:2]),directed = F)
edge.attributes(graph.Corr.1)$cor <- flattCorr$cor

embe <- lapply(prot.annotation, function(ann){
  
  Corr.Compartments <- t(sapply(ann$genes_found, function(x){
    ind.sign <- edge.attributes(induced.subgraph(graph.Corr,x))$p < 0.05
    if (sum(ind.sign)){
      return(cbind(mean(abs(edge.attributes(induced.subgraph(graph.Corr.1,x))$cor[ind.sign])),
                   sum(ind.sign),length(ind.sign)))
    } else {
      return(cbind(0,sum(ind.sign),length(ind.sign)))
    }
  }))
  rownames(Corr.Compartments) <- ann$description
  colnames(Corr.Compartments) <- c("mean terms","n_sign_corr","length terms")
  
  return(Corr.Compartments)
})

vectCorr.signratio <- c(sign.notPPI,N.corrnotPPI)
f.test.sign <- lapply(embe, function(x){
  f.i <- apply(x,1,function(y){
    cont.table <- cbind(y[2:3],vectCorr.signratio)
    return(fisher.test(cont.table))
  })
})

embe1 <- embe
for (i in names(embe)){
  embe1[[i]] <- as.data.frame(embe1[[i]])
  embe1[[i]]$pvalue.ft  <- sapply(f.test.sign[[i]], function(x) x$p.value) 
  embe1[[i]]$pvalue.adj <- p.adjust(embe1[[i]]$pvalue.ft)
}

embe.sign <- lapply(embe1, function(x){
  x[x$pvalue.adj < 0.05,]
})


n.start.terms <- sapply(embe, nrow)
n.sign.terms  <- sapply(embe.sign, nrow)

perc.sign <- round(n.sign.terms/n.start.terms*100, 2)
