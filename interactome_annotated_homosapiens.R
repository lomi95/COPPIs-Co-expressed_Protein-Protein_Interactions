source("loading_function.R")
interactome.hs <- read.delim("9606.protein.links.detailed.v12.0.txt",sep = " ")
score_interactome_type = c("experimental","database")
score_interactome_type_threshold = c(0.15,0.3)*1000
interactome.hs <- filter_interactome(unique(interactome.hs),
                                     score_interactome_type,
                                     score_interactome_type_threshold)


all.interacting.proteins <- unique(c(interactome.hs$protein1,
                                     interactome.hs$protein2))

vect.ann <- seq(1,length(all.interacting.proteins),by = 1999)

all.annotation <- list()
all.mapping <- list()
for (i in 1:(length(vect.ann)-1)){
  all.annotation[[i]] <- rba_string_annotations(
    all.interacting.proteins[vect.ann[i]:(vect.ann[i+1]-1)],9606,verbose = F, split_df = F)
  all.mapping[[i]] <- rba_string_map_ids(
    all.interacting.proteins[vect.ann[i]:(vect.ann[i+1]-1)],9606,verbose = F)
}
all.annotation[[i+1]] <- rba_string_annotations(
  all.interacting.proteins[vect.ann[i+1]:length(all.interacting.proteins)],9606, split_df = F)
all.mapping[[i+1]] <- rba_string_map_ids(
  all.interacting.proteins[vect.ann[i+1]:length(all.interacting.proteins)],9606)

merge_annotation <- function(df1,df2){
  df.m <- merge(df1,df2, by = "term",all = T)
  
  col.original <- c("category","description")
  col.genes    <- c("inputGenes","preferredNames")
  
  new_col <- list()
  for (i in col.original){
    col.i <- df.m[,grep(i,colnames(df.m))]
    new_col[[i]] <- apply(col.i,1,function(x){
      x[which(!is.na(x))[1]]
    })
  }
  union.not.na <- function(x,y){
    as.vector(na.omit(union(x,y)))
  }
  for (i in col.genes){
    col.i <- df.m[,grep(i,colnames(df.m))]
    new_col[[i]] <- mapply(union.not.na, col.i[,1], col.i[,2])
  }
  new_col.df <- data.frame(term = df.m$term,
                           category = new_col$category,
                           description = new_col$description)
  new_col.df$inputGenes <- new_col$inputGenes
  new_col.df$preferredNames <- new_col$preferredNames
  new_col.df$number_of_genes <- sapply(new_col$preferredNames,length)
  return(new_col.df)
}
map_id.merged1 <- Reduce(rbind,all.mapping)
annotation.df <- Reduce(merge_annotation,all.annotation)
annotation.df1 <- annotation.df[annotation.df$number_of_genes>1,]

# save(annotation.df, file = "annotation_df.RDa")

### TROPPO PESANTE ###
# annotation.df$graph.ppi <- lapply(annotation.df$inputGenes,function(x){
#   induced.subgraph(graph.interactome.hs,x)
# })

# unire compartments e components?



Corum.downloaded <- read.xlsx("CORUM download 2022_09_12.xlsx")
Corum.organisms <- split(Corum.downloaded,Corum.downloaded$Organism)
Corum.human <- Corum.organisms$Human
gene_name.complexes <- str_split(Corum.human$`subunits(Gene.name)`, ";")

gene_name <- unique(unlist(str_split(Corum.human$`subunits(Gene.name)`, ";")))

mapGene.name1 <- rba_string_map_ids(gene_name, 9606)
mapGene.name  <- mapGene.name1[na.omit(match(all.interacting.proteins,mapGene.name1$stringId)),]
Corum.light <- data.frame(term = Corum.human$ComplexID,
                          category = rep("CORUM", length(Corum.human$ComplexID)),
                          description = Corum.human$ComplexName)
inputGenes <- lapply(gene_name.complexes, function(x){
  mapGene.name$stringId[na.omit(match(x,mapGene.name$queryItem))]
})


preferredName <- lapply(gene_name.complexes, function(x){
  mapGene.name$preferredName[na.omit(match(x,mapGene.name$queryItem))]
})
Corum.light$inputGenes      <- inputGenes
Corum.light$preferredNames  <- preferredName
Corum.light$number_of_genes <- sapply(preferredName, length)

Corum.light1 <- Corum.light[Corum.light$number_of_genes>1,]


annotation.df2 <- rbind.data.frame(annotation.df1,
                                   Corum.light1)

map_id.merged <- rbind(map_id.merged1, mapGene.name)



ppi.int <- get_ppi_interactions(interactome.hs,annotation.df2,map_id.merged)
#save(ppi.int, file = "interactome_annotations_hs.RDa")
