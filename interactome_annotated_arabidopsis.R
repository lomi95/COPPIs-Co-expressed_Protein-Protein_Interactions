interactome.arth <- read.delim("3702.protein.links.detailed.v12.0.txt",sep = " ")
score_interactome_type = c("experimental","database")
score_interactome_type_threshold = c(0.15,0.3)*1000
interactome.arth <- filter_interactome(unique(interactome.arth),
                                       score_interactome_type,
                                       score_interactome_type_threshold)


all.interacting.proteins <- unique(c(interactome.arth$protein1,
                                     interactome.arth$protein2))

vect.ann <- seq(1,length(all.interacting.proteins),by = 1999)

all.annotation <- list()
all.mapping <- list()
for (i in 1:(length(vect.ann)-1)){
  all.annotation[[i]] <- rba_string_annotations(
    all.interacting.proteins[vect.ann[i]:(vect.ann[i+1]-1)],3702,verbose = F, split_df = F)
  all.mapping[[i]] <- rba_string_map_ids(
    all.interacting.proteins[vect.ann[i]:(vect.ann[i+1]-1)],3702,verbose = F)
}
all.annotation[[i+1]] <- rba_string_annotations(
  all.interacting.proteins[vect.ann[i+1]:length(all.interacting.proteins)],3702, split_df = F)
all.mapping[[i+1]] <- rba_string_map_ids(
  all.interacting.proteins[vect.ann[i+1]:length(all.interacting.proteins)],3702)

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
map_id.merged <- Reduce(rbind,all.mapping)
annotation.df <- Reduce(merge_annotation,all.annotation)
annotation.df <- annotation.df[annotation.df$number_of_genes>1,]

# save(annotation.df, file = "annotation_df.RDa")

### TROPPO PESANTE ###
# annotation.df$graph.ppi <- lapply(annotation.df$inputGenes,function(x){
#   induced.subgraph(graph.interactome.arth,x)
# })

# unire compartments e components?


ppi.int_arth <- get_ppi_interactions(interactome.arth,annotation.df,map_id.merged)
#save(ppi.int_arth, file = "interactome_annotations_arth.RDa")
