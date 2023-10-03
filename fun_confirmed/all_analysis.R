all_analysis <- function(dataset, # input: matrix (colnames = groups, rows = proteins)
                         name_genes,
                         tax_ID,
                         names_of_groups = NULL,
                         pos.vectors_groups = NULL, #list
                         ignore_case = T,
                         just_significant_cor = F,
                         interactome.01 = NULL,
                         score_interactome_threshold = 0,
                         score_interactome_type = c("escore","dscore"),
                         score_interactome_type_threshold = c(0,0), # vector score for type
                         cor_groups = NULL,
                         compute_weights = T,
                         cor.test = "spearman",
                         significance = 0.05,
                         p.adj = "BH",
                         sig.all = NULL,
                         split_categories = T,
                         categories = c("KEGG","Process","RCTM","WikiPathways"),
                         min_edge = 1, # numero di interazioni minime per essere considerato un pathway
                         max_edge = Inf,
                         filter_percentage = 0
){
  # Parameter Control:
  #   The code starts by performing various checks on the input parameters 
  #   to ensure they meet certain criteria. For instance, it checks whether 
  #   the provided correlation test type is either "spearman" or "pearson," 
  #   whether the p.adj method is valid, whether the specified categories 
  #   are valid, and more.
  
  if (length(as.vector(name_genes))==1){
    if (is.character(name_genes)){
      name_genes <- match.arg(name_genes,colnames(dataset))
      genes.1 <- dataset[,name_genes]
      dataset <- dataset[,-match(name_genes,colnames(dataset))]
      message("The genes column starts with: ", str_c(head(genes.1),collapse = " "),"\n")
      
    } else if (is.numeric(name_genes)){
      genes.1 <- dataset[,name_genes]
      dataset <- dataset[,-name_genes]
      message("The genes column starts with: ", str_c(head(genes.1),collapse = " "),"\n")
    } else {
      stop("Something's wrong with Gene name class")
    }
  } else if (length(as.vector(name_genes)) == nrow(dataset)){
    if (length(as.vector(name_genes)) != length(unique(as.vector(name_genes)))){
      message("Warning : There may be duplicated genes")
    }
    genes.1 <- name_genes
    message("The genes column starts with: ",str_c(head(genes.1),collapse = " "),"\n")
  } else if (length(as.vector(name_genes)) != nrow(dataset)){
    stop("The length of genes differs from the rows of the dataset\n")
  }
  
  map_genes <- rba_string_map_ids(genes.1,tax_ID)
  if (length(setdiff(genes.1,map_genes$queryItem))){
    message(str_c(setdiff(genes.1,map_genes$queryItem), collapse = " "), " were not found and dropped\n") 
  }

  
  
  dataset <- t(dataset[match(map_genes$queryItem,genes.1),])
  
  genes.2 <- genes.1[match(map_genes$queryItem,genes.1)]
  genes <- map_genes$preferredName[match(genes.2,map_genes$queryItem)]
  colnames(dataset) <- genes
  cor.test_available <- c("spearman","pearson")
  cor.test <- match.arg(cor.test, cor.test_available)
  
  p.adj <- match.arg(p.adj,p.adjust.methods)
  
  
  categories_available <- c("all","COMPARTMENTS","Component","DISEASES","Function",
                            "HPO","InterPro","KEGG","Keyword","NetworkNeighborAL",
                            "Pfam","PMID","Process","RCTM",
                            "SMART","TISSUES","WikiPathways")
  categories <- match.arg(categories,categories_available, several.ok = T)
  
  if (sum("all"==categories)){
    categories <- NULL
  }
  
  if (is.null(pos.vectors_groups)){
    if (is.null(names_of_groups)){
      stop("\nNeither names_of_groups nor pos.vector_groups are specified")
    }
    names(names_of_groups) <- names_of_groups
    list_groups <- lapply(names_of_groups,function(x){
      dataset[grep(x,rownames(dataset), ignore.case = ignore_case),]
    })
    if (sum(sapply(list_groups,length)==0)){
      not_found <- names_of_groups[which(sapply(list_groups,length)==0)]
      stop(str_c("\n",str_c(not_found,collapse = " ")," were not found in the dataset rownames"))
    }
    
  } else {
    if (!is.list(pos.vectors_groups)){
      stop("pos.vector_groups must be a list with the row indexes of the groups")
    }
    if (!is.null(names_of_groups)){
      message("\nWarning! Both names_of_groups and pos.vectors_groups were specified, pos.vectors_groups will be used")
    }
    if (is.null(names(pos.vectors_groups))){
      names_of_groups <- str_c("Group_",seq(1,length(pos.vectors_groups)))
      names(names_of_groups) <- names_of_groups
    } else {
      names_of_groups <- names(pos.vectors_groups)
      names(names_of_groups) <- names_of_groups
    }
    list_groups <- lapply(pos.vectors_groups,function(x){
      dataset[x,]
    })
    names(list_groups) <- names_of_groups
  }
  
  score_interactome_avaliable <- c("score","nscore","fscore","pscore",
                                   "ascore","escore","dscore","tscore")
  
  score_interactome_type <- match.arg(score_interactome_type,
                                      score_interactome_avaliable,
                                      several.ok = T)
  
  # Data Processing:
  #   
  #   The code processes the input dataset and separates it into groups based 
  #   on the specified group names or row indices (pos.vectors_groups).
  #   It then filters and extracts interactions from the protein interaction 
  #   network (interactome) based on certain interaction scores.
  #   Enrichment analysis is performed on the proteins using the STRING database, 
  #   and the enrichment results are filtered based on specified categories.
  
  if (is.null(interactome.01)){
    interactome.01 <- rba_string_interactions_network(genes,
                                                      species = tax_ID,
                                                      required_score = score_interactome_threshold)
  }
  
  # filtrare l'interattoma prendendo solo interazioni con score di esperimento e database
  interactome <- filter_score_or(unique(interactome.01),
                                 score_interactome_type,
                                 score_interactome_type_threshold)
  
  
  # enrichment su proteine
  
  if (split_categories){
    prot_enr <- rba_string_enrichment(genes,species = tax_ID)
    if (!is.null(categories)){
      prot_enr <- prot_enr[categories]
    }
    if (sum(is.na(names(prot_enr)))){
      message("\nWarning : there are no process enriched in\n ", 
              str_c(categories[which(is.na(names(prot_enr)))],collapse = "\n"))
      prot_enr <- prot_enr[-which(is.na(names(prot_enr)))]
    }
  } else {
    prot_enr <- list(All_categories = rba_string_enrichment(genes,
                                                            species = tax_ID,
                                                            split_df = split_categories))
  }
  
  # Correlation Calculation:
  #   
  #   If cor_groups is not provided, the code calculates protein correlation 
  #   within each group using either Spearman or Pearson correlation test. 
  #   It also applies p-value p.adj.
  #   The calculated correlations are stored in cor_groups.
  
  if (is.null(cor_groups)){
    start_time <- Sys.time() 
    message("\nComputing protein correlation for each group")
    start_time <- Sys.time()
    cor_groups <- cor_fun(list_groups,
                          test_type = cor.test,
                          significance = significance,
                          p.adj = p.adj,
                          compute_weights = compute_weights)
    time_employed <- format_time_difference(difftime(Sys.time(),start_time,units = "s"))
    message(str_c("Correlation computed in", time_employed,sep = " "))
  } 
  # Pathway Extraction and Filtering:
  #   The code extracts pathways enriched with proteins showing specific
  #   correlations (sig.all). The pathways are filtered based on the number 
  #   of edges (filter_edge) and percentage (filter_percentage).
  sd_groups <- sapply(cor_groups$cor_cond,function(x) sd(abs(x$ALL$weights)))
  
  if (is.null(sig.all)){
    start_time <- Sys.time() 
    message("\nExtraction of upcorrelated pathways")
    sig.all <- lapply(prot_enr,pipe_for_term, 
                      cor_groups1 = cor_groups$cor_cond,
                      sd_groups = sd_groups,
                      graph_cor.groups = cor_groups$cor_graphs, 
                      interac = interactome,
                      p.adj = p.adj,
                      significance = significance,
                      min_edge = min_edge, 
                      max_edge = max_edge,
                      filter_percentage = filter_percentage,
                      js = just_significant_cor)
    time_employed <- format_time_difference(difftime(Sys.time(),start_time,units = "s"))
    message(str_c("Extraction completed in", time_employed,sep = " "))
  }
  
  # Graph Construction:
  #   
  #   The code constructs pairwise graphs between pathways and computes a graph 
  #   that compares all pathways. The graphs represent relationships between 
  #   pathways based on protein interactions.
  start_time <- Sys.time() 
  message("\nCreating pathways-pathways graphs")
  gpp.all <- lapply(sig.all, function(x){
    if (!is.null(x$sig_path.p)){
      gpp.pair <- lapply(x$sig_path.p, function(y){
        if (length(y[,2])>0){
          gpp.pairwise <- gpp_fun(as.data.frame(x$tab_pathprot)[y[,1],],y)
        } else {
          gpp.pairwise <- list(gpp = make_empty_graph(directed = F),
                               el_gpp = matrix(nrow=0,ncol=7))
        }
      })
      if (length(names_of_groups)>2){
        gpp.complete <- gpp.compare(gpp.pair,names_of_groups,ignore_case)
        return(list(pairwise = gpp.pair,
                    complete = gpp.complete))
      } else {
        return(list(pairwise = gpp.pair))
      }
    }
  })
  time_employed <- format_time_difference(difftime(Sys.time(),start_time,units = "s"))
  message(str_c("Pathways-pathways graphs created in", time_employed,sep = " "))
  
  
  return(list(SIG.all = sig.all,
              cor.all = cor_groups,
              gpp.all = gpp.all,
              Enrichment  = prot_enr,
              Interactome = interactome))
}
