all_analysis <- function(dataset, # input: matrix (rownames = groups, colnames = proteins)
                         names_of_groups = NULL,
                         pos.vectors_groups = NULL, #list
                         ignore_case = T,
                         interactome.01 = NULL,
                         tax_ID,
                         score_interactome_threshold = 0,
                         score_interactome_type = c("escore","dscore"), # vector type of score - lista da far vedere
                         score_interactome_type_threshold = c(0,0), # vector score for type
                         cor_groups = NULL,
                         compute_weights = T,
                         test_type = "spearman",
                         significance = 0.05,
                         correction = "BH",
                         sig.all = NULL,
                         categories = NULL,
                         min_edge = -Inf, # numero di interazioni minime per essere considerato un pathway
                         max_edge = Inf,
                         funzione_peso_gpp = gpp_fun,
                         filter_percentage = 0,
                         filter_score = F,
                         filter_gpp = F,
                         node.comp_thr = 50, 
                         quantile.node_thr = 0.15,
                         quantile.edge_thr = 0.15
){
  
  # Parameter Control:
  #   The code starts by performing various checks on the input parameters 
  #   to ensure they meet certain criteria. For instance, it checks whether 
  #   the provided correlation test type is either "spearman" or "pearson," 
  #   whether the correction method is valid, whether the specified categories 
  #   are valid, and more.
  
  if (!identical(funzione_peso_gpp, gpp_fun) && !identical(funzione_peso_gpp, fun_Wcor)) {
    stop("The weight function for pathway pathway graphs must be gpp_fun or fun_Wcor, 
         for more information read the documentation.\n
         The choice of this function doesn't influence the significant pathways")
  }
  
  
  categ_avaliable <- c("COMPARTMENTS","Component","DISEASES","Function",
                       "HPO","InterPro","KEGG","Keyword","NetworkNeighborAL",
                       "Pfam","PMID","Process","RCTM",
                       "SMART","TISSUES","WikiPathways") 
  
  if (sum(test_type==c("spearman","pearson"))==0){
    stop("The correlation test avaliable are spearman and pearson")
  }
  
  if (sum(correction==p.adjust.methods) != 1){
    stop(str_c("\nYou can use the following p.adjust.methods\n",paste(p.adjust.methods,collapse = ", ")))
  }
  
  if (!is.null(categories)){
    if (length(intersect(categories,categ_avaliable))!=length(categories)){
      not_found <- setdiff(categories,intersect(categories,categ_avaliable))
      stop(str_c("\n",paste(not_found, collapse = " ")," were not found", "\nThe categories avaliable are:\n", 
                 paste(categ_avaliable,collapse ="\n")))
    }
    if (length(categories) > length(categ_avaliable)){
      not_found <- setdiff(categories,categ_avaliable)
      stop(str_c("\n", paste(not_found,collapse = " ")," is not avaliable\n",
                 "The categories avaliable are:\n ", 
                 paste(categ_avaliable,collapse ="\n")))
    }
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
      stop(str_c("\n",not_found," were not found in the dataset rownames"))
    }
    
  } else {
    if (!is.list(pos.vectors_groups)){
      stop("pos.vector_groups must be a list with the row indexes of the groups")
    }
    if (!is.null(names_of_groups)){
      print("Warning! Both names_of_groups and pos.vectors_groups were specified, pos.vectors_groups will be used")
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
  if (!is.null(score_interactome_type)){
    if (length(intersect(score_interactome_type,score_interactome_avaliable))!=length(score_interactome_type)){
      not_found <- setdiff(score_interactome_type,intersect(score_interactome_type,score_interactome_avaliable))
      stop(str_c("\n",paste(not_found, collapse = " ")," were not found", "\nThe score_interactome_type avaliable are:\n", 
                 paste(score_interactome_avaliable,collapse ="\n")))
    }
    if (length(score_interactome_type) > length(score_interactome_avaliable)){
      not_found <- setdiff(score_interactome_type,score_interactome_avaliable)
      stop(str_c("\n", paste(not_found,collapse = " ")," is not avaliable\n",
                 "The score_interactome_type avaliable are:\n ", 
                 paste(score_interactome_avaliable,collapse ="\n")))
    }
  }
  
  
  # Data Processing:
  #   
  #   The code processes the input dataset and separates it into groups based 
  #   on the specified group names or row indices (pos.vectors_groups).
  #   It then filters and extracts interactions from the protein interaction 
  #   network (interactome) based on certain interaction scores.
  #   Enrichment analysis is performed on the proteins using the STRING database, 
  #   and the enrichment results are filtered based on specified categories.
  
  if (is.null(interactome.01)){
    interactome.01 <- rba_string_interactions_network(colnames(dataset),
                                                      species = tax_ID,
                                                      required_score = score_interactome_threshold)
  }
  
  
  # filtrare l'interattoma prendendo solo interazioni con score di esperimento e database
  interactome <- filter_score_or(unique(interactome.01),
                                 score_interactome_type,
                                 score_interactome_type_threshold)
  
  
  # enrichment su proteine (separato per categorie)
  prot_enr <- rba_string_enrichment(colnames(dataset),species = tax_ID,split_df = T)
  
  if (!is.null(categories)){
    prot_enr <- prot_enr[categories]
  }
  
  
  # Correlation Calculation:
  #   
  #   If cor_groups is not provided, the code calculates protein correlation 
  #   within each group using either Spearman or Pearson correlation test. 
  #   It also applies p-value correction.
  #   The calculated correlations are stored in cor_groups.
    
  if (is.null(cor_groups)){
    print("Computing protein correlation")
    cor_groups <- cor_fun(list_groups,
                          test_type = test_type,
                          significance = significance,
                          p.adj = correction,
                          compute_weights = compute_weights)
  } 
  
  # Pathway Extraction and Filtering:
  #   The code extracts pathways enriched with proteins showing specific
  #   correlations (sig.all). The pathways are filtered based on the number 
  #   of edges (filter_edge) and percentage (filter_percentage).
  sd_groups <- sapply(cor_groups$cor_cond,function(x) sd(abs(x$ALL$weights)))
  if (is.null(sig.all)){
    print("Extraction of upcorrelated pathways")
    sig.all <- lapply(prot_enr,pipe_for_term, 
                      cor_groups1 = cor_groups$cor_cond,
                      sd_groups = sd_groups,
                      graph_cor.groups = cor_groups$cor_graphs, 
                      interac = interactome,
                      p.adj = correction,
                      significance = significance,
                      min_edge = min_edge, 
                      max_edge = max_edge,
                      filter_percentage = filter_percentage)
  }
  for (i in names(sig.all)){
    attributes(sig.all[[i]])$Tag <- i
  }
  
  
  if (filter_score){
    # filtering pathways by score
    sig.all <- lapply(sig.all,function(x){
      x$sig_path.p <- lapply(x$sig_path.p,function(y){
        s.num <- as.numeric(y[,3])
        y <- y[s.num > mean(s.num),]
        return(y)
      })
      return(x)
    })
  }
  
  # Graph Construction:
  #   
  #   The code constructs pairwise graphs between pathways and computes a graph 
  #   that compares all pathways. The graphs represent relationships between 
  #   pathways based on protein interactions.
  
  print("Creating pathways-pathways graph")
  gpp.all <- lapply(sig.all, function(x){
    gpp.pair <- lapply(x$sig_path.p, function(y){
      if (length(y[,2])>0){
        gpp.pairwise <- funzione_peso_gpp(as.data.frame(x$tab_pathprot)[y[,1],],y, 
                                          filter_gpp,
                                          node.comp_thr = node.comp_thr, 
                                          quantile.node_thr = quantile.node_thr,
                                          quantile.edge_thr = quantile.edge_thr)
      } else {
        gpp.pairwise <- list(gpp = make_empty_graph(directed = F),
                             el_gpp = matrix(nrow=0,ncol=7))
      }
    })
    gpp.complete <- gpp.compare(gpp.pair,names_of_groups,ignore_case)
    
    return(list(pairwise = gpp.pair,
                complete = gpp.complete))
  })
  
  
  return(list(SIG.all = sig.all,
              cor.all = cor_groups,
              gpp.all = gpp.all))
}
