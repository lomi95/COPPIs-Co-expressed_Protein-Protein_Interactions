pipeline <- function(dataset,
                     genes_id,
                     names_of_groups,
                     name_analysis,
                     summary_interactome,
                     gCOR.groups = NULL,
                     gCOR.categories = NULL,
                     signifCorr = 0.05, 
                     correctionCorr = "BH",
                     min_n.corr = NULL,
                     significance.coppi = 0.05,
                     perc.ned.nppi = 0.01,
                     perc.sign = 0.1,
                     correction.coppi = "BH",
                     corr.test = "spearman",
                     compute_weights = T,
                     min_edges = 1,
                     max_edges = Inf,
                     categories = c("Component","Process","RCTM","WikiPathways"),
                     threshold_edge_cyto = 0.6,
                     Cytoscape = T,
                     thr_degree = 0.75,
                     thr_corr = NULL){
  if (!is.null(min_n.corr) && compute_weights){
    message("Warning: there is a threshold for correlation pairs \n\t the compute_weights setting will be set to FALSE")
  }
  
  if (Cytoscape) {
    tryCatch(cytoscapePing(), error = function(err) {
      message("Open Cytoscape, or restart it")
    })
  }
  rownames(dataset) <- genes_id
  if (max_edges<min_edges){
    message("max edges is greater than min edges, we are switching them")
    m1 <- max_edges
    max_edges <- min_edges
    min_edges <- m1
  }
  
  ginter <- graph_from_edgelist(as.matrix(ppi.int$interactome[,3:4]),directed = F)
  start.time.annotation <- Sys.time()
  prot.annotation <- lapply(summary_interactome$annotation.interactome,function(x){
    compute_annotation(genes_id,x)
  })
  time_employed <- format_time_difference(difftime(Sys.time(),
                                                   start.time.annotation,units = "s"))
  message("Annotation computed in ", time_employed,sep = " ")
  
  names(prot.annotation) <- str_to_title(names(prot.annotation))
  categories <- str_to_title(categories)
  categories.avaliable <- names(prot.annotation)
  categ.notfound <- setdiff(categories,categories.avaliable)
  if (length(categ.notfound)){
    if (length(categ.notfound) == length(categories)){
      message("Error:\n",str_c(setdiff(categories,
                                       categories.avaliable), ", "), 
              "were not found in the avaliable categories")
      message("The categories avaliable are:\n", str_c(categories.avaliable, ", "))
      return(Categories)
    } else {
      message("Warning:\n",str_c(setdiff(categories,
                                         categories.avaliable), ", "), 
              "were not found in the avaliable categories")
      message("The categories avaliable are:\n", 
              str_c(setdiff(categories.avaliable,
                            categories), ", "))
      categories <- intersect(categories, categories.avaliable)
    }
  }
  
  dataset.t <- t(dataset)
  ## correlazione
  names(names_of_groups) <- names_of_groups
  list.groups <- lapply(names_of_groups, function(x){
    d <- dataset.t[grep(x,rownames(dataset.t), ignore.case = T),]
    if (nrow(d)){
      return(d)
    } else {
      message("No columns named ",x," were found in the dataset")
    }
  })
  list.groups <- list.groups[!sapply(list.groups, is.null)]
  
  if (length(list.groups) < 2){
    message("Error: The number of groups are less than 2, 
            please check the colnames(dataset)
            or the input names_of_groups")
    return(Groups)
  }
  start.time.correlation <- Sys.time()
  if (is.null(gCOR.categories)){
    CorrNULL <- matrix(1, nrow = ncol(list.groups[[1]]), ncol = ncol(list.groups[[1]]),
                       dimnames = list(colnames(list.groups[[1]]),colnames(list.groups[[1]])))
    flattNULL <- flattenCorrMatrix(CorrNULL,CorrNULL,CorrNULL)
    flattNULL$cor_features <- paste(flattNULL$row,"and",flattNULL$column)
    
    # intersection with interactome of categories
    gNULL.categories   <- lapply(summary_interactome$categories.interactome, function(y){
      
      gNULL <- graph_from_edgelist(as.matrix(flattNULL[,1:2]),directed = F)
      edge.attributes(gNULL)$features <- flattNULL$cor_features
      
      gNULL.inters <- intersection(gNULL, y, keep.all.vertices = F)
      ind.ppi <- match(sort(edge.attributes(gNULL.inters)$features),flattNULL$cor_features)
      gNULL.el <- flattNULL[ind.ppi,]
      
      return(list(el = gNULL.el,
                  length.interactome = length(E(y))))
    })
    names(gNULL.categories) <- str_to_title(names(gNULL.categories))
    if (is.null(gCOR.groups)){
      gCOR.groups.l <- lapply(list.groups,function(x){
        Corr <- rcorr(x,type = corr.test)
        flattCorr <- flattenCorrMatrix(Corr$r,Corr$P,Corr$n)
        flattCorr$cor_features <- paste(flattCorr$row,"and",flattCorr$column)
        flattCorr$cor[is.na(flattCorr$cor)] <- 0
        flattCorr$p[is.na(flattCorr$p)] <- 1
        
        if (!is.null(min_n.corr)){
          flattCorr$p[flattCorr$n < min_n.corr,] <- 1
          compute_weights <- T
        }
        
        if (is.null(correctionCorr)){
          flattCorr$p.adj <- flattCorr$p
        } else {
          flattCorr$p.adj <- p.adjust(flattCorr$p, correctionCorr)
        }
        ind_sig <- flattCorr$p.adj <= signifCorr
        
        if (compute_weights){
          ### normalizzazione per togliere gli 1 dal pvalue, per rendere possibile il logaritmo
          # L'obiettivo è trasformare l'intervallo del pvalue da (Significativià ; 1 ] a (0 ; 1)
          
          # Il nuovo massimo del pvalue sarà la media tra 1 e il secondo massimo
          if (max(flattCorr$p.adj)==1){
            new_max <- mean(c(1,max(flattCorr$p.adj[-which(flattCorr$p.adj==1)])))
          } else {
            new_max <- max(flattCorr$p.adj)
          }
          
          # calcolo parametri normalizzazione: (pvalue + epsilon ) * k = p_norm1
          # ( signifCorr + eps ) * k = signifCorr
          # ( max(pvalue)  + eps ) * k = new_max
          
          # risolvendo questo sistema si ottiene
          eps <- ( new_max - max(flattCorr$p.adj) ) * signifCorr / ( signifCorr - new_max )
          k   <- ( signifCorr - new_max ) / ( signifCorr - max(flattCorr$p.adj) )
          
          # applichiamo la normalizzazione solo ai pvalue non significativi
          ind_notsig <- which(flattCorr$p.adj > signifCorr)
          norm1 <- (flattCorr$p.adj[ind_notsig] + eps) * k
          
          # adesso applichiamo -log10 per "invertire" i pvalue bassi in pesi alti e viceversa 
          lg    <- -log10(norm1)
          
          ## normalizziamo lg tra 0 e 1, calcolo parametri ( lg + eps_log) * k_log = p_norm2
          # ( MAX_lg  + eps_log ) * k_log = 1
          # ( min(lg) + eps_log ) * k_log = min(lg)
          MAX_lg <- -log10(signifCorr)
          
          # risolvendo il sistema si ottiene 
          eps_log <- ( min(lg) - MAX_lg * min(lg) ) / ( min(lg) - 1 )
          k_log   <- ( min(lg) - 1 ) / ( min(lg) - MAX_lg)
          norm2   <- (lg + eps_log) * k_log
          
          #moltiplichiamo la correlazione per il pnorm3
          flattCorr$weights <- flattCorr$cor
          flattCorr$weights[ind_notsig] <- flattCorr$cor[ind_notsig]*norm2^3
          
        } else {
          flattCorr$weights <- flattCorr$cor
          flattCorr$weights[!ind_sig] <- 0
        }
        
        Sign <- flattCorr[ind_sig,] #storing significative ones
        Sign <- Sign[order(Sign$cor_features),]
        
        ftw <- flattCorr[,c(1,2,8)]
        ftw[,3] <- 10 - 10 *abs(ftw[,3])
        gcor <- graph_from_edgelist(as.matrix(ftw[,1:2]),directed = F)
        edge.attributes(gcor)$weights <- ftw[,3]
        
        g.f <- intersection(gcor,ginter) ## grafo correlazioni e interazioni, con pesi per cytoscape 10 * 10 - w
        g.fel <- cbind(as_edgelist(g.f),edge.attributes(g.f)$weights)
        
        return(list(gCor.groups   = flattCorr,
                    gInter.groups = g.f))
        
      })
      gCOR.groups   <- lapply(gCOR.groups.l, function(x) x$gCor.groups)
    } else {
      gCOR.groups.l <- list(gCor.groups   = gCOR.groups,
                            gInter.groups = NULL)
    }
    
    gCOR.categories <- lapply(gNULL.categories, function(x){
      ## aggiungere trycatch con mesasggio di funzione gcor.groups standard
      flattCateg <- lapply(gCOR.groups, function(y){
        gcor.el <- y[match(x$el$cor_features,y$cor_features),]
        if (is.null(correctionCorr)){
          gcor.el$p.adj <- gcor.el$p
        } else {
          gcor.el$p.adj <- p.adjust(gcor.el$p, correctionCorr)
        }
        ind_sig <- gcor.el$p.adj <= signifCorr
        
        if (compute_weights){
          ### normalizzazione per togliere gli 1 dal pvalue, per rendere possibile il logaritmo
          # L'obiettivo è trasformare l'intervallo del pvalue da (Significativià ; 1 ] a (0 ; 1)
          
          # Il nuovo massimo del pvalue sarà la media tra 1 e il secondo massimo
          if (max(gcor.el$p.adj)==1){
            new_max <- mean(c(1,max(gcor.el$p.adj[-which(gcor.el$p.adj==1)])))
          } else {
            new_max <- max(gcor.el$p.adj)
          }
          
          # calcolo parametri normalizzazione: (pvalue + epsilon ) * k = p_norm1
          # ( signifCorr + eps ) * k = signifCorr
          # ( max(pvalue)  + eps ) * k = new_max
          
          # risolvendo questo sistema si ottiene
          eps <- ( new_max - max(gcor.el$p.adj) ) * signifCorr / ( signifCorr - new_max )
          k   <- ( signifCorr - new_max ) / ( signifCorr - max(gcor.el$p.adj) )
          
          # applichiamo la normalizzazione solo ai pvalue non significativi
          ind_notsig <- which(gcor.el$p.adj > signifCorr)
          norm1 <- (gcor.el$p.adj[ind_notsig] + eps) * k
          
          # adesso applichiamo -log10 per "invertire" i pvalue bassi in pesi alti e viceversa 
          lg    <- -log10(norm1)
          
          ## normalizziamo lg tra 0 e 1, calcolo parametri ( lg + eps_log) * k_log = p_norm2
          # ( MAX_lg  + eps_log ) * k_log = 1
          # ( min(lg) + eps_log ) * k_log = min(lg)
          MAX_lg <- -log10(signifCorr)
          
          # risolvendo il sistema si ottiene 
          eps_log <- ( min(lg) - MAX_lg * min(lg) ) / ( min(lg) - 1 )
          k_log   <- ( min(lg) - 1 ) / ( min(lg) - MAX_lg)
          norm2   <- (lg + eps_log) * k_log
          
          #moltiplichiamo la correlazione per il pnorm3
          gcor.el$weights <- gcor.el$cor
          gcor.el$weights[ind_notsig] <- gcor.el$cor[ind_notsig]*norm2^3
          
        } else {
          gcor.el$weights <- gcor.el$cor
          gcor.el$weights[!ind_sig] <- 0
        }
        
        Sign <- gcor.el[ind_sig,] #storing significative ones
        Sign <- Sign[order(Sign$cor_features),]
        
        return(list(All_correlation  = gcor.el,
                    Sign_correlation = Sign))
      })
      gCor.categ.all <- lapply(flattCateg, function(y){
        g.Cor <- graph_from_edgelist(as.matrix(y$All_correlation[,1:2]), directed = F)
        edge.attributes(g.Cor)$weight <- y$All_correlation$weights
        return(g.Cor)
      })
      gCor.categ.sign <- lapply(flattCateg, function(y){
        g.Cor <- graph_from_edgelist(as.matrix(y$Sign_correlation[,1:2]), directed = F)
        edge.attributes(g.Cor)$weight <- y$Sign_correlation$weights
        return(g.Cor)
      })
      
      Other.PPI <- x$length.interactome - nrow(x$el)
      sd_groups <- sapply(flattCateg, function(x){
        sd(c(x$All_correlation$weights,rep(1,Other.PPI)))/2
      })
      return(list(flattCorr = flattCateg,
                  gCor.categ.all = gCor.categ.all,
                  gCor.categ.sign = gCor.categ.sign,
                  Sd_groups = sd_groups,
                  Other.PPI = Other.PPI))
    })
    time_employed <- format_time_difference(difftime(Sys.time(),
                                                     start.time.correlation,units = "s"))
    message("Correlations computed in ", time_employed,sep = " ")
  }
  
  
  names(categories) <- categories
  
  results <- lapply(categories, function(x){
    
    return(CoPPI(annotations = prot.annotation[[x]], 
                 correlations = gCOR.categories[[x]],
                 perc.ned.nppi = perc.ned.nppi,
                 perc.sign = perc.sign,
                 significance.coppi = significance.coppi,
                 correction.coppi = correction.coppi,
                 min_edges = min_edges, 
                 max_edges = max_edges,
                 names_of_groups = names_of_groups))
  })
  
  dir.create(name_analysis)
  setwd(name_analysis)
  tryCatch({
    lapply(results,save_results, name_analysis)
    cytoscapePing()
    closeSession(save.before.closing = F)
    all_excels <- dir()[grep(".xlsx",dir())]
    names(all_excels) <- all_excels
    
    lapply(all_excels,Load2Cytoscape,
           filter_similarity =  threshold_edge_cyto)
    setCytoStyle()
    saveSession(str_c(name_analysis,"_",threshold_edge_cyto,".cys"))
    message("Output saved in Cytoscape")
  }, error = function(err){
    message("Warning in Cytoscape:")
    message(err)
    message("\n")
    setwd("..")
  })
  
  Centralities <- tryCatch({
    Centralities <- Compute_centralities(gCOR.groups = gCOR.groups,
                                         interactome = summary_interactome$interactome,
                                         signifCorr = signifCorr,
                                         correctionCorr = correctionCorr,
                                         compute_weights = compute_weights,
                                         thr_degree = thr_degree, 
                                         thr_corr = thr_corr)
    tryCatch({
      dir.create("Centralities")
      setwd("Centralities")
      save_excel.centralities(Centralities, name_analysis)
      setwd("..")},
      error = function(err){
        message("Warning in Centralities:")
        message(err)
        message("\n")
        setwd("..")
        message("The centralities are computed but not saved")
      })
    },
    error = function(err){
      message("Warning in centralities:")
      message(err)
      message("\n")
    })
  
  
  CoPPI_results <- list(prot.annotation = prot.annotation,
                        gCOR.groups = gCOR.groups.l,
                        gCOR.categories = gCOR.categories,
                        resultsCoPPI = results,
                        centralities = Centralities,
                        parameters = list(name_analysis = name_analysis,
                                          signifCorr = signifCorr, 
                                          correctionCorr = correctionCorr,
                                          min_n.corr = min_n.corr,
                                          perc.ned.nppi = perc.ned.nppi,
                                          perc.sign = perc.sign,
                                          significance.coppi = significance.coppi,
                                          correction.coppi = correction.coppi,
                                          corr.test = corr.test,
                                          compute_weights = compute_weights,
                                          min_edges = min_edges,
                                          max_edges = max_edges,
                                          categories = categories,
                                          threshold_edge_cyto = threshold_edge_cyto,
                                          thr_degree = thr_degree,
                                          thr_corr = thr_corr))
  save(CoPPI_results, file = str_c(CoPPI_results$parameters$name_analysis,".RDa"))
  setwd("..")
  return(CoPPI_results)
}
