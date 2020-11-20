





GeneralExpOverview <- function(ClassDiscovery_List,
                               ClassComparison_list,
                               ControlSample = "_WT",
                               returnHeatmaps = T,
                               factorVector) {
  
  
  ## functions needed
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  VolcanoPlots <- function(in_mat,
                           ControlSample,
                           FactorCols = F,
                           factor_out,
                           omics_test,
                           returnPlotsPNGs = F,
                           outID) {
    
    if (FactorCols == T) {
      
      factor_out <- colnames(in_mat)
      
    } else {
      
      colnames(in_mat) <- factor_out
      
    }
    
    #print(factor_out)
    
    in_omics <- omics_test
    
    volcanoes_treatments <- unique(factor_out)
    
    
    ### non-variable data
    
    test <- unique(colnames(in_mat[,grep(ControlSample, colnames(in_mat))]))
    
    if (length(test) > 1) {
      
      #### more than one cluster
      
      for (i in 1:length(test)) {
        
        test_names <- strsplit(test[i], ControlSample)[[1]]
        
        out_factor <- c()
        
        for (i in 1:length(test_names)) {
          
          runner_test <- strsplit(test_names[i], split = "")[[1]]
          
          if (length(runner_test) > 0) {
            
            out_factor <- paste0(out_factor, paste(runner_test, collapse = ""))
            
          }
          
        }
        
        control <- as.numeric(rowMeans(in_mat[,grep(out_factor,
                                                    colnames(in_mat))]))
        
        list_plots <- list()
        
        rm_ctrl <- grep(paste0(ControlSample,out_factor), volcanoes_treatments)
        
        for (i in 1:length(volcanoes_treatments[-rm_ctrl])) {
          
          runner <- volcanoes_treatments[-rm_ctrl][i]
          
          runner_means <- as.numeric(rowMeans(in_mat[,grep(runner,
                                                           colnames(in_mat))]))
          
          Abscissa <- log2(runner_means/control)
          
          Ordinate <- -log10(as.numeric(in_omics[,intersect(grep(out_factor,
                                                                 colnames(in_omics)),
                                                            grep(ControlSample,
                                                                 colnames(in_omics)))]))
          
          data <- cbind.data.frame(Abscissa, Ordinate)
          
          x <- ggplot(data, aes(Abscissa,
                                Ordinate,
                                label = apply(list2df(strsplit(rownames(in_mat),
                                                               "_"))[1,],2,as.character))) +
            geom_text_repel() +
            geom_point(color = 'red') +
            theme_classic(base_size = 12) + 
            geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red") +
            ggtitle(paste0(outID, runner)) +
            labs(x = paste0("Log2(", runner, "/", ControlSample, ")"),
                 y = paste0("-log10(Padj-Values [", runner, "])"))
          
          
          list_plots[[i]] <- x
          
          if (returnPlotsPNGs == T) {
            
            ggsave(list_plots[[i]], file=paste0("plot_", outID, runner,".png"),
                   width = 44.45, height = 27.78, units = "cm", dpi=300)
            
          }
        }
        
        return(list_plots)
        
      }
      
    } else {
      
      #### one cluster (working)
      
      control <- as.numeric(rowMeans(in_mat[,grep(ControlSample,
                                                  colnames(in_mat))]))
      
      list_plots <- list()
      
      rm_ctrl <- grep(ControlSample, volcanoes_treatments)
      
      for (i in 1:length(volcanoes_treatments[-rm_ctrl])) {
        
        runner <- volcanoes_treatments[-rm_ctrl][i]
        
        runner_means <- as.numeric(rowMeans(in_mat[,grep(runner,
                                                         colnames(in_mat))]))
        
        Abscissa <- log2(runner_means/control)
        
        Ordinate <- -log10(as.numeric(in_omics[,grep("factor_glm", colnames(in_omics))]))
        
        data <- cbind.data.frame(Abscissa, Ordinate)
        
        x <- ggplot(data, aes(Abscissa,
                              Ordinate,
                              label = apply(list2df(strsplit(rownames(in_mat),
                                                             "_"))[1,],2,as.character))) +
          geom_text_repel() +
          geom_point(color = 'red') +
          theme_classic(base_size = 12) + 
          geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red") +
          ggtitle(runner) +
          labs(x = paste0("Log2(", runner, "/", ControlSample, ")"),
               y = paste0("-log10(Padj-Values [", runner, "])"))
        
        
        list_plots[[i]] <- x
        
        if (returnPlotsPNGs == T) {
          
          ggsave(list_plots[[i]], file=paste0("plot_", outID, runner,".png"),
                 width = 44.45, height = 27.78, units = "cm", dpi=300)
          
        }
      }
      
      return(list_plots)
      
      #### here closes loop of one cluster if
      
    }
    
  }
  
  
  
  ## main
  
  ### assembling in_mat from the out_list of class comparison function
  
  clustNo <- c()
  
  for (i in 1:length(ClassDiscovery_List)) {
    
    if (!is.null(ClassDiscovery_List[[i]])) {
      
      clustNo[i] <- dim(ClassDiscovery_List[[i]])[2]
      
    }
  }
  
  max_ClustNo <- max(na.omit(clustNo))
  
  ClustIndexes <- unique(na.omit(clustNo))
  
  #### for loop to assemble matrices from features with common numbers of k partitions
  
  names_list_out <- c()
  
  list_out <- list()
  
  list_omics_out <- list()
  
  for (i in 1:length(ClustIndexes)) {
    
    ### iterating through the different ClustNos.
    
    runner_list <- ClassDiscovery_List[which(clustNo == ClustIndexes[i])]
    
    ### producing in_mat
    
    out_mat <- matrix(NA,
                      nrow = length(melt(runner_list[1])$value),
                      ncol = length(runner_list))
    
    ### producing omics_test
    
    mat_width <- length(ClassComparison_list[which(names(ClassComparison_list) == names(runner_list)[1])][[1]])
    
    out_omics_mat <- matrix(NA,
                            nrow = length(runner_list),
                            ncol = mat_width)
    
    ### producing vector with feature names
    
    namei <- c()
    
    for (j in 1:length(runner_list)) {
      
      out_mat[,j] <- melt(runner_list[j])$value
      
      runner_omics <- ClassComparison_list[which(names(ClassComparison_list) == names(runner_list)[j])][[1]]
      
      out_omics_mat[j,] <- runner_omics
      
      namei[j] <- unique(melt(runner_list[j])$L1)
      
    }
    
    colnames(out_mat) <- namei
    rownames(out_mat) <- paste0(rep(factorVector, ClustIndexes[i]),
                                melt(runner_list[j])$Var2)
    
    colnames(out_omics_mat) <- names(runner_omics)
    rownames(out_omics_mat) <- namei
    
    names_list_out[i] <- paste0("ClustNo_", ClustIndexes[i])
    
    list_out[[i]] <- out_mat
    
    list_omics_out[[i]] <- out_omics_mat
    
    ## returning heatmap
    
    if (returnHeatmaps == T) {
      
      png(filename = paste0(names_list_out[i], ".png"))
      heatmap(t(out_mat))
      dev.off()
      
    }
    
  }
  
  names(list_out) <- names_list_out
  
  names(list_omics_out) <- names_list_out
  
  Hawaii <- list()
  
  for (i in 1:length(list_out)) {
    
    test_volcanoes <- VolcanoPlots(in_mat = t(list_out[[i]]),
                                   FactorCols = T, 
                                   factor_out = factorVector,
                                   ControlSample = ControlSample,
                                   omics_test = list_omics_out[[i]],
                                   returnPlotsPNGs = T,
                                   outID = names(list_out)[i])
    
    Hawaii <- c(Hawaii, test_volcanoes)
    
  }
  
  return(Hawaii)
  
}


