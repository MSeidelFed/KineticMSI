




ClassDiscovery_kMSI <- function(FilesPath,
                                method.dist = "abscor",
                                nboot = 100,
                                return_SigClustHist = T,
                                fun_to_clust = c("Enrichment", "Spatial"),
                                logiTransformation = F,
                                pattern = "MeanEnrichment.csv",
                                alpha = 0.9,
                                returnObject = c(NULL, "ClassComparisonInput", "minDataset"),
                                CompareSampledSet = T) {
  
 
  
  original_alpha = alpha
  
  ## functions needed for the main

  
  #### Tukey function
  
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  
  #### handling errors
  
  show_condition <- function(code) {
    tryCatch(code,
             error = function(c) "error",
             message = function(c) "message"
    )
  }
  
  ### list2df
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  } 
  
  ### EnrichmentSorted_RandomSampling_MSI
  
  EnrichmentSorted_RandomSampling_MSI <- function(FilesPath,
                                                  pattern = "MeanEnrichment.csv",
                                                  feature_choice = 0) {
    #### choosing a lipid
    
    reps <- list.files(path = FilesPath, pattern = pattern, all.files = FALSE,
                       full.names = T, recursive = T,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
    if (feature_choice == 0) {
      
      file_test <- read.csv(file = reps[1], header = T, row.names = 1)
      
      lipid_Nr <- unique(rownames(file_test))
      
      print(lipid_Nr)
      
      lipid_choice <- readline(prompt="Enter index of lipid of interest: ")
      
    } else {lipid_choice = feature_choice}
    
    #### defining the minimum set length
    
    out_vector <- c()
    
    for (i in 1:length(reps)) {
      
      file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
      
      out_vector[i] <- dim(file_runner)[2]
      
    }
    
    min_dataset <- min(out_vector)
    
    #### sampling randomly but sorting by enrichment percentage
    
    out_single_lipid <- matrix(NA, length(reps), min_dataset)
    
    for (i in 1:length(reps)) {
      
      file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
      
      out_single_lipid[i,] <- sort(as.numeric(sample(file_runner[as.numeric(lipid_choice),],
                                                     min_dataset,
                                                     replace = F))
                                   , T)
      
    }
    
    rownames(out_single_lipid) <- reps
    
    return(out_single_lipid)
    
  }
  
  
  ### SpatialRandomSampling_MSI
  
  SpatialRandomSampling_MSI <- function(FilesPath,
                                        pattern = "MeanEnrichment.csv",
                                        feature_choice = 0) {
    
    #### choosing a lipid
    
    reps <- list.files(path = FilesPath, pattern = pattern, all.files = FALSE,
                       full.names = T, recursive = T,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
    if (feature_choice == 0) {
      
      file_test <- read.csv(file = reps[1], header = T, row.names = 1)
      
      lipid_Nr <- unique(rownames(file_test))
      
      print(lipid_Nr)
      
      lipid_choice <- readline(prompt="Enter index of lipid of interest: ")
      
    } else  { lipid_choice = feature_choice }
    
    #### defining the minimum set length
    
    out_vector <- c()
    
    for (i in 1:length(reps)) {
      
      file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
      
      out_vector[i] <- dim(file_runner)[2]
      
    }
    
    min_dataset <- min(out_vector)
    
    #### sampling randomly but conserving spatial relationships among pixels
    
    out_single_lipid <- matrix(NA, length(reps), min_dataset)
    
    for (i in 1:length(reps)) {
      
      file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
      
      lipid_runner <- as.matrix(sample(file_runner[as.numeric(lipid_choice),],
                                       min_dataset,
                                       replace = F))
      
      df2sort <- as.data.frame(t(rbind(apply(list2df(strsplit(colnames(lipid_runner),
                                                              "\\."))[2,],
                                             2,
                                             as.numeric),
                                       lipid_runner)))
      
      colnames(df2sort) <- c("pixelNr", "Enr")
      
      SortedEnrichment <- df2sort[order(df2sort$pixelNr),2]
      
      out_single_lipid[i,] <- SortedEnrichment
      
    }
    
    rownames(out_single_lipid) <- reps
    
    return(out_single_lipid)
    
  }
  
  
  
  
  ## Main
  
  ### grabbing data
  
  reps <- list.files(path = FilesPath, pattern = pattern, all.files = FALSE,
                     full.names = T, recursive = T,
                     ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  ### getting the common lipid_Nr for all datasets
  
  lipid_Nr <- list()
  
  lipid_Nr_vec <- c()
  
  for (i in 1:length(reps)) {
    
    file_test <- read.csv(file = reps[i], header = T, row.names = 1)
    
    lipid_Nr[[i]] <- unique(rownames(file_test))
    
    lipid_Nr_vec <- c(lipid_Nr_vec, unique(rownames(file_test)))
    
  }
  
  lipid_Nr = Reduce(intersect, lipid_Nr)
  
  lipid_Nr_vec = unique(lipid_Nr_vec)
  
  if (length(which((lipid_Nr == lipid_Nr_vec) == F)) > 0) {
    
    for (i in 1:length(reps)) {
      
      file_test <- read.csv(file = reps[i], header = T, row.names = 1)
      
      out_dir <- paste0(strsplit(reps[i],
                                 "/")[[1]][-c(length(strsplit(reps[i],
                                                              "/")[[1]]))], collapse = "/")
      
      file_name <- paste0(strsplit(strsplit(reps[i],
                                            "/")[[1]][length(strsplit(reps[i],
                                                                      "/")[[1]])],
                                   "\\.")[[1]][1], "_SharedFeatures",
                          ".csv")
      
      write.csv(file_test[lipid_Nr,], file = paste0(out_dir, "/", file_name),
                quote = F)
      
    }
    
    pattern = "_SharedFeatures.csv"
    
  }
  
  print(pattern)
  
  
  
  if (CompareSampledSet == T) {
    
    ## getting rowmeans of all files to compare to the sampled ones
    
    ### full matrices
    
    CompleteSet_Means <- c()
    CompleteSet_SD <- c()
    
    for (i in 1:length(reps)) {
      
      runner_file <- read.csv(reps[i], header = T, row.names = 1)
      colnames(runner_file) <- NULL
      
      CompleteSet_Means <- c(CompleteSet_Means, rowMeans(runner_file))
      CompleteSet_SD <- c(CompleteSet_SD, rowSds(as.matrix(runner_file)))
      
    }
    
    ### sampled matrices
    
    out_vector <- c()
    
    for (i in 1:length(reps)) {
      
      file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
      
      out_vector[i] <- dim(file_runner)[2]
      
    }
    
    min_dataset <- min(out_vector)
    
    out_mat <- matrix(NA, nrow = 0, ncol = min_dataset)
    
    for (i in 1:length(reps)) {
      
      runner_file <- read.csv(reps[i], header = T, row.names = 1)
      
      runner_sampled <- as.matrix(sample(runner_file, min_dataset, replace = F))
      
      out_mat <- rbind(out_mat, runner_sampled)
      
    }
    
    MinSet_Means <- rowMeans(out_mat)
    MinSet_SD <- rowSds(out_mat)
    
    pdf(file = "ComparedSampledSet.pdf")
    
    par(mfrow = c(2,1))
    
    plot(CompleteSet_Means/MinSet_Means, main = "Means")
    plot(CompleteSet_SD/MinSet_SD, main = "SD")
    
    dev.off()
    
  }
  
  
  ### all the rep files need to be in the same directory
  
  #### finding a significant number of partitions in joined datasets across metabolites and samples
  
  pdf(file = "Bootstrapped_Dendrograms.pdf")
  
  par(mfrow = c(1,1))
  
  list_out <- list()
  
  list_out2 <- list()
  
  sig_clust_Nr <- c()
  
  lipid_names <- c()
  
  for (i in 1:length(lipid_Nr)) {
    
    alpha = original_alpha
    
    if (fun_to_clust == "Enrichment") {
      
      runner <- EnrichmentSorted_RandomSampling_MSI(FilesPath = FilesPath,
                                                    pattern = pattern, feature_choice = i)
      
    } else if (fun_to_clust == "Spatial") {
      
      runner <- SpatialRandomSampling_MSI(FilesPath = FilesPath,
                                          pattern = pattern, feature_choice = i)
      
    } else {
      
      cat("Error: select an appropriate sorting function")
      
    }
    
    colnames(runner) <- paste0(rep("X", dim(runner)[2]), c(1:dim(runner)[2]))
    
    runner_WO_zeros <- runner[,as.numeric(colMeans(runner)) != 0]
    
    ## transforming to logit the successful if set to TRUE
    
    if (logiTransformation == T) {
      
      logit_transformation <- function(x){log(x /(1 - x))}
      
      ### Replace zero with a very small value so it can fit in a model
      
      out_mat_DG_wzero <- runner_WO_zeros
      
      logit_mat_DG <- matrix(NA, nrow = dim(out_mat_DG_wzero)[1], ncol = dim(out_mat_DG_wzero)[2])
      
      for (k in 1:dim(out_mat_DG_wzero)[1]) {
        
        for (l in 1:dim(out_mat_DG_wzero)[2]) {
          
          if (out_mat_DG_wzero[k,l] > 0) {
            
            logit_mat_DG[k,l] <- logit_transformation(out_mat_DG_wzero[k,l])
            
          } else {
            
            logit_mat_DG[k,l] <- 0
            
          }
          
        }
        
      }
      
      rownames(logit_mat_DG) <- rownames(runner_WO_zeros)
      colnames(logit_mat_DG) <- colnames(runner_WO_zeros)
      
      runner_WO_zeros = logit_mat_DG
      
      runner_WO_zeros[is.na(runner_WO_zeros)] <- 0
      
      runner_WO_zeros[logit_mat_DG == "NaN"] <- 0
      
      runner_WO_zeros[logit_mat_DG == "Inf"] <- 0
      
      runner_WO_zeros[logit_mat_DG == "-Inf"] <- 0
      
    }
    
    
    if(show_condition(pvclust(runner_WO_zeros, nboot = 10,
                              quiet = T, method.dist = method.dist)[1]) == "error") {
      # getting rid of metabolites that produce an error in the pvclust function
      
      cat("...\n")
      cat(paste0("Error: ", lipid_Nr[i], " - ", i, "\n"))
      cat("...\n")
      
    } else {
      
      cat(paste0("Successful: ", lipid_Nr[i]," - ", i, "\n"))
    
      ## bootstrapping the clustering on the final object
      
      HCA_boot_lipid <- pvclust(runner_WO_zeros,
                                method.hclust = "average",
                                method.dist = method.dist,
                                nboot = nboot, 
                                quiet = T)
      
      if (dim(runner_WO_zeros)[2] > 0.2 * dim(runner)[2]) { 
        # plot works on matrices with 20% or more columns with values different than zero
        
        table_lipids <- pvpick(HCA_boot_lipid, alpha = alpha, pv="au")
        
        #### changing the alpha in dendrograms with one significant cluster
        
        count = 0
        
        alpha_out <- alpha
        
        if (length(table_lipids$edges) == 1) {
          
          cat("k before relaxation: ", length(table_lipids$edges), "\n")
          
        }
        
        ## lower limit until AU = .77
        
        while (length(table_lipids$edges) == 1) {
          
          count = count + 1
          
          if (count < 7) {
            
            alpha = alpha-((1-alpha)/8)
            
            table_lipids <- pvpick(HCA_boot_lipid, alpha = alpha, pv="au")
            
            alpha_out = alpha
            
          } else {
            
            break
            
          }
          
        }
        
        ## upper limit until AU = .99
        
        while (length(table_lipids$edges) == 1) {
          
          count = count + 1
          
          if (count < 30) {
            
            alpha = alpha+((1-alpha)/8)
            
            table_lipids <- pvpick(HCA_boot_lipid, alpha = alpha, pv="au")
            
            alpha_out = alpha
            
          } else {
            
            break
            
          }
          
        }
        
        cat("Final k : ", length(table_lipids$edges), "\n")
        
        cat("Final AU-P: ", alpha_out, "\n")
        
        #### plotting dendrograms with significant clusters
        
        plot(x = HCA_boot_lipid, print.pv = "au", print.num = T, cex = 1, lwd = 2, main = lipid_Nr[i])
        
        pvrect(HCA_boot_lipid, alpha = alpha, pv="au")
        
        #### building matrices using significant clusters as factors
        
        if (length(table_lipids$edges) < 10 & length(table_lipids$edges) > 0) {
          ##### Metabolites with more than 10 parental partitions are neglected
          
          sig_clust_Nr[i] <- length(table_lipids$edges)
          
          lipid_names[i] <- lipid_Nr[i]
          
          #### mining out the cluster to build means
          
          ##### establishing the maximum cluster size
          
          out_vector <- c()
          
          for (k in 1:length(table_lipids$clusters)) {
            
            out_vector[k] <- length(table_lipids$clusters[[k]])
            
          }
          
          max_cluster_length <- max(out_vector)
          
          ##### grabbing clusters with colnames of clustered pixels into a data frame
          
          clusters_p95 <- sapply(1:length(as.matrix(table_lipids$clusters)),function(j) as.matrix(table_lipids$clusters)[[j]][1:max_cluster_length])
          
          #### subsetting the runner_WO_zeros object and getting the means
          
          out_mat <- matrix(NA, nrow = dim(runner_WO_zeros)[1], ncol = sig_clust_Nr[i])
          
          for (n in 1:dim(clusters_p95)[2]) {
            
            out_mat[,n] <- rowMeans(runner_WO_zeros[,na.omit(clusters_p95[,n])])
            
          }
          
          colnames(out_mat) <- table_lipids$edges
          
          list_out[[i]] <- out_mat
          list_out2[[i]]  <- runner
          
          cat("...\n")
          
        } else {
          
          cat(paste0("ommited - too few or too many significant K classes ", lipid_Nr[i]," - ", i, "\n"))
          cat("...\n")
          
        }
        
      } else {
        
        cat(paste0("Error in dendrogram plot (too many Zeros): ", lipid_Nr[i]," - ", i, "\n"))
        cat("...\n") 
        
      }}
  }
  
  names(list_out) <- lipid_names
  names(list_out2) <- lipid_names
  
  dev.off()
  
  ### heavy-tailed distribution of cluster number
  
  if (return_SigClustHist == T) {
    
    png(filename = "SigClustHist.png")
    plot(density(na.omit(as.numeric(sig_clust_Nr))),
         xlab = "Significant Cluster No.",
         main = "Cluster No. across features")
    dev.off()
    
  }
  
  if(returnObject == "ClassComparisonInput") {
    
    return(list_out)
    
  } else if (returnObject == "minDataset") {
    
    return(list_out2)
    
  } else if (is.null(returnObject)) {
    
    return(list_out)
    
  }
   
}

