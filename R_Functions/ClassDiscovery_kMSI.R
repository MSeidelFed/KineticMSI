




ClassDiscovery_kMSI <- function(FilesPath,
                                method.dist = "abscor",
                                nboot = 100,
                                return_SigClustHist = T,
                                fun_to_clust = c("Enrichment", "Spatial"),
                                pattern = "MeanEnrichment.csv",
                                alpha = 0.9) {
  
 
  
  original_alpha = alpha
  
  ## functions needed for the main
  
  ## distribution ggplot function
  
  Class_Distribution <- function(in_mat,
                                 Treatments, 
                                 plots = NULL, 
                                 returnTable = F) {
    
    mydata4.1 <- as.data.frame(t(in_mat))
    colnames(mydata4.1) <- seq(1, dim(mydata4.1)[2])
    
    mydata4.2 <- cbind(Treatments, mydata4.1)
    colnames(mydata4.2) <- c("X1_1", seq(1,dim(mydata4.1)[2]))
    
    mydata4.3 <- melt(mydata4.2, id.vars= seq(1), measure.vars= seq(2, (dim(mydata4.1)[2]+1)))
    
    mydata4.4 <- mydata4.3[order(mydata4.3$X1_1),]
    
    if(plots == T) {
      
      print(ggplot(mydata4.4, aes(x = as.numeric(value), y = X1_1)) + 
              xlim(summary(as.numeric(mydata4.4$value))["Min."] - summary(as.numeric(mydata4.4$value))["Mean"],
                   summary(as.numeric(mydata4.4$value))["Mean"] + summary(as.numeric(mydata4.4$value))["Mean"]) +
              ggtitle ("Test") +
              geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
              theme_ridges(font_size = 10, grid = TRUE) +
              theme(axis.title.y = element_blank()))
      
      
      print(ggplot(mydata4.4, aes(x = as.numeric(value), fill = X1_1)) +
              geom_density(alpha = 0.2, position = "identity") +
              xlim(summary(as.numeric(mydata4.4$value))["Min."] - summary(as.numeric(mydata4.4$value))["Mean"],
                   summary(as.numeric(mydata4.4$value))["Mean"] + summary(as.numeric(mydata4.4$value))["Mean"]) +
              ggtitle ("Test"))
      
    }
    
    
    if(returnTable == T) {
      
      return(mydata4.4)
      
    }
    
  }
  
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
  
  ### all the rep files need to be in the same directory
  
  #### finding a significant number of partitions in joined datasets across metabolites and samples
  
  pdf(file = "Bootstrapped_Dendrograms.pdf")
  
  par(mfrow = c(1,1))
  
  list_out <- list()
  
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
    
    if(show_condition(pvclust(runner_WO_zeros, nboot = 10,
                              quiet = T, method.dist = method.dist)[1]) == "error") {
      # getting rid of metabolites that produce an error in the pvclust function
      
      cat("...\n")
      cat(paste0("Error: ", lipid_Nr[i], " - ", i, "\n"))
      cat("...\n")
      
    } else {
      
      cat(paste0("Successful: ", lipid_Nr[i]," - ", i, "\n"))
      
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
        
        #### plotting distributions
        
        full_file_names <- list2df(strsplit(apply(list2df(strsplit(rownames(runner_WO_zeros), "/"))[5,], 2, as.character), "_"))
        
        rep_name <- paste0(apply(full_file_names[2,], 2, as.character),
                           "_",
                           apply(full_file_names[3,], 2, as.character))
        
        Class_Distribution(in_mat = t(runner_WO_zeros), Treatments = as.factor(rep_name), plots = T, returnTable = F)
        
        #plot(density(runner_WO_zeros), main = lipid_Nr[i], xlab = "Enrichment Proportion")
        
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
  
  dev.off()
  
  ### heavy-tailed distribution of cluster number
  
  if (return_SigClustHist == T) {
    
    png(filename = "SigClustHist.png")
    plot(density(na.omit(as.numeric(sig_clust_Nr))),
         xlab = "Significant Cluster No.",
         main = "Cluster No. across features")
    dev.off()
    
  }
  
  return(list_out)
   
}

