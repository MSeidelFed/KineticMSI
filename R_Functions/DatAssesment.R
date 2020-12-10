










DatAssesment <- function(FilesPath,
                             pattern = "MeanEnrichment.csv",
                             CompareSampledSet = T,
                             returnObject = c(NULL, "ClassComparisonInput", "minDataset"),
                             factorVector,
                             logiTransformation = F) {
  
  
   ### list2df
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  } 
  
  ## distribution ggplot function
  
  Class_Distribution <- function(in_mat,
                                 Treatments, 
                                 plots = NULL, 
                                 returnTable = F,
                                 factorVector) {
    
    mydata4.1 <- as.data.frame(t(in_mat))
    colnames(mydata4.1) <- seq(1, dim(mydata4.1)[2])
    
    mydata4.2 <- cbind(Treatments, mydata4.1)
    colnames(mydata4.2) <- c("X1_1", seq(1,dim(mydata4.1)[2]))
    
    mydata4.3 <- melt(mydata4.2, id.vars= seq(1), measure.vars= seq(2, (dim(mydata4.1)[2]+1)))
    
    mydata4.4 <- mydata4.3[order(mydata4.3$X1_1),]
    
    unique_factors <- unique(factorVector)
    
    out_vec <- c()
    
    for (i in 1:length(unique_factors)) {
      
      runner <- grep(unique(factorVector)[i], mydata4.4$X1_1)
      
      out_vec <- c(out_vec, rep(unique_factors[i], length(runner)))
      
    }
    
    if(plots == T) {
      
      print(ggplot(mydata4.4, aes(x = as.numeric(value), y = X1_1, color = out_vec)) + 
              xlim(summary(as.numeric(mydata4.4$value))["Min."] - summary(as.numeric(mydata4.4$value))["Mean"],
                   summary(as.numeric(mydata4.4$value))["Mean"] + summary(as.numeric(mydata4.4$value))["Mean"]) +
              ggtitle ("Test") +
              geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
              theme_ridges(font_size = 10, grid = TRUE) +
              theme(axis.title.y = element_blank()))
      
      
      print(ggplot(mydata4.4, aes(x = as.numeric(value), fill = X1_1, color = out_vec)) +
              geom_density(alpha = 0.2, position = "identity") +
              xlim(summary(as.numeric(mydata4.4$value))["Min."] - summary(as.numeric(mydata4.4$value))["Mean"],
                   summary(as.numeric(mydata4.4$value))["Mean"] + summary(as.numeric(mydata4.4$value))["Mean"]) +
              ggtitle ("Test"))
      
    }
    
    
    if(returnTable == T) {
      
      return(mydata4.4)
      
    }
    
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
  
  
  pdf(file = "Data_assesment.pdf")
  
  par(mfrow = c(1,1))
  
  list_out <- list()
  
  list_out2 <- list()
  
  sig_clust_Nr <- c()
  
  lipid_names <- c()
  
  for (i in 1:length(lipid_Nr)) {
    
    cat("...\n")
    cat(paste0("Processing: ", lipid_Nr[i]," - ", i, "\n"))
    cat("...\n")
    
    runner <- SpatialRandomSampling_MSI(FilesPath = FilesPath,
                                          pattern = pattern, feature_choice = i)
      
    colnames(runner) <- paste0(rep("X", dim(runner)[2]), c(1:dim(runner)[2]))
    
    runner_WO_zeros <- runner
    
    ## transforming to logit if set to TRUE
    
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
    
    
    list_out[[i]] <- runner_WO_zeros
    lipid_names[i] <- lipid_Nr[i]
    
    #### plotting distributions
    
    full_file_names <- list2df(strsplit(apply(list2df(strsplit(rownames(runner_WO_zeros), "/"))[5,], 2, as.character), "_"))
    
    rep_name <- paste0(apply(full_file_names[2,], 2, as.character),
                       "_",
                       apply(full_file_names[3,], 2, as.character))
    
    Class_Distribution(in_mat = t(runner_WO_zeros), Treatments = as.factor(rep_name),
                       plots = T, returnTable = F, factorVector = factorVector)
    
    #plot(density(runner_WO_zeros), main = lipid_Nr[i], xlab = "Enrichment Proportion")
    
    
    ## second output (class comparison input & distribution test)
    
    out_mat2 <- as.matrix(rowMeans(runner_WO_zeros))
    
    list_out2[[i]] <- out_mat2
    
    plot_name <- testing_distributions(Distribution_test_mat = out_mat2)
    
    plot_mat <- as.matrix(cbind(out_mat2,out_mat2,out_mat2,out_mat2,out_mat2,
                                out_mat2,out_mat2,out_mat2,out_mat2,out_mat2,
                                out_mat2,out_mat2,out_mat2,out_mat2,out_mat2,
                                out_mat2,out_mat2,out_mat2,out_mat2,out_mat2))
    
    plotting_distributions(test_mat = plot_mat,
                           transparency = 0, vector_colors = "black",
                           MainPlotName = paste0("Suggested link ",
                                                 "(",
                                                 lipid_Nr[i],
                                                 ") ",
                                                 plot_name))
        
  }
  
  names(list_out) <- lipid_names
  names(list_out2) <- lipid_names
  
  dev.off()
  
  if(returnObject == "ClassComparisonInput") {
    
    return(list_out2)
    
  } else if (returnObject == "minDataset") {
    
    return(list_out)
    
  } else if (is.null(returnObject)) {
    
    return(list_out2)
    
  }
  
}
