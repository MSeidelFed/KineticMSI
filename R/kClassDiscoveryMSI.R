#' A function to subset consolidated KineticMSI datasets into groups of related pixels
#'
#' This function allows users to subset consolidated KineticMSI matrices (see our kAssesmentMSI.R function for details) into subsets that are validated internally by the data structure. The function clusters consolidated data matrices using a hierarchical clustering algorithm (HCA) with bootstrapping, which is a dependency from the R package pvclust. Then using a user defined significance threshold the function grabs an optimized number of significant clusters that feature AU-P values above the threshold. The optimization increases the threshold when many significant clusters are found until it optimizes the outcome to the minimum cluster number, or if there are any significant clusters then the optimization parameter lowers the threshold until obtaining significant partitions. The outcome from this tuning process is returned to the console. The function returns to the R environment a list of abundances and coordinates matrices per molecular feature that can be used to both map the clusters onto the original MSI images using kReconstructMSI.R or perform class comparison using the clustered subsets using kSubSetClassComparisonMSI.R. Additionally the function returns to the working directory a number of diagnostic plots that allow users to better contextualize the partitions in their datasets.
#' @param path parent directory were all enrichment files are contained, digs within recursive folders.
#' @param PatternEnrichment defaults to "MeanEnrichment". Defines a character vector used to grab input csv enrichment files that can be later subset.
#' @param DistMethod argument inherited from method.dist from pvclust.
#' @param nboot argument inherited from pvclust.
#' @param alpha defaults to 0.9, allows users to define the threshold AU-P value above which clusters are selected as significantly supported by the data. For details see pvclust package.
#' @param SigClustHist defaults to TRUE. Allows to return a histogram featuring the densities of significant cluster frequencies across molecular features.
#' @param SubSetRepsIntensities defaults to FALSE. Allows to subset the MSI file list found in path.
#' @param CompareSampledSet defaults to TRUE. Compares the subsets and entire sets by plotting the ratios between mean and standard deviations across molecular features to allow evaluation of biased sampling. Optimally the ratios should oscillate around 1 with minimum variance.
#' @param returnObject can be either "RowMeansDataset" or "minDatasetPlusCoords" and returns either a list with input matrices per molecular feature for the following mean class comparison (only the rowmeans per feature) or the minimum data subset to the minimum amount of pixels per entity.
#' @param fun2clust allows to decide whether pixels across files are sorted before subsetting using their magnitude, i.e., "Enrichment" or their acquisition time, i.e., "AcqTime". 
#' @param ZeroAction allows to define what happens with zeros within the entity data matrices. If null, the zeros are fully preserved, if "remove" the rows that only contain zeros are deleted and if "replace" all zeros are replaced by normally distributed randomly generated numbers in the scale of 10^-12 to 10^-13.
#' @param logiTransformation allows users to perform logit transformation before the class distribution assessment.
#' @keywords MSI Replicates Subsetting Significant Clustering Tracer Dynamics
#' @export
#' @examples
#' 
#' ...


kClassDiscoveryMSI <- function(path,
                               PatternEnrichment = "MeanEnrichment",
                               DistMethod = "abscor",
                               nboot = 100,
                               alpha = 0.9,
                               SigClustHist = TRUE,
                               SubSetRepsIntensities = FALSE,
                               CompareSampledSet = TRUE,
                               returnObject = c(NULL, "RowMeansDataset", "minDatasetPlusCoords"),
                               fun2clust = c("Enrichment", "AcqTime"),
                               ZeroAction = c(NULL, "remove", "replace"),
                               logiTransformation = FALSE) {
  
  
  
  
  ## functions needed for the main
  
  ### numeric.factor
  
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  
  ### handling errors
  
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
  
  
  
  ## Main
  
  original_alpha = alpha
  
  ### grabbing data
  
  reps <- list.files(path = path, pattern = PatternEnrichment, all.files = FALSE,
                     full.names = T, recursive = T,
                     ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  ### Sub-setting files
  
  if (SubSetRepsIntensities == T) {
    
    reps = SubSetInputFiles(reps = reps)
    
  }
  
  ### Producing files with shared features if necessary
  
  reps1 = ProduceNewFiles(reps = reps)
  
  if(!is.null(reps1)) {
    
    reps = reps1
    
  }
  
  ### getting the common features for all datasets to perform analyses
  
  feature_Nr <- list()
  
  feature_Nr_vec <- c()
  
  for (i in 1:length(reps)) {
    
    file_test <- read.csv(file = reps[i], header = T, row.names = 1)
    
    feature_Nr[[i]] <- unique(rownames(file_test))
    
    feature_Nr_vec <- c(feature_Nr_vec, rownames(file_test))
    
  }
  
  feature_Nr = Reduce(intersect, feature_Nr)
  
  ### comparing sampled and full datasets if set to TRUE
  
  if (CompareSampledSet == T) {
    
    #### getting rowmeans of all files to compare to the sampled ones
    
    #### full matrices
    
    CompleteSet_Means <- c()
    CompleteSet_SD <- c()
    
    for (i in 1:length(reps)) {
      
      runner_file <- read.csv(reps[i], header = T, row.names = 1)
      colnames(runner_file) <- NULL
      
      CompleteSet_Means <- c(CompleteSet_Means, rowMeans(runner_file))
      CompleteSet_SD <- c(CompleteSet_SD, rowSds(as.matrix(runner_file)))
      
    }
    
    #### sampled matrices
    
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
  
  ### finding a significant number of partitions in joined datasets across metabolites and samples
  
  pdf(file = "Bootstrapped_Dendrograms.pdf")
  
  par(mfrow = c(1,1))
  
  list_out <- vector(mode = "list", length = length(feature_Nr))
  
  list_out2 <- vector(mode = "list", length = length(feature_Nr))
  
  sig_clust_Nr <- c()
  
  feature_names <- c()
  
  for (i in 1:length(feature_Nr)) {
    
    #### restarting alpha for each feature
    
    alpha = original_alpha
    
    #### clustering pixels after a defined sorting...
    
    if (fun2clust == "Enrichment") {
      
      runner_sampling <- EnrichmentSortedSamplingMSI(path = reps,
                                                     RepsOrPath = "reps", 
                                                     PatternEnrichment = PatternEnrichment,
                                                     featureChoice = i)
      
      runner_coords <- runner_sampling[[2]]
      
      runner <- runner_sampling[[1]]
      
    } else if (fun2clust == "AcqTime") {
      
      runner_sampling <- AcqTimeRandomSamplingMSI(path = reps,
                                                  RepsOrPath = "reps",
                                                  PatternEnrichment = PatternEnrichment,
                                                  featureChoice = i)
      
      runner_coords <- runner_sampling[[2]]
      
      runner <- runner_sampling[[1]]
      
    } else {
      
      stop("Error: select an appropriate sorting function")
      
    }
    
    colnames(runner) <- paste0(rep("X", dim(runner)[2]), c(1:dim(runner)[2]))
    
    colnames(runner_coords) <- paste0(rep("X", dim(runner)[2]), c(1:dim(runner)[2]))
    
    #### Zero action taken on the datasets
    
    if(is.null(ZeroAction)) {
      
      if(length(which(colSds(runner) == 0)) > 0) {
        
        stop("ERROR: there are columns in the assembled matrices with NULL sd, please select an appropriate ZeroAction parameter in the function call")
        
      }
      
      runner_WO_zeros <- as.matrix(runner)
      
      runner_coords_WO_zeros <- as.matrix(runner_coords)
      
      rownames(runner_coords_WO_zeros) <- rownames(runner_WO_zeros)
      
    } else if(ZeroAction == "remove") {
      
      runner_WO_zeros <- as.matrix(runner[,as.numeric(colSds(runner)) != 0])
      
      runner_coords_WO_zeros <- as.matrix(runner_coords[,as.numeric(colSds(runner)) != 0])
      
      rownames(runner_coords_WO_zeros) <- rownames(runner_WO_zeros)
      
    } else if(ZeroAction == "replace") {
        
      runner[runner == 0] <- abs(rnorm(length(which(runner == 0)),
                                       mean = 0.0000000000001,
                                       sd = 0.00000000001))
      
      runner_WO_zeros <- as.matrix(runner)
      
      runner_coords_WO_zeros <- as.matrix(runner_coords)
      
      rownames(runner_coords_WO_zeros) <- rownames(runner_WO_zeros)
      
    }
    
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
                              quiet = T, method.dist = DistMethod)[1]) == "error") {
      
      #### getting rid of metabolites that produce an error in the pvclust function
      
      #### need to produce NULLs here to coincide with the feature Nr.
      
      feature_names[i] <- feature_Nr[i]
      
      cat("...\n")
      cat(paste0("Error: ", feature_Nr[i], " - ", i, "\n"))
      cat("...\n")
      
    } else {
      
      cat(paste0("Successful: ", feature_Nr[i]," - ", i, "\n"))
      
      ## bootstrapping the clustering on the final object
      
      HCA_boot_feature <- pvclust(runner_WO_zeros,
                                  method.hclust = "average",
                                  method.dist = DistMethod,
                                  nboot = nboot, 
                                  quiet = T)
      
      if (dim(runner_WO_zeros)[2] > 0.2 * dim(runner)[2]) { 
        
        #### plot only works on matrices with 20% or more columns with values different than zero
        
        table_features <- pvpick(HCA_boot_feature, alpha = alpha, pv="au")
        
        #### changing the alpha in dendrograms with one significant cluster
        
        count = 0
        
        alpha_out <- alpha
        
        Initial_k = length(table_features$edges)
        
        if (length(table_features$edges) == 1 | length(table_features$edges) > 4) {
          
          cat("k before relaxation: ", length(table_features$edges), "\n")
          
        }
        
        ## lower limit until AU = .77
        
        while (length(table_features$edges) == 1 | length(table_features$edges) > 4) {
          
          count = count + 1
          
          if (count < 7) {
            
            alpha = alpha-((1-alpha)/8)
            
            table_features <- pvpick(HCA_boot_feature, alpha = alpha, pv="au")
            
            alpha_out = alpha
            
          } else {
            
            break
            
          }
          
        }
        
        ## upper limit until AU = .99
        
        count = 0
        
        while (length(table_features$edges) == 1 | length(table_features$edges) > 4) {
          
          count = count + 1
          
          if (count < 30) {
            
            alpha = alpha+((1-alpha)/8)
            
            table_features <- pvpick(HCA_boot_feature, alpha = alpha, pv="au")
            
            alpha_out = alpha
            
          } else {
            
            break
            
          }
          
        }
        
        if (Initial_k < length(table_features$edges)) {
          
          alpha_out = original_alpha
          
          table_features <- pvpick(HCA_boot_feature, alpha = alpha_out, pv="au")
          
        }
        
        cat("Final k : ", length(table_features$edges), "\n")
        
        cat("Final AU-P: ", alpha_out, "\n")
        
        #### plotting dendrograms with significant clusters
        
        plot(x = HCA_boot_feature, print.pv = "au", print.num = T, cex = 1, lwd = 2, main = feature_Nr[i])
        
        pvrect(HCA_boot_feature, alpha = alpha, pv="au")
        
        #### building matrices using significant clusters as factors
        
        if (length(table_features$edges) < 10 & length(table_features$edges) > 0) {
          
          ##### Metabolites with more than 10 parental partitions are neglected
          
          sig_clust_Nr[i] <- length(table_features$edges)
          
          feature_names[i] <- feature_Nr[i]
          
          #### mining out the cluster to build means
          
          ##### establishing the maximum cluster size and the pixel proportions
          
          out_vector <- c()
          
          for (k in 1:length(table_features$clusters)) {
            
            out_vector[k] <- length(table_features$clusters[[k]])
            
          }
          
          max_cluster_length <- max(out_vector)
          
          ##### grabbing clusters with colnames of clustered pixels into a data frame
          
          clusters_p95 <- sapply(1:length(as.matrix(table_features$clusters)),function(j) as.matrix(table_features$clusters)[[j]][1:max_cluster_length])
          
          #### subsetting the runner_WO_zeros object and getting the means
          
          ##### producing return objects
          
          out_mat <- matrix(NA, nrow = dim(runner_WO_zeros)[1], ncol = sig_clust_Nr[i])
          
          out_list_coords <- list()
          
          out_list <- list()
          
          for (n in 1:dim(clusters_p95)[2]) {
            
            out_mat[,n] <- rowMeans(runner_WO_zeros[,na.omit(clusters_p95[,n])])
            
            out_list_coords[[n]] <- runner_coords_WO_zeros[,na.omit(clusters_p95[,n])]
            
            out_list[[n]] <- runner_WO_zeros[,na.omit(clusters_p95[,n])]
            
          }
          
          colnames(out_mat) <- table_features$edges
          
          names(out_list) <- table_features$edges
          
          names(out_list_coords) <- table_features$edges
          
          result_list <- list(out_list, out_list_coords)
          
          names(result_list) <- c("Values", "Coordinates")
          
          ## return row means per cluster
          
          list_out[[i]] <- out_mat
          
          ## return minDataSet; values plus coords
          
          list_out2[[i]]  <- result_list
          
          #### plot distributions with relevant stats
          
          total_pixel_No_clusters = sum(out_vector)
          
          cluster_means = colMeans(out_mat)
          
          out_pixel_proportion <- c()
          
          for (k in 1:length(table_features$clusters)) {
            
            out_pixel_proportion[k] <- paste0("Pixel Percentage: ",
                                              round((out_vector[k]/total_pixel_No_clusters)*100,
                                                    digits = 2),
                                              "% - Mean: ",
                                              round(cluster_means[k], 3))
            
          }
          
          plot(density(runner_WO_zeros), main = feature_Nr[i],
               xlab = "Enrichment Proportion",
               ylim = c(0, max(density(runner_WO_zeros)[["y"]])+5))
          
          legend("topleft", legend=c(paste0("k: ", length(table_features$edges)),
                                     out_pixel_proportion)) # col=c("red", "blue"), lty=1:2, cex=0.8
          
          cat("...\n")
          
        } else {
          
          feature_names[i] <- feature_Nr[i]
          
          cat(paste0("ommited - too few or too many significant K classes ", feature_Nr[i]," - ", i, "\n"))
          cat("...\n")
          
        }
        
      } else {
        
        feature_names[i] <- feature_Nr[i]
        
        cat(paste0("Error in dendrogram plot (too many Zeros): ", feature_Nr[i]," - ", i, "\n"))
        cat("...\n") 
        
      }
    }
  }
  
  names(list_out) <- feature_names
  names(list_out2) <- feature_names
  
  dev.off()
  
  ### heavy-tailed distribution of cluster number
  
  if (SigClustHist == T) {
    
    png(filename = "SigClustHist.png")
    plot(density(na.omit(as.numeric(sig_clust_Nr))),
         xlab = "Significant Cluster No.",
         main = "Cluster No. across features")
    dev.off()
    
  }
  
  if(is.null(returnObject) | returnObject == "RowMeansDataset") {
    
    return(list_out)
    
  } else if (returnObject == "minDatasetPlusCoords") {
    
    return(list_out2)
    
  } 
  
}

