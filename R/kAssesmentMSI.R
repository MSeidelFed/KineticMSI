#' A function to perform an initial step of replicated data quality assessment on KineticMSI datsets
#'
#' This function allows KineticMSI users to assess the quality in terms of reproducibility and mean distribution from measured molecular features across several independent treatments. The function first subsets all datasets to a common vector of molecular features and a common number of pixels. If molecular features are not the same across datatsets, the function produces new sets with the "Shared Features", which are named with this extension in the original IsoCorrectoR folders. Pixels are previously sorted per sample according to their magnitude, e.g., in the case of enrichment percentage from high to low enrichment. Afterwards a plot reflecting the ratio of means from subset and entire sets is produced in order to evaluate whether the subsetting procedure is skewing the data or whether the change trends will be fully preserved. An alternative to sorting by magnitude is sorting by acquisition time, which allows the users later on to evaluate differences that might be constrained in space or dependent on the instrument performance during a single sample. The function returns distribution plots in a PDF file across treatments that allow interpreting shifts in data that would otherwise remain unnoticed. Finally the function also returns either a list with the minimum datasets for all matrices or the compressed matrices used for mean distribution assessment and later for mean class comparison across samples.
#' @param path parent directory were all enrichment files are contained, digs within recursive folders.
#' @param PatternEnrichment defaults to "MeanEnrichment". Defines a character vector used to grab input csv enrichment files that can be later subset.
#' @param SubSetRepsIntensities defaults to FALSE. Allows to subset the MSI file list found in path.
#' @param CompareSampledSet defaults to TRUE. Compares the subsets and entire sets by plotting the ratios between mean and standard deviations across molecular features to allow evaluation of biased sampling. Optimally the ratios should oscillate around 1 with minimum variance.
#' @param returnObject can be either "ClassComparisonInput" or "minDataset" and returns either a list with input matrices per molecular feature for the following mean class comparison function or the minimum data set subset to the minimum amount of pixels.
#' @param factorVector character vector that needs to define the treatments using the same nomenclature and naming scheme as for the input files (follow the exemplary KineticMSI data for details).
#' @param fun2clust allows to decide whether pixels across files are sorted before subseting using their magnitude, i.e., "Enrichment" or their acquisition time, i.e., "AcqTime". 
#' @param logiTransformation Allows users to perform logit transformation before the class distribution assessment.
#' @param ScalingFactorXaxisDensityPlot defaults to NULL = 1 and defines the boundary of the x-axis in the ggplots that assess the distribution of the data.
#' @keywords MSI Replicates Assessment Tracer Dynamics
#' @export
#' @examples
#' 
#' ...


kAssesmentMSI <- function(path,
                          PatternEnrichment = "MeanEnrichment",
                          SubSetRepsIntensities = FALSE,
                          CompareSampledSet = TRUE,
                          returnObject = c(NULL, "ClassComparisonInput", "minDataset"),
                          factorVector,
                          fun2clust = c("Enrichment", "AcqTime"),
                          logiTransformation = F,
                          ScalingFactorXaxisDensityPlot = NULL) {
  
  
  ## needed functions
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  ## Main
  
  ### Importing all files
  
  reps <- list.files(path = path, pattern = PatternEnrichment,
                     all.files = FALSE, full.names = TRUE, recursive = TRUE,
                     ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  
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
    ##### full matrices
    
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
    
    #### comparison plot
    
    pdf(file = "ComparedSampledSet.pdf")
    
    par(mfrow = c(2,1))
    
    plot(CompleteSet_Means/MinSet_Means, main = "Means")
    plot(CompleteSet_SD/MinSet_SD, main = "SD")
    
    dev.off()
    
  }
  
  ### Assessing kMSI data
  
  pdf(file = "Data_assesment.pdf")
  
  par(mfrow = c(1,1))
  
  list_out <- list()
  
  list_out2 <- list()
  
  sig_clust_Nr <- c()
  
  lipid_names <- c()
  
  for (i in 1:length(feature_Nr)) {
    
    cat("...\n")
    cat(paste0("Processing: ", feature_Nr[i]," - ", i, "\n"))
    cat("...\n")
    
    if (fun2clust == "Enrichment") {
      
      runner <- EnrichmentSortedSamplingMSI(path = reps,
                                            RepsOrPath = "reps", 
                                            PatternEnrichment = PatternEnrichment,
                                            featureChoice = i)
      
    } else if (fun2clust == "AcqTime") {
      
      runner <- AcqTimeRandomSamplingMSI(path = reps,
                                         RepsOrPath = "reps",
                                         PatternEnrichment = PatternEnrichment,
                                         featureChoice = i)
      
    } else {
      
      cat("Error: select an appropriate sorting function")
      
    }
    
    colnames(runner[[1]]) <- paste0(rep("X", dim(runner[[1]])[2]), c(1:dim(runner[[1]])[2]))
    
    runner_WO_zeros <- as.matrix(runner[[1]][,as.numeric(colSds(runner[[1]])) != 0])
    
    #### transforming to logit if set to TRUE
    
    if (logiTransformation == T) {
      
      logit_transformation <- function(x){log(x /(1 - x))}
      
      ##### Replace zero with a very small value so it can fit in a model
      
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
    
    #### generating the first output (minimum data set)
    
    list_out[[i]] <- runner_WO_zeros
    lipid_names[i] <- feature_Nr[i]
    
    #### plotting distributions
    
    split_file_names <- list2df(strsplit(rownames(runner_WO_zeros), "/"))
    
    
    full_file_names <- list2df(strsplit(apply(split_file_names[dim(split_file_names)[1],], 2, as.character), "_"))
    
    rep_name <- paste0(apply(full_file_names[2,], 2, as.character),
                       "_",
                       apply(full_file_names[3,], 2, as.character))
    
    if (dim(runner_WO_zeros)[2] < 1) {
      
      NULL
      
    } else {
      
      ClassDistribution(inMat = t(runner_WO_zeros), Treatments = as.factor(rep_name),
                        plots = T, returnTable = F, factorVector = factorVector,
                        PlotMain = feature_Nr[i], ScalingFactorXaxis = ScalingFactorXaxisDensityPlot)
      
    }
    
    #### generating the second output (class comparison input & distribution test)
    
    out_mat2 <- as.matrix(rowMeans(runner_WO_zeros))
    
    list_out2[[i]] <- out_mat2
    
    if (dim(runner_WO_zeros)[2] < 1) {
      
      NULL
      
    } else {
      
      normality_test <- paste0(shapiro.test(out_mat2)[["method"]],
                               ": p = ",
                               round(shapiro.test(out_mat2)[["p.value"]], 2))
      
      plot_name <- unique(testing_distributions(Distribution_test_mat = cbind(as.numeric(out_mat2[,1]),
                                         as.numeric(out_mat2[,1]))))
      
      plot_mat <- as.matrix(cbind(out_mat2,out_mat2,out_mat2,out_mat2,out_mat2,
                                  out_mat2,out_mat2,out_mat2,out_mat2,out_mat2,
                                  out_mat2,out_mat2,out_mat2,out_mat2,out_mat2,
                                  out_mat2,out_mat2,out_mat2,out_mat2,out_mat2))
      
      plotting_distributions(test_mat = plot_mat,
                             transparency = 0, vector_colors = "black",
                             MainPlotName = c(paste0("Suggested link ",
                                                     "(",
                                                     feature_Nr[i],
                                                     ") ",
                                                     plot_name),
                                              normality_test))
      
    }
    
  }
  
  names(list_out) <- lipid_names
  names(list_out2) <- lipid_names
  
  dev.off()
  
  ### returning output
  
  if (returnObject == "minDataset") {
    
    return(list_out)
    
  } else if(is.null(returnObject) | returnObject == "ClassComparisonInput") {
    
    return(list_out2)
    
  }
  
}
