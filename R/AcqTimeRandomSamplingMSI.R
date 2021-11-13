#' A function to sort previous to unification independent kMSI samples according to pixel acquisition time
#'
#' This function allows KineticMSI functions to build a matrix with independent samples by unifying them.
#' @param path parent directory were all enrichment files are contained, digs within recursive folders.
#' @param PatternEnrichment defaults to "MeanEnrichment". Inherits from parent function.
#' @param RepsOrPath defines whether the input is a vector with paths to individual reps or already a path itself.
#' @param featureChoice Allows selecting the index of a feature of interest to be sorted.
#' @keywords MSI Subsetting Tracer Dynamics
#' @export
#' @examples
#' 
#' ...


AcqTimeRandomSamplingMSI <- function(path,
                                      PatternEnrichment = "MeanEnrichment.csv",
                                      RepsOrPath = c("reps", "path"),
                                      featureChoice = 0) {
  
  ### list2df
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  } 
  
  #### choosing a feature
  
  if (RepsOrPath == "reps") {
    
    reps <- path
    
  } else if (RepsOrPath == "path"){
    
    reps <- list.files(path = path, pattern = PatternEnrichment, all.files = FALSE,
                       full.names = T, recursive = T,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
  }
  
  #### optional argument that activates itself if a feature index is not defined in the call
  
  if (featureChoice == 0) {
    
    file_test <- read.csv(file = reps[1], header = T, row.names = 1)
    
    feature_Nr <- unique(rownames(file_test))
    
    print(feature_Nr)
    
    featureChoice <- readline(prompt="Enter index of feature of interest: ")
    
  }
  
  #### defining the minimum set length
  
  out_vector <- c()
  
  for (i in 1:length(reps)) {
    
    file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
    
    out_vector[i] <- dim(file_runner)[2]
    
  }
  
  min_dataset <- min(out_vector)
  
  #### sampling randomly but conserving acquisition time relationships among pixels
  
  out_single_feature <- matrix(NA, length(reps), min_dataset)
  
  out_single_feature_coords <- matrix(NA, length(reps), min_dataset)
  
  for (i in 1:length(reps)) {
    
    file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
    
    feature_runner <- as.matrix(sample(file_runner[as.numeric(featureChoice),],
                                     min_dataset,
                                     replace = F))
    
    df2sort <- as.data.frame(t(rbind(apply(list2df(strsplit(colnames(feature_runner),
                                                            "\\."))[2,],
                                           2,
                                           as.numeric),
                                     feature_runner)))
    
    colnames(df2sort) <- c("pixelNr", "Enr")
    
    SortedEnrichment <- df2sort[order(df2sort$pixelNr),2]
    
    SortedCoords <- paste0("Spot.", df2sort[order(df2sort$pixelNr),1])
    
    out_single_feature[i,] <- SortedEnrichment
    
    out_single_feature_coords[i,] <- SortedCoords
    
  }
  
  rownames(out_single_feature) <- reps
  
  return(list(out_single_feature,
              out_single_feature_coords))
  
}
