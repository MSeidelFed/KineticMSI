#' A function to sort previous to unification independent kMSI samples according to enrichment percentage or any other kind of proxy for tracer incorporation
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


EnrichmentSortedSamplingMSI <- function(path,
                                        PatternEnrichment = "MeanEnrichment.csv",
                                        RepsOrPath = c("reps", "path"),
                                        featureChoice = 0) {

  
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
  
  #### sampling randomly but sorting by enrichment percentage
  
  out_single_feature <- matrix(NA, length(reps), min_dataset)
  
  out_single_feature_coords <- matrix(NA, length(reps), min_dataset)
  
  for (i in 1:length(reps)) {
    
    file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
    
    Sampled_set <- apply(sample(file_runner[as.numeric(featureChoice),],
                                min_dataset,
                                replace = F), 2, as.numeric)
    
    Sampled_set[which(is.na(Sampled_set))] <- 0
    
    Sorted_set <- sort(Sampled_set, T)
    
    out_single_feature[i,] <- Sorted_set
    
    out_single_feature_coords[i,] <- names(Sorted_set)
    
  }
  
  rownames(out_single_feature) <- reps
  
  return(list(out_single_feature,
              out_single_feature_coords))
  
}

