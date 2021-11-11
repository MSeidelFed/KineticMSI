#' A function to produce new input files at multiple kineticMSI steps if the molecular species differ among the provided input files
#'
#' This function allows KineticMSI functions to produce new input files given uneven number of features.
#' @param reps must be a list of replicate files.
#' @keywords KineticMSI Preprocessing troubleshooting
#' @export
#' @examples
#' 
#' ...

ProduceNewFiles <- function(reps) {
  
  ## functions needed

  ReplaceOldFiles <- function(reps,
                              featureNr) {
    
    reps_new <- c() 
    
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
      
      reps_new[i] <- paste0(out_dir, "/", file_name)
      
      write.csv(file_test[featureNr,], file = reps_new[i],
                quote = F)
      
    }
    
    return(reps_new)
    
  }
  
  ## Main

  ### getting the common feature_Nr for all datasets
  
  feature_Nr <- list()
  
  feature_Nr_vec <- c()
  
  for (i in 1:length(reps)) {
    
    file_test <- read.csv(file = reps[i], header = T, row.names = 1)
    
    feature_Nr[[i]] <- unique(rownames(file_test))
    
    feature_Nr_vec <- c(feature_Nr_vec, rownames(file_test))
    
  }
  
  feature_Nr = Reduce(intersect, feature_Nr)
  
  feature_Nr_vec = unique(feature_Nr_vec)
  
  ### if common features differ either in No. or order then produce new files, otherwise nothing happens
  
  if (length(feature_Nr) != length(feature_Nr_vec)) {
    
    reps = ReplaceOldFiles(reps = reps, featureNr = feature_Nr)
    
  } else if (length(which((feature_Nr == feature_Nr_vec) == F)) > 0) {
    
    reps = ReplaceOldFiles(reps = reps, featureNr = feature_Nr)
    
  } 
  
}
