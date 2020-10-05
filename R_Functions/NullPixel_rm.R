





NullPixel_rm <- function(MeasurementFile_dir,
                         return_csv){
  
  ### functions needed
  
  NullPixel_deplet <- function(MeasurementFile_dir) {
    
    list2df <- function(x) 
    { 
      MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
      DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
      colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
      DF 
    }
    
    MeasurementFile <- read.csv(file = MeasurementFile_dir, header = T, row.names = 1)
    
    lipids <- unique(apply(X = list2df(strsplit(x = rownames(MeasurementFile),
                                                split = "_"))[1,],
                           MARGIN = 2,
                           FUN = as.character))
    
    list_isotopologs <- list()
    
    for (i in 1:length(lipids)) {
      
      isotopolog_envelope <- MeasurementFile[grep(pattern = lipids[i],
                                                  x = rownames(MeasurementFile)),]
      empty_pixels <- NULL
      
      for (j in 1:dim(isotopolog_envelope)[2]) {
        
        if (isotopolog_envelope[1,j] == 0) {
          
          isotopolog_envelope[,j] = NA
          
        } else if (sum(isotopolog_envelope[2:dim(isotopolog_envelope)[1],j]) == 0) {
          
          isotopolog_envelope[,j] = NA
          
        } else { runner = 0 }
        
        isotopolog_envelope[,j] = isotopolog_envelope[,j]
      }
      list_isotopologs[[i]] <- isotopolog_envelope
    }
    
    ### assembling a matrix from the different length isotopologs
    
    out_df <- matrix(NA, nrow = 0, ncol = dim(MeasurementFile)[2])
    
    for (i in 1:length(list_isotopologs)) {
      
      runner <- list_isotopologs[[i]]
      
      out_df <-  rbind(out_df, runner)
      
    }
    return(out_df)
  }
  
  ## Main
  
  reps_test <- list.files(path = MeasurementFile_dir, pattern = ".csv", all.files = FALSE,
                          full.names = T, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  zero_rm_MatList <- list()
  
  if (length(reps_test) > 0) {
    
    
    reps <- list.files(path = MeasurementFile_dir, pattern = ".csv", all.files = FALSE,
                       full.names = T, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
    for (i in 1:length(reps)) {
      
      zero_rm_MatList[[i]] <- NullPixel_deplet(reps[i])
      
    }
    
  } else {
    
    zero_rm_MatList[[1]] <- NullPixel_deplet(MeasurementFile_dir = MeasurementFile_dir)
    
  }
  
  
  
  if (return_csv == T) {
    
    for (i in 1:length(zero_rm_MatList)) {
      
      test_print <- cbind(rownames(zero_rm_MatList[[i]]), zero_rm_MatList[[i]])
      colnames(test_print) <- c("Measurements/Samples", colnames(zero_rm_MatList[[i]]))
      
      file_name_spl <- strsplit(reps[i], split = "\\.")[[1]]
      
      file_name <- paste0(file_name_spl[1], "_rm0.", file_name_spl[2])
      
      write.csv(x = test_print, file = file_name, row.names = F, quote = F)
      
    }
    
  } else {
    
  }
  
  return(zero_rm_MatList)
  
}

