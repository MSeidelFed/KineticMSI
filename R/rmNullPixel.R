#' A function to crop input kMSI datasets
#'
#' This function allows you to remove MSI pixels that would impair interpretation of true 0% enrichment in the downstream calculations. The function generates corrected csv files and a list with the corrected matrices as a return object in the R environment.The function takes an entire directory and it grabs all csv files within the provided directory. The function grabs each isotopologue envelope and sets to NA all of those pixels that would produce a misinterpretation of the NIA correction leading to misinterpreted enrichment percentages.
#' @param MeasurementFileDir directory where the input files are stored.
#' @param pattern defaults to "csv". It is used to define the pattern on which input files are looked for in the MeasurementFileDir
#' @param SubSetReps defaults to FALSE. Allows to subset the file list found in MeasurementFileDir.
#' @param csvReturn defaults to TRUE. Returns corrected csv files.
#' @param OnlyDeletePixelsWOIsotopologs defaults to FALSE. When TRUE, it allows to only remove MSI pixels that lack isotopologues, while preserving those pixels that only have isotopologues, which  may be an indication of full stable isotope incorporation.
#' @param verbose defaults to FALSE. When TRUE it returns to the console the progression across the input files. Thus the parameter is meant to allow users to spot errors in the input files. 
#' @param verboseFeature defaults to FALSE. When TRUE it returns to the console the progression across molecular variables. Thus the parameter is meant to allow users to spot errors in the input files.
#' @param rmDataStore needs to be defined. Specifies the output directory. Either a new directory "NewDir" or the same directory as the input "InputDir".
#' @param outdir only useful when rmDataStore = "NewDir". Defines the name of the new directory in which output is stored.
#' @keywords KineticMSI Preprocessing
#' @export
#' @examples
#' 
#' ...


rmNullPixel <- function(MeasurementFileDir,
                        pattern = "csv",
                        SubSetReps = F,
                        csvReturn = T,
                        OnlyDeletePixelsWOIsotopologs = F,
                        verbose = F,
                        verboseFeature = F,
                        rmDataStore = c("NewDir", "InputDir"),
                        outdir = "rmOutput"){
  
  
  ### functions needed
  
  #### list2df
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  #### Null Pixel Deplete
  
  NullPixel_deplet <- function(MeasurementFileDir,
                               verbose = F,
                               OnlyDeleteIsotopologPeak = OnlyDeleteIsotopologs) {
    
    MeasurementFile <- read.csv(file = MeasurementFileDir, header = T, row.names = 1)
    
    lipids <- unique(apply(X = list2df(strsplit(x = rownames(MeasurementFile),
                                                split = "_"))[1,],
                           MARGIN = 2,
                           FUN = as.character))
    
    list_isotopologs <- list()
    
    count = 0
    
    for (i in 1:length(lipids)) {
      
      count = count + 1
      
      if (verbose == T) {
        
        cat("...", "\n")
        cat("Processing Feature No. :", i, "i.e.,", lipids[i])
        cat("...", "\n")
        
      }
      
      isotopolog_envelope <- MeasurementFile[grep(pattern = lipids[i], x = rownames(MeasurementFile)),]
      
      count_NA_0 = 0
      
      for (j in 1:dim(isotopolog_envelope)[2]) {
        
        ## Conserves pixel intact if there are isotopologes and user is not interested in depleting pixels without monoisotopic
        
        if (length(which(is.na(isotopolog_envelope))) > 0 &
            count_NA_0 == 0) {
          
          stop("ERROR: your input files contain NAs; NAs are not allowed for this function")
          
        } else if (isotopolog_envelope[1,j] == 0 &
            OnlyDeleteIsotopologPeak == T &
            sum(isotopolog_envelope[2:dim(isotopolog_envelope)[1],j]) != 0) {
          
          isotopolog_envelope[,j] = isotopolog_envelope[,j]
          
          ## Sets to NA pixel if user is interested in depleting pixels without monoisotopic
          
        } else if (isotopolog_envelope[1,j] == 0 &
                   OnlyDeleteIsotopologPeak == F &
                   sum(isotopolog_envelope[2:dim(isotopolog_envelope)[1],j]) != 0) {
          
          isotopolog_envelope[,j] = NA
          
        ## Sets to NA pixels without monoisotopic nor isotopologs to prevent false enrichment interpretation
          
        } else if (isotopolog_envelope[1,j] == 0 &
                   sum(isotopolog_envelope[2:dim(isotopolog_envelope)[1],j]) == 0) {
          
          isotopolog_envelope[,j] = NA
          
        ## Sets to NA pixels without isotopologs to prevent false enrichment interpretation
          
        } else if (sum(isotopolog_envelope[2:dim(isotopolog_envelope)[1],j]) == 0) {
          
          isotopolog_envelope[,j] = NA
          
        } else { runner = 0 
        
        }
        
        count_NA_0 = count_NA_0 + 1
      
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
  
  
  ### Main
  
  reps_test <- list.files(path = MeasurementFileDir, pattern = pattern,
                          all.files = FALSE, full.names = TRUE, recursive = TRUE,
                          ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  
  if (length(reps_test) > 0) {
    
    
    reps <- list.files(path = MeasurementFileDir, pattern = pattern,
                       all.files = FALSE, full.names = TRUE, recursive = TRUE,
                       ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
    
    if (SubSetReps == T) {
      
      reps = SubSetInputFiles(reps = reps)
      
    }
    
    count = 0
    
    zero_rm_MatList <- list()
    
    for (i in 1:length(reps)) {
      
      count = count + 1
      
      if (verbose == T) {
        
        cat("...", "\n")
        cat("Processing File No. :", i, "i.e.,", reps[i])
        cat("...", "\n")
        
      }
      
      zero_rm_MatList[[i]] <- NullPixel_deplet(MeasurementFileDir = reps[i],
                                               verbose = verboseFeature, 
                                               OnlyDeleteIsotopologPeak = OnlyDeletePixelsWOIsotopologs)
      
    }
    
  } else {
    
    stop("ERROR: no files found under the given pattern parameter in the function call")
    
  }
  
  
  if (csvReturn == T) {
    
    for (i in 1:length(zero_rm_MatList)) {
      
      test_print <- cbind(rownames(zero_rm_MatList[[i]]), zero_rm_MatList[[i]])
      colnames(test_print) <- c("Measurements/Samples", colnames(zero_rm_MatList[[i]]))
      
      file_name_spl <- strsplit(reps[i], split = "/")[[1]] 
      
      rep_name <- strsplit(file_name_spl[[length(strsplit(reps[i], split = "/")[[1]])]], "\\.")[[1]][1]
      
      if (rmDataStore == "NewDir") {
        
        dir.create(path = outdir, showWarnings = F)
        
        file_name <- paste0(outdir, "/", rep_name,  "_rm0.csv")
        
        
      } else if (rmDataStore == "InputDir"){
        
        file_name <- paste0(paste0(file_name_spl[-length(file_name_spl)], collapse = "/"), "/", rep_name,  "_rm0.csv")
        
      }
      
      write.csv(x = test_print, file = file_name, row.names = F, quote = F)
      
    }
    
  } else {
    
  }
  
  return(zero_rm_MatList)
  
}

