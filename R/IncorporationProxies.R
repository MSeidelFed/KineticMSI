#' A function to produce files that describe through different proxies the tracer dynamics within kineticMSI datasets
#'
#' This function allows to calculate across MSI pixels various values that reflect different aspects of the tracer dynamics. The function tests if the molecular features are shared across all datasets, if these are not shared, the function produces files with the common features before carrying on with the calculations. This is to prevent errors in the joined steady state pool files that are generated. The function outputs the isotope incorporation proxies as csv files to the same input directories.
#' @param ParentDir directory where the input files are stored. This is also the address that will contain the isotope incorporation proxies as new csv files.
#' @param SteadyStatePoolsDir directory where the output steady state pool files are stored after running the function.
#' @param ColSumNorm defaults to FALSE. Allows to normalize the values of the calculated *de novo* synthesized pools by the sum of all isotopologue pools, i.e., sum across feature of M0(A0) to Mn(An).
#' @keywords KineticMSI Tracer Incorporation Dynamics
#' @export
#' @examples
#' 
#' ...


IncorporationProxies <- function(ParentDir,
                                 SteadyStatePoolsDir = NULL,
                                 ColSumNorm = FALSE) {
  
  dir.create(path = SteadyStatePoolsDir, showWarnings = F)
  
  ## functions
  
  split_ <- function(x) {strsplit(x = x, split = "_")[1]}
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  ## main

  ### Managing input files

    Enrichment_files <- list.files(path = ParentDir, pattern = "Corrected.csv",
                                   all.files = FALSE, full.names = TRUE, recursive = TRUE,
                                   ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
    
    RawData <- list.files(path = ParentDir, pattern = "RawData.csv",
                          all.files = FALSE, full.names = TRUE, recursive = TRUE,
                          ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
    
      
    ProduceNewFiles(reps = Enrichment_files)
      
    ProduceNewFiles(reps = RawData)
      
    Corrected_Enrichment_files <- list.files(path = ParentDir, pattern = "Corrected_SharedFeatures.csv",
                                             all.files = FALSE, full.names = TRUE, recursive = TRUE,
                                             ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
      
    if (length(Corrected_Enrichment_files) > 0) {
        
      Enrichment_files <- Corrected_Enrichment_files
        
      RawData <- list.files(path = ParentDir, pattern = "RawData_SharedFeatures.csv",
                            all.files = FALSE, full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
        
    }
      
    
    MolecularSpecies <- unique(apply(list2df(sapply(X = rownames(read.csv(Enrichment_files[1],
                                                                          header = T, row.names = 1)),
                                                    FUN = split_))[1,],
                                     MARGIN = 2,
                                     FUN = as.character))
    
    
    ### Producing the isotopic proxies
    
    SteadyStatePoolsM1M0 <- matrix(NA, nrow = length(Enrichment_files),
                                   ncol = length(MolecularSpecies))
    
    SteadyStatePoolsM1Mn <- matrix(NA, nrow = length(Enrichment_files),
                                   ncol = length(MolecularSpecies))
    
    treatment_names <- c()
    
    for (i in 1:length(Enrichment_files)) {
      
      Enrichment_file_runner <- Enrichment_files[i]
      
      Raw_file_runner <- RawData[i]
      
      ### runner indexes
      
      File_i <- read.csv(file = Enrichment_file_runner, header = T, row.names = 1)
      
      Raw_File_i <- read.csv(file = Raw_file_runner, header = T, row.names = 1)
      
      MolecularSpecies_i <- unique(apply(list2df(sapply(X = rownames(File_i) , FUN = split_))[1,],
                                         MARGIN = 2,FUN = as.character))
      
      name_i <- strsplit(strsplit(Enrichment_file_runner, split = "/")[[1]][length(strsplit(Enrichment_file_runner, split = "/")[[1]])],
                         split = "\\.")[[1]][1]
      
      out_path_i <- paste0(strsplit(Enrichment_file_runner,
                                    split = "/")[[1]][-length(strsplit(Enrichment_file_runner,
                                                                       split = "/")[[1]])], collapse = "/")
      
      ### generating new files
      
      M0s <- File_i[grep(pattern = "_0", x = rownames(File_i)),]
      
      write.csv(M0s, paste0(out_path_i, "/", name_i, "M0s.csv"))
      
      M1s <- File_i[grep(pattern = "_1", x = rownames(File_i)),]
      
      write.csv(M1s, paste0(out_path_i, "/", name_i, "M1s.csv"))
      
      M1M0ratio <- M1s/M0s
      
      write.csv(M1M0ratio, paste0(out_path_i, "/", name_i, "M1M0ratios.csv"))
      
      M1fraction <- M1s/M0s+M1s
      
      write.csv(M1fraction, paste0(out_path_i, "/", name_i, "M1fraction.csv"))
      
      #### grabing isotopologues per feature from the raw data and building steady state pools either from M0+M1 or from M0+...+Mn or de novo synthesized pools
      
      means_isotopologues_Raw <- c()
      
      means_M1M0 <- c()
      
      denovo <- matrix(NA, nrow = length(MolecularSpecies_i), ncol = dim(File_i)[2])
      
      for (j in 1:length(MolecularSpecies_i)) {
        
        feature_i_Raw <- Raw_File_i[grep(pattern = MolecularSpecies_i[j], x = rownames(File_i)),]
        
        feature_i_corr <- File_i[grep(pattern = MolecularSpecies_i[j], x = rownames(File_i)),]
        
        ## mean steady state pool m0:mn
        
        M0_L <- rbind(feature_i_Raw[1,], feature_i_corr[2:dim(feature_i_corr)[1],])
        
        means_isotopologues_Raw[j] <- mean(colSums(M0_L), na.rm = T)
        
        ## mean steady state pool m0:m1
        
        M1_L <- rbind(feature_i_Raw[1,], feature_i_corr[2,])
        
        means_M1M0[j] <- mean(colSums(M1_L), na.rm = T)
        
        ## Colsum mean steady state pool m1:mn (de novo synthesized) normalized or not
        
        ### files will be reconstructed at each iteration like in before the loop...
        
        if (ColSumNorm == T) {
          
          denovo[j,] <- colSums(feature_i_corr[2:dim(feature_i_corr)[1],]) / colSums(M0_L) #colsums of full isotopic envelopes
          
        } else {
          
          denovo[j,] <- colSums(feature_i_corr[2:dim(feature_i_corr)[1],])
          
        }
        
      }
      
      rownames(denovo) <- MolecularSpecies
      colnames(denovo) <- colnames(File_i)
      
      if (ColSumNorm == T) {
        
        write.csv(denovo, paste0(out_path_i, "/", name_i, "_denovo_synthesized_pool_Norm.csv"))
        
      } else {
        
        write.csv(denovo, paste0(out_path_i, "/", name_i, "_denovo_synthesized_pool.csv"))
        
      }
      
      ## steady state pool files
      
      SteadyStatePoolsM1Mn[i,] <- means_isotopologues_Raw
      
      SteadyStatePoolsM1M0[i,] <- means_M1M0
      
      treatment_names[i] <- name_i
      
    }
    
    colnames(SteadyStatePoolsM1Mn) <- MolecularSpecies
    colnames(SteadyStatePoolsM1M0) <- MolecularSpecies
    rownames(SteadyStatePoolsM1Mn) <- treatment_names
    rownames(SteadyStatePoolsM1M0) <- treatment_names
    
    if (is.null(SteadyStatePoolsDir)) {
      
      write.csv(SteadyStatePoolsM1Mn, file = "SteadyStatePoolsM0Mn.csv")
      
      write.csv(SteadyStatePoolsM1M0, file = "SteadyStatePoolsM0M1.csv")
      
    } else {
      
      write.csv(SteadyStatePoolsM1Mn, file = paste0(SteadyStatePoolsDir,
                                                    "/", "SteadyStatePoolsM0Mn.csv"))
      
      write.csv(SteadyStatePoolsM1M0, file = paste0(SteadyStatePoolsDir,
                                                    "/", "SteadyStatePoolsM0M1.csv"))
    }
  
}
