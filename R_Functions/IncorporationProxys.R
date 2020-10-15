




IncorporationProxys <- function(Parent_dir, SteadyStatePools_dir = NULL) {
  
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
  
  if (length(strsplit(Parent_dir,split = "/")[[1]]) == 2) {
    
    Enrichment_files <- list.files(path = Parent_dir, pattern = "Corrected.csv",
                                   full.names = T, recursive = T)
    
    RawData <- list.files(path = Parent_dir, pattern = "RawData.csv",
                          full.names = T, recursive = T)
    
    MolecularSpecies <- unique(apply(list2df(sapply(X = rownames(read.csv(Enrichment_files[1],
                                                                          header = T, row.names = 1)),
                                                    FUN = split_))[1,],
                                     MARGIN = 2,
                                     FUN = as.character))
    
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
      
      name_i <- strsplit(strsplit(Enrichment_file_runner, split = "/")[[1]][5],
                         split = "\\.")[[1]][1]
      
      treatment_i <- strsplit(strsplit(strsplit(Enrichment_file_runner, split = "/")[[1]][5],
                                       split = "\\.")[[1]][1], split = "_")[[1]][2]
      
      rep_i <- strsplit(strsplit(strsplit(Enrichment_file_runner, split = "/")[[1]][5],
                                 split = "\\.")[[1]][1], split = "_")[[1]][3]
      
      out_path_i <- paste0(Parent_dir, "/", strsplit(Enrichment_file_runner,
                                                     split = "/")[[1]][4], "/")
      
      ### generating new files
      
      M0s <- File_i[grep(pattern = "_0", x = rownames(File_i)),]
      
      write.csv(M0s, paste0(out_path_i, name_i, "M0s.csv"))
      
      M1s <- File_i[grep(pattern = "_1", x = rownames(File_i)),]
      
      write.csv(M1s, paste0(out_path_i, name_i, "M1s.csv"))
      
      M1M0ratio <- M1s/M0s
      
      write.csv(M1M0ratio, paste0(out_path_i, name_i, "M1M0ratios.csv"))
      
      M1fraction <- M1s/M0s+M1s
      
      write.csv(M1fraction, paste0(out_path_i, name_i, "M1fraction.csv"))
      
      ## need to grab isotopologues per lipid from the raw data and then build steady state pools 
      ## either from M0+M1 or from M0+...+Mn
      
      means_isotopologues_Raw <- c()
      
      means_M1M0 <- c()
      
      for (j in 1:length(MolecularSpecies_i)) {
        
        lipid_i_Raw <- Raw_File_i[grep(pattern = MolecularSpecies_i[j], x = rownames(File_i)),]
        
        lipid_i_corr <- File_i[grep(pattern = MolecularSpecies_i[j], x = rownames(File_i)),]
        
        M0_L <- rbind(lipid_i_Raw[1,], lipid_i_corr[2:dim(lipid_i_corr)[1],])
        
        means_isotopologues_Raw[j] <- mean(colSums(M0_L), na.rm = T)
        
        means_M1M0[j] <- mean(colSums(lipid_i_Raw[1:2,]), na.rm = T)
        
      }
      
      SteadyStatePoolsM1Mn[i,] <- means_isotopologues_Raw
      
      SteadyStatePoolsM1M0[i,] <- means_M1M0
      
      treatment_names[i] <- paste0(treatment_i, "_", rep_i)
      
    }
    
    colnames(SteadyStatePoolsM1Mn) <- MolecularSpecies
    colnames(SteadyStatePoolsM1M0) <- MolecularSpecies
    rownames(SteadyStatePoolsM1Mn) <- treatment_names
    rownames(SteadyStatePoolsM1M0) <- treatment_names
    
    if (is.null(SteadyStatePools_dir)) {
      
      write.csv(SteadyStatePoolsM1Mn, file = "SteadyStatePoolsM0Mn.csv")
      
      write.csv(SteadyStatePoolsM1M0, file = "SteadyStatePoolsM0M1.csv")
      
    } else {
      
      write.csv(SteadyStatePoolsM1Mn, file = paste0(SteadyStatePools_dir,
                                                    "/", "SteadyStatePoolsM0Mn.csv"))
      
      write.csv(SteadyStatePoolsM1M0, file = paste0(SteadyStatePools_dir,
                                                    "/", "SteadyStatePoolsM0M1.csv"))
    }
    
    
  } else {
    
    print("IsoCorrectoR directories must be three levels below, e.g.: ParentDir/xx/20...")
    
  }
  
  
}
