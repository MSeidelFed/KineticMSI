




IncorporationProxys <- function(Parent_dir) {
  
  if (length(strsplit(Parent_dir,split = "/")[[1]]) == 2) {
    
    Enrichment_files <- list.files(path = Parent_dir, pattern = "Corrected.csv",
                                   full.names = T, recursive = T)
    
    for (i in 1:length(Enrichment_files)) {
      
      Enrichment_file_runner <- Enrichment_files[i]
      
      ### runner indexes
      
      File_i <- read.csv(file = Enrichment_file_runner, header = T, row.names = 1)
      
      name_i <- strsplit(strsplit(Enrichment_file_runner, split = "/")[[1]][5], split = "\\.")[[1]][1]
      
      out_path_i <- paste0(Parent_dir, "/", strsplit(Enrichment_file_runner, split = "/")[[1]][4], "/")
      
      ### generating new files
      
      M0s <- File_i[grep(pattern = "_0", x = rownames(File_i)),]
      
      write.csv(M0s, paste0(out_path_i, name_i, "M0s.csv"))
      
      M1s <- File_i[grep(pattern = "_1", x = rownames(File_i)),]
      
      write.csv(M1s, paste0(out_path_i, name_i, "M1s.csv"))
      
      M1M0ratio <- M1s/M0s
      
      write.csv(M1M0ratio, paste0(out_path_i, name_i, "M1M0ratios.csv"))
      
      M1fraction <- M1s/M0s+M1s
      
      write.csv(M1fraction, paste0(out_path_i, name_i, "M1fraction.csv"))
      
    }
      
    } else {
      
      print("IsoCorrectoR directories must be three levels below, e.g.: ParentDir/xx/20...")
      
    }
    
    
}

