### AcqTimeRandomSampling_MSI

AcqTimeRandomSampling_MSI <- function(FilesPath,
                                      pattern = "MeanEnrichment.csv",
                                      feature_choice = 0) {
  
  #### choosing a lipid
  
  reps <- list.files(path = FilesPath, pattern = pattern, all.files = FALSE,
                     full.names = T, recursive = T,
                     ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  if (feature_choice == 0) {
    
    file_test <- read.csv(file = reps[1], header = T, row.names = 1)
    
    lipid_Nr <- unique(rownames(file_test))
    
    print(lipid_Nr)
    
    lipid_choice <- readline(prompt="Enter index of lipid of interest: ")
    
  } else  { lipid_choice = feature_choice }
  
  #### defining the minimum set length
  
  out_vector <- c()
  
  for (i in 1:length(reps)) {
    
    file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
    
    out_vector[i] <- dim(file_runner)[2]
    
  }
  
  min_dataset <- min(out_vector)
  
  #### sampling randomly but conserving spatial relationships among pixels
  
  out_single_lipid <- matrix(NA, length(reps), min_dataset)
  
  for (i in 1:length(reps)) {
    
    file_runner <- read.csv(file = reps[i], header = T, row.names = 1)
    
    lipid_runner <- as.matrix(sample(file_runner[as.numeric(lipid_choice),],
                                     min_dataset,
                                     replace = F))
    
    df2sort <- as.data.frame(t(rbind(apply(list2df(strsplit(colnames(lipid_runner),
                                                            "\\."))[2,],
                                           2,
                                           as.numeric),
                                     lipid_runner)))
    
    colnames(df2sort) <- c("pixelNr", "Enr")
    
    SortedEnrichment <- df2sort[order(df2sort$pixelNr),2]
    
    out_single_lipid[i,] <- SortedEnrichment
    
  }
  
  rownames(out_single_lipid) <- reps
  
  return(out_single_lipid)
  
}
