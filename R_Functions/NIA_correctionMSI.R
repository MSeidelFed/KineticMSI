






NIA_correctionMSI <- function(rm0_dir,
                              pattern = "_rm0",
                              ElementFile_dir,
                              MoleculeFile_dir,
                              out_dir){
  
  
  reps <- list.files(path = rm0_dir, pattern = pattern, all.files = FALSE,
                     full.names = T, recursive = FALSE,
                     ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  for (i in 1:length(reps)) {
    
    Enrichment <- IsoCorrectoR::IsoCorrection(MeasurementFile = reps[i],
                                              ElementFile = ElementFile_dir,
                                              MoleculeFile = MoleculeFile_dir,
                                              CorrectTracerImpurity = T,
                                              CorrectTracerElementCore = T,
                                              CalculateMeanEnrichment = T,
                                              UltraHighRes = F, 
                                              DirOut = out_dir,
                                              FileOut = strsplit(strsplit(reps[i],
                                                                          split = "/")[[1]][2],
                                                                 split = "\\.")[[1]][1],
                                              FileOutFormat = "csv",
                                              ReturnResultsObject = T,
                                              CorrectAlsoMonoisotopic = T)
    
  }
  
}
