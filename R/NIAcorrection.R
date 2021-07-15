#' A function to correct natural isotopic abundances (NIA) from kineticMSI datasets or multiple csv files with the right format
#'
#' This function allows you to correct isotopologue envelopes for NIA inheriting all specifications from the R package [IsoCorrectoR](https://www.bioconductor.org/packages/release/bioc/html/IsoCorrectoR.html). The function calculates the percentage of enrichment of a defined stable isotope as well as other important values that reflect tracer dynamics within enrichment experiments. The function takes an entire directory and grabs all the csv files contained within in a recursive manner. Subsequently, the function generates csv files with each of the relevant returns in the *kinetic*MSI context. Each column in the input table "MeasurementFile.csv" belongs to a single coordinate on the original image where the isotopologues could be measured and mined out.
#' @param MeasurementFileDir directory where the input files are stored.
#' @param pattern defaults to "_rm0", that is, it grabs the files produced using rmNullPixel. It is used to define the pattern on which input files are looked for in the MeasurementFileDir
#' @param SubSetReps defaults to FALSE. Allows to subset the file list found in MeasurementFileDir.
#' @param ElementFileDir directory where the element input file is stored (see IsoCorrectoR documentation for input requirements).
#' @param MoleculeFileDir directory where the molecule input file is stored (see IsoCorrectoR documentation for input requirements).
#' @param kCorrectTracerImpurity defaults to TRUE.see IsoCorrectoR documentation for usage details.
#' @param kCorrectTracerElementCore defaults to TRUE. see IsoCorrectoR documentation for usage details.
#' @param kCalculateMeanEnrichment defaults to TRUE. see IsoCorrectoR documentation for usage details.
#' @param kCorrectAlsoMonoisotopic defaults to TRUE. see IsoCorrectoR documentation for usage details.
#' @param kUltraHighRes defaults to FALSE. see IsoCorrectoR documentation for usage details.
#' @param kCalculationThreshold defaults to IsoCorrectoR predefined value. see IsoCorrectoR documentation for usage details.
#' @param kCalculationThreshold_UHR defaults to IsoCorrectoR predefined value. see IsoCorrectoR documentation for usage details.
#' @param verbose defaults to FALSE. When TRUE it returns to the console the progression across the input files. Thus the parameter is meant to allow users to spot errors in the input files.
#' @param outdir defines the new directory that will contain IsoCorrectoR outputs, if the direction does not exist it will be automatically created. The new directory contains subdirectories named with the exact time of the run according to IsoCorrectoR conventions.
#' @keywords KineticMSI NIA correction IsoCorrectoR
#' @export
#' @examples
#' 
#' ...


NIAcorrection <- function(MeasurementFileDir,
                          pattern = "_rm0",
                          SubSetReps = F,
                          ElementFileDir,
                          MoleculeFileDir,
                          kCorrectTracerImpurity = TRUE,
                          kCorrectTracerElementCore = TRUE,
                          kCalculateMeanEnrichment = TRUE,
                          kCorrectAlsoMonoisotopic = TRUE,
                          kUltraHighRes = FALSE,
                          kCalculationThreshold = 10^-8,
                          kCalculationThreshold_UHR = 8,
                          verbose = FALSE,
                          outdir){
  
  
  reps <- list.files(path = MeasurementFileDir, pattern = pattern,
                     all.files = FALSE, full.names = TRUE, recursive = TRUE,
                     ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  
  if (SubSetReps == T) {
    
    reps = SubSetInputFiles(reps = reps)
    
  }
  
  for (i in 1:length(reps)) {
    
    IsoCorrectoR::IsoCorrection(MeasurementFile = reps[i],
                                ElementFile = ElementFileDir,
                                MoleculeFile = MoleculeFileDir,
                                CorrectTracerImpurity = kCorrectTracerImpurity,
                                CorrectTracerElementCore = kCorrectTracerElementCore,
                                CalculateMeanEnrichment = kCalculateMeanEnrichment,
                                UltraHighRes = kUltraHighRes,
                                ReturnResultsObject = FALSE,
                                CorrectAlsoMonoisotopic = kCorrectAlsoMonoisotopic, 
                                DirOut = outdir,
                                FileOut = strsplit(strsplit(reps[i],
                                                            split = "/")[[1]][2],
                                                   split = "\\.")[[1]][1],
                                FileOutFormat = "csv", 
                                CalculationThreshold = kCalculationThreshold,
                                CalculationThreshold_UHR = kCalculationThreshold_UHR,
                                verbose = verbose)
    
  }
  
}
