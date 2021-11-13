#' A function to subset input files
#'
#' This function allows KineticMSI functions to subset lists of files at multiple steps.
#' @param reps must be a list of replicate files.
#' @param SubSetReps defaults to TRUE. Allows to subset the file list found in MeasurementFileDir.
#' @keywords KineticMSI Preprocessing
#' @export
#' @examples
#' 
#' ...




SubSetInputFiles <- function(reps,
                             SubSetReps = T){
  
  for (rep in 1:length(reps)) {
    
    cat("\n", reps[rep])
    
  }
  
  while (SubSetReps == T) {
    
    grepPattern <- readline("Please enter new pattern to subset current reps (text): ")
    
    reps = reps[grep(grepPattern, reps)]
    
    cat("Final Reps: ")
    cat("\n")
    
    for (rep in 1:length(reps)) {
      
      cat("\n", reps[rep])
      
    }
    
    SubSetReps <- readline("Do you want to subset reps even more ?? (T or F): ")
    
  }
  
  return(reps)
  
}