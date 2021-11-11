#' A function to manually subset consolidated data matrices into alike pixel groups using their distribution
#'
#' This function is interactive and needs to be run directly from the console to enhance the interactive mode. Thee function provides manual partition based on data distribution for molecular features that failed to be successfully subset using pvclust. The function will build a distribution density plot for each indexed molecular feature. The user then will be able to select how many partitions each feature has according to the histogram. Typically multi-modal distributions feature multiple peaks that can be captured this way. Subsequently the user will define and x-lim before which all pixels will be grabbed to integrate them into a subset; this procedure will be iterated for how many partitions the user previously defined. Finally, the function will return the mean partitions, just as the class discovery function before it.
#' @param kAssesmentOutput output from the previous KineticMSI function. Namely, kAssessmentMSI.R. The object must be a list of matrices, one matrix for each molecular feature of interest measured across replicates and treatments (i.e., rows in each matrix). kAssesmentOutput must come from "minDataset" in the parameters of the call.
#' @param indexVector defines the index of molecular features inside the kAssesment object that will be manually partitioned.
#' @param ZeroAction allows to define what happens with zeros within the entity data matrices. If null, the zeros are fully preserved, if "remove" the rows that only contain zeros are deleted and if "replace" all zeros are replaced by normally distributed randomly generated numbers in the scale of 10^-12 to 10^-13.
#' @keywords MSI Replicates Class Discovery Tracer Dynamics
#' @export
#' @examples
#' 
#' ...


kUserAssistedPartitionMSI <- function(kAssesmentOutput,
                                      indexVector,
                                      ZeroAction = c(NULL, "remove", "replace")) {
  
  ## Functions
  
  ### list2df
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  ## Main
  
  return_list <- list()
  
  names_return_list <- c()
  
  for (i in 1:length(indexVector)) {
    
    runner_df <- kAssesmentOutput[[indexVector[i]]]
    
    runner_name <- names(kAssesmentOutput)[indexVector[i]]
    
    ### Zero action taken on the datasets: make decision for full df or df wo zeros
    
    if(is.null(ZeroAction)) {
      
      if(length(which(colSds(runner_df) == 0)) > 0) {
        
        stop("ERROR: there are columns in the assembled matrices with NULL sd, please select an appropriate ZeroAction parameter in the function call")
        
      }
      
      runner_df <- as.matrix(runner_df)
      
    } else if(ZeroAction == "remove") {
      
      runner_df <- as.matrix(runner_df[,as.numeric(colSds(runner_df)) != 0])
      
    } else if(ZeroAction == "replace") {
      
      runner_df[runner_df == 0] <- rnorm(length(which(runner_df == 0)),
                                         mean = 0.0000000000001,
                                         sd = 0.00000000001)
      
      runner_df <- as.matrix(runner_df)
      
    }
    
    print(plot(density(colMeans(runner_df)), main = runner_name,
               xlab = "Enrichment Proportion",
               ylim = c(0, max(density(colMeans(runner_df))[["y"]])+5)))
      
    cat("\n")
    return_decision <- readline(prompt="How many partitions are there ? (ANS must be integer): ")
    cat("\n")
      
    colmeans_partitions <- colMeans(runner_df)
      
    out_vec <- c()
      
    for (j in 1:return_decision) {
        
      ### here ask for every partition and take the proportion of pixels of all of them
        
      cat("\n")
      return_decision2 <- readline(prompt = paste0("Please provide x-axis partition # ", j, " (ANS must be real number): "))
      cat("\n")
        
      ## test out of range partitions
        
      test_partitions <- which(colmeans_partitions <= as.numeric(return_decision2))
        
      if (length(test_partitions) > 0) {
          
        out_vec[j] <- return_decision2
          
      } else {
        
        cat("\n")
        return_decision2 <- readline(prompt = paste0("Try again, partition out of range! Please provide x-axis partition # ", j, " (ANS must be real number): "))
        cat("\n")
          
        out_vec[j] <- return_decision2
          
      }
    }
      
    total_pixel_Nr = length(colmeans_partitions)
      
    string_peak_j <- c()
      
    return_df <- matrix(NA,
                        nrow = length(rownames(runner_df)),
                        ncol = length(out_vec))
      
    for (j in 1:length(out_vec)) {
        
      ## pixel subset within each partitioned peak
        
      pixel_runner <- which(colmeans_partitions <= as.numeric(out_vec[j]))
        
      pixel_proportion = round((length(pixel_runner)/total_pixel_Nr)*100, 2)
        
      pixel_mean = round(mean(colmeans_partitions[pixel_runner]), 2)
        
      string_peak_j[j] <- paste0("Peak ", j, ",",
                                 " Pixel %: ", pixel_proportion, ",",
                                 " Mean: ", pixel_mean)
        
      ### return object
        
      return_df[,j] <- rowMeans(runner_df[,pixel_runner])
        
      ## replacement to grab each segment at a time with <=
        
      colmeans_partitions <- colmeans_partitions[-c(pixel_runner)]
        
    }
      
    png(filename = paste0(runner_name, "_AssistedPartition.png"))
      
    plot(density(colMeans(runner_df)), main = runner_name,
         xlab = "Enrichment Proportion",
         ylim = c(0, max(density(colMeans(runner_df))[["y"]])+(10*length(out_vec)))) +
        
     abline(v = out_vec, lty = 2, col = "red")
      
     legend("topleft", legend = string_peak_j)
      
     dev.off()
      
    ## return objects
      
    colnames(return_df) <- string_peak_j
    rownames(return_df) <- rownames(runner_df)
      
    return_list[[i]] <- return_df
      
    names_return_list[i] <- runner_name
    
  }
  
  names(return_list) <- names_return_list
  
  return(return_list)
  
}

