



UserAssistedPartition <- function(DatAssesmentObject,
                                  EvalIndividualReps = F,
                                  index_vector,
                                  RemoveZeros = F) {
  
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
  
  for (i in 1:length(index_vector)) {
    
    runner_df <- DatAssesmentObject[[index_vector[i]]]
    
    runner_name <- names(DatAssesmentObject)[index_vector[i]]
    
    ## make decission for full df or df wo zeros
    
    if (RemoveZeros == T) {
      
      test_zero_means <- which(colMeans(runner_df) == 0)
      
      if (length(test_zero_means) > 0) {
        
        runner_df[,-c(test_zero_means)]
        
      } else {
        
        cat("\n")
        cat("No columns with mean Zero were found", runner_name)
        cat("\n")
        
      }
      
    }
    
    if (EvalIndividualReps == T) {
      
      for (j in 1:dim(runner_df)[1]) {
        
        file_name_index <- length(strsplit(apply(list2df(strsplit(rownames(runner_df), "\\."))[1,], 2, as.character), "/")[1]$V1)
        
        file_names <- apply(list2df(strsplit(apply(list2df(strsplit(rownames(runner_df), "\\."))[1,], 2, as.character), "/"))[file_name_index,], 2, as.character)
        
        name_j <- paste0(runner_name, "_" , file_names[j])
        
        plot(density(runner_df[j,]), main = name_j,
             xlab = "Enrichment Proportion",
             ylim = c(0, max(density(runner_df[j,])[["y"]])+5))
        
      }
      
      
    } else {
      
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
          
        }
      }
      
      total_pixel_Nr = length(colmeans_partitions)
      
      string_peak_j <- c()
      
      for (j in 1:length(out_vec)) {
        
        ## pixel subset within each partitioned peak
        
        pixel_runner <- which(colmeans_partitions <= as.numeric(out_vec[j]))
        
        pixel_proportion = round((length(pixel_runner)/total_pixel_Nr)*100, 2)
        
        pixel_mean = round(mean(colmeans_partitions[pixel_runner]), 2)
        
        string_peak_j[j] <- paste0("Peak ", j, ",",
                                   " Pixel %: ", pixel_proportion, ",",
                                   " Mean: ", pixel_mean)
        
        ## replacement to grab each segment at a time with <=
        
        colmeans_partitions <- colmeans_partitions[-c(pixel_runner)]
        
      }
      
      png(filename = paste0(runner_name, "_AssistedPartition.png"))
      
      plot(density(colMeans(runner_df)), main = runner_name,
           xlab = "Enrichment Proportion",
           ylim = c(0, max(density(colMeans(runner_df))[["y"]])+(10*length(out_vec)))) +
        
        abline(v = out_vec, lty = 2, col = "red")
      
      legend("topleft", legend=string_peak_j)
      
      dev.off()
      
    }
    
  }
  
}

