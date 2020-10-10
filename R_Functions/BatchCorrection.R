




BatchCorrection <- function(array_dir = "Data/Steady_state_pools/SteadyStatePoolsM1Mn.csv",
                            Treatments_dir = "Data/Treatments.csv",
                            Duplicate = F){
  
  ## functions
  
  Class_Distribution <- function(in_mat,
                                 Treatments, 
                                 plots = NULL, 
                                 returnTable = F) {
    
    mydata4.1 <- as.data.frame(t(in_mat))
    colnames(mydata4.1) <- seq(1, dim(mydata4.1)[2])
    
    mydata4.2 <- cbind(Treatments, mydata4.1)
    colnames(mydata4.2) <- c("X1_1", seq(1,dim(mydata4.1)[2]))
    
    mydata4.3 <- melt(mydata4.2, id.vars= seq(1), measure.vars= seq(2, (dim(mydata4.1)[2]+1)))
    
    mydata4.4 <- mydata4.3[order(mydata4.3$X1_1),]
    
    if(plots == T) {
      
      print(ggplot(mydata4.4, aes(x = as.numeric(value), y = X1_1)) + 
              xlim(-summary(as.numeric(mydata4.4$value))["Mean"],
                   summary(as.numeric(mydata4.4$value))["Mean"]) +
              ggtitle ("Test") +
              geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
              theme_ridges(font_size = 10, grid = TRUE) +
              theme(axis.title.y = element_blank()))
      
      
      print(ggplot(mydata4.4, aes(x = as.numeric(value), fill = X1_1)) +
              geom_density(alpha = 0.2, position = "identity") +
              xlim(-summary(as.numeric(mydata4.4$value))["Mean"],
                   summary(as.numeric(mydata4.4$value))["Mean"]) +
              ggtitle ("Test"))
      
    }
    
    if(returnTable == T) {
      
      return(mydata4.4)
      
    }
    
  }
  
  ## data
  
  ### building array
  
  BatchCor_array <- read.csv(file = array_dir, header = T, row.names = 1)
  
  TreatmentFile <- read.csv(file = Treatments_dir, header = T)
  
  if (Duplicate == T) {
    
    array <-  cbind(rbind(Condition = TreatmentFile$Sample,
                          Batch = TreatmentFile$Batch,
                          Treatment = rownames(TreatmentFile),
                          t(SteadyStatePools_LM0Mn)),
                    rbind(Condition = TreatmentFile$Sample,
                          Batch = TreatmentFile$Batch,
                          Treatment = rownames(TreatmentFile),
                          t(SteadyStatePools_LM0Mn)))
    
    colnames(array) <- c(rownames(BatchCor_array), rownames(BatchCor_array))
    
  } else {
    
    array <-  rbind(Condition = TreatmentFile$Sample,
                    Batch = TreatmentFile$Batch,
                    Treatment = rownames(TreatmentFile),
                    t(SteadyStatePools_LM0Mn))
    
    colnames(array) <- rownames(BatchCor_array)
    
  }
  
  ### preparing the data in the correct format
  
  array_intensities <- array[4:dim(array)[1],]
  colnames(array_intensities) <- colnames(array)
  
  Numeric_array <-apply(X = array_intensities, MARGIN = 2, FUN = as.numeric)
  
  Numeric_array[is.na(Numeric_array)] <- 0
  rownames(Numeric_array) <- rownames(array_intensities)
  colnames(Numeric_array) <- as.character(colnames(array_intensities))
  
  Numeric_array <- Numeric_array[as.numeric(rowMeans(Numeric_array)) != 0,]

  batch <- array["Batch",]
  
  ## main
  
  Class_Distribution(in_mat = Numeric_array, 
                     Treatments = batch, 
                     plots = T, returnTable = F)
  
  batch_necessity <- readline(prompt="Do you want to proceed with Batch correction? (Y/N)")
  
  
  if (batch_necessity == "Y") {
    
    NC_Normalized <- ComBat(dat =  Numeric_array, batch = batch,
                            mod = NULL,  par.prior = T, prior.plots = T)
    
    colnames(NC_Normalized) <- colnames(Numeric_array)
    rownames(NC_Normalized) <- rownames(Numeric_array)
    
    Class_Distribution(in_mat = NC_Normalized, 
                       Treatments = batch, 
                       plots = T, returnTable = F)
    
    return(NC_Normalized)
    
  } else {
    
    stop("Exit")
    
  }
  
}


