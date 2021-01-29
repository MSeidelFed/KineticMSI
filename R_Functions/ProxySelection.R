




ProxySelection <- function(array_dir,
                           Treatments_dir,
                           Duplicate = F,
                           BatchCorr = F,
                           DirSteadyStatePoolsNL,
                           Factor = "Genotype",
                           Treatment_means = T) {
  
  ## functions
  
  ### batch correction function
  
  BatchCorrection <- function(array_dir,
                              Treatments_dir,
                              Duplicate = F,
                              Treatment_means){
    
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
    
    TreatmentFile <- read.csv(file = Treatments_dir, header = T, row.names = 1)
    
    if (Duplicate == T) {
      
      array <-  cbind(rbind(Condition = TreatmentFile$Sample,
                            Batch = TreatmentFile$Batch,
                            Treatment = rownames(TreatmentFile),
                            t(BatchCor_array)),
                      rbind(Condition = TreatmentFile$Sample,
                            Batch = TreatmentFile$Batch,
                            Treatment = rownames(TreatmentFile),
                            t(BatchCor_array)))
      
      colnames(array) <- c(rownames(BatchCor_array), rownames(BatchCor_array))
      
    } else {
      
      array <-  rbind(Condition = TreatmentFile$Sample,
                      Batch = TreatmentFile$Batch,
                      Treatment = rownames(TreatmentFile),
                      t(BatchCor_array))
      
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
  
  ## Main
  
  if (BatchCorr == T) {
    
    NC_Transformed_M0M1 <- BatchCorrection(array_dir = array_dir,
                                           Treatments_dir = Treatments_dir, 
                                           Duplicate = F,
                                           Treatment_means = Treatment_means)
    
    NC_Transformed_M0Mn <- BatchCorrection(array_dir = array_dir,
                                           Treatments_dir = Treatments_dir, 
                                           Duplicate = F,
                                           Treatment_means = Treatment_means)
    
    
  }
  
  
  ### Importing non-labelled pools
  
  SteadyStatePools_NL <- t(read.csv(file = DirSteadyStatePoolsNL,
                                    header = T, row.names = 1))
  
  TreatmentFile <- read.csv(file = Treatments_dir, header = T, row.names = 1)
  
  ### LvsNL pools
  
  DF_L = NC_Transformed_M0M1
  
  DF_NL = SteadyStatePools_NL
  
  Treatment = as.factor(TreatmentFile[,Factor])
  
  test_col <- colnames(DF_L) == colnames(DF_NL) ## all names are not equal
  
  test_row <- rownames(DF_L) == rownames(DF_NL) ## all names are not equal
  
  
  if (length(unique(test_row)) > 1) {
    
    if (which.min(c(length(rownames(DF_L)), length(rownames(DF_NL)))) == 1) {
      
      DF_NL <- DF_NL[rownames(DF_L),]
      
    } else {
      
      DF_L <- DF_L[rownames(DF_NL),]
      
    }
    
  }
  
  if (length(unique(test_col)) > 1) {
    
    print("Treatments must be the same in both matrices")
    
  }
  
  ratio_L_NL <- DF_L / DF_NL
  
  ratio_L_NL[is.na(ratio_L_NL)] <- 0
  ratio_L_NL[ratio_L_NL == Inf] <- 0
  ratio_L_NL[ratio_L_NL == -Inf] <- 0
  
  if (Treatment_means == T) {
    
    out_mat <- matrix(NA, nrow = length(unique(Treatment)), ncol = dim(ratio_L_NL)[1])
    
    for (i in 1:dim(ratio_L_NL)[1]) {
      
      out_mat[,i] <- aggregate(ratio_L_NL[i,] ~ Treatment, FUN = mean)[,2]
      
    }
    
    rownames(out_mat) <- aggregate(ratio_L_NL[i,] ~ Treatment, FUN = mean)[,1]
    colnames(out_mat) <- rownames(ratio_L_NL)
    
    ratio_L_NL <- t(out_mat)
    
    type = gsub("s\\d+_", "", unique(Treatment))
    
  } else {
    
    type = gsub("s\\d+_", "", Treatment)
    
  }
  
  ha = HeatmapAnnotation(df = data.frame(type = type))
  
  h<- Heatmap(matrix = ratio_L_NL,
              name = "Intensities",
              km = 1,
              col = colorRamp2(c(0.5, 0.75, 1, 1.25, 1.5),
                               c("yellow", "white", "white", "white", "grey")),
              show_column_names = F,
              show_row_names = T,
              cluster_columns = F,
              cluster_rows = T,
              top_annotation = ha)
  print(h)
  
  return(ratio_L_NL[row_order(h),])
  
  
}