







EnrichmentProportions_comparison <- function(FilesPath,
                                             pattern = "MeanEnrichment.csv",
                                             ProportionOperator = c("equal", "less", "greater"),
                                             ProportionLimit = 0,
                                             factorVector,
                                             k = 5,
                                             n_boot = 100,
                                             ClustMethod = "average",
                                             returnProprotionsHeatmap = T) {
  
  ## main
  
  reps <- list.files(path = FilesPath, pattern = pattern, all.files = FALSE,
                     full.names = T, recursive = T,ignore.case = FALSE,
                     include.dirs = FALSE, no.. = FALSE)
  
  out_mat <- matrix(NA, nrow = 85, ncol = length(reps))
  
  for (j in 1:length(reps)) {
    
    test <- read.csv(file = reps[j], header = T, row.names = 1)
    
    out_vec <- c()
    
    #count = 0
    
    for (i in 1:dim(test)[1]) {
      
      #count = count +1
      
      #print(count)
      
      if (ProportionOperator == "equal") {
        
        out_vec[i] <- length(which(test[i,] == ProportionLimit))/dim(test)[2]
        
      } else if (ProportionOperator == "greater") {
        
        out_vec[i] <- length(which(test[i,] > ProportionLimit))/dim(test)[2]
        
      } else if (ProportionOperator == "less") {
        
        out_vec[i] <- length(which(test[i,] < ProportionLimit))/dim(test)[2]
        
      }
      
    }
    
    out_mat[,j] <- out_vec
    
  }
  
  rownames(out_mat) <- rownames(test)
  colnames(out_mat) <- factorVector
  
  out_df <- as.data.frame(t(out_mat))
  
  P_values <- c()
  
  for (i in 1:dim(out_df)[2]) {
    
    P_values[i] <- summary(glm(out_df[,i] ~ factorVector))$coefficients[2,"Pr(>|t|)"]
    
  }
  
  test_mean <- aggregate(. ~ as.factor(factorVector), out_df, mean)
  
  test_num <- as.matrix(test_mean[,2:dim(test_mean)[2]])
  rownames(test_num) <- test_mean[,1]
  
  Padj_values <- as.matrix(p.adjust(P_values, "bonferroni"))
  rownames(Padj_values) <- colnames(out_df)
  
  P_values <- as.matrix(P_values)
  rownames(P_values) <- colnames(out_df)
  
  ## set to 1 NaN (come from no 0 in the matrices - signature region lipids)
  
  Padj_values[is.na(Padj_values)] <- 1
  
  P_values[is.na(P_values)] <- 1
  
  ## set to null names without sig padj values 
  
  rownames(Padj_values)[which(Padj_values > 0.1)] <- "."
  
  rownames(P_values)[which(P_values > 0.1)] <- "."
  
  if (length(unique(rownames(Padj_values))) > 1) {
    
    P_heatmap = Padj_values
    
  } else {
    
    P_heatmap = P_values
    
  }
  
  type = gsub("s\\d+_", "", unique(factorVector))
  
  ha = HeatmapAnnotation(df = data.frame(type = type))
  
  if (returnProprotionsHeatmap == T) {
    
    pdf("Proportions_Heatmap.pdf")
    
    Scal_Heatmap <- Heatmap(matrix = as.matrix(t(test_num)),
                            name = "Intensities",
                            km = k, ## five K-means
                            row_km_repeats = n_boot, ## repeats to get consensus
                            col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                                             c("white",
                                               "yellow",
                                               "darkgoldenrod1",
                                               "violet",
                                               "purple")),
                            top_annotation = ha ,
                            show_column_names = F,
                            cluster_columns = F,
                            cluster_rows = T,
                            clustering_method_rows = ClustMethod,
                            clustering_distance_rows = "euclidean") +
      
      Heatmap(as.matrix(P_heatmap),
              col = colorRamp2(c(0, 0.1, 0.11, 1),
                               c("black","grey","white","white")),
              show_row_names = T)
    
    print(Scal_Heatmap)
    
    dev.off()
    
  }
  
  
  return_mat <- cbind(t(out_df),
                      P_values = P_values[,1],
                      Padj_values = Padj_values[,1])
  
  return(return_mat)
  
  
}
