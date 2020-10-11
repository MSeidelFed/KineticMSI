









LvsNLpools <- function(DF_L,
                       DF_NL,
                       Treatment,
                       Treatment_means = T) {
  
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
              col = colorRamp2(as.numeric(summary(melt(ratio_L_NL)$value)[-3]),
                               c("yellow", "white", "white", "white", "grey")),
              show_column_names = F,
              show_row_names = T,
              cluster_columns = F,
              cluster_rows = T,
              top_annotation = ha)
  print(h)
  
}
