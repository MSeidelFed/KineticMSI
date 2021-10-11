#' A function that allows class comparison in replicated KineticMSI datasets
#'
#' This function allows KineticMSI users to use classical class comparison algorithms such as ANOVA + Tukey HSD or generalized linear models (GLMs) + FDR correction. The input for the function comes from the previous kAssessmentMSI.R output. The function outputs graphics embedded in PDF files detailing the results from the class comparison algorithms. Additionally, a table with all the results is returned to the environment when the ouput is assigned to an object. 
#' @param kAssesmentOutput output from the previous KineticMSI function. Namely, kAssessmentMSI.R. The object must be a list of matrices, one matrix for each molecular feature of interest measured across replicates and treatments (i.e., rows in each matrix).
#' @param factorVector character vector that needs to define the treatments using the same nomenclature and order in the naming scheme as the row names in matrices within the input files (follow the exemplary KineticMSI data for details).
#' @param PDFname defaults to "test". Defines the name of the output PDF files that will be produce as outcome in the current working directory.
#' @param returnGLMplots defaults to TRUE. Defines whether GLM plots are returned.
#' @param patternGLMplot defaults to "Q values". Defines whether Q or P values are used to return significance on top of the boxplots within the PDF containing the GLM results.
#' @param ylabTukey defaults to NULL. Character vector that defines the text on y-axes from Tukey HSD plots.
#' @param xlabTukey defaults to NULL. Character vector that defines the text on x-axes from Tukey HSD plots.
#' @param ylabGLM defaults to NULL. Character vector that defines the text on y-axes from GLM plots.
#' @param xlabGLM defaults to NULL. Character vector that defines the text on x-axes from GLM plots.
#' @keywords MSI Replicates Class Comparison Tracer Dynamics
#' @export
#' @examples
#' 
#' ...


kClassComparisonMSI <- function(kAssesmentOutput,
                                factorVector,
                                PDFname = "test",
                                returnGLMplots = TRUE,
                                patternGLMplot = c("Q values", "P values"),
                                ylabTukey = NULL,
                                xlabTukey = NULL,
                                ylabGLM = NULL,
                                xlabGLM = NULL) {
  
  ## functions needed
  
  GLMplot_CustomizedMSI <- function(Variable,
                                    factor_glm,
                                    P_values,
                                    patternGLMplot = patternGLMplot,
                                    MainTitle = "",
                                    ylabGLMs = ylabGLM,
                                    xlabGLMs = xlabGLM){
    
    colorVector = rand_color(length(levels(factor_glm)))
    
    FactorNumber <- length(levels(factor_glm))
    
    
    #### Using the glm P-values
    
    color_plot <- rep(colorVector, (length(P_values)/FactorNumber))
    
    #### assigning stars to P-values, ° = +inf:0.1, * = 0.1:0.05, ** 0.05:0.01, *** 0.01:0
    
    stars_out <- c()
    
    for (i in 1:length(P_values)) {
      
      if (P_values[i] < 0.1 & P_values[i] > 0.05) {
        runner = "*"
      } else if (P_values[i] < 0.05 & P_values[i] > 0.01) {
        runner = "**"
      } else if (P_values[i] < 0.01) {
        runner = "***"
      } else {
        runner = "°"
      }
      
      stars_out[i] <- runner 
      
    }
    
    stars_out[1] <- "C"
    
    #### plotting
    
    boxplot(Variable ~ factor_glm,
            ylim = c(0, (max(Variable) + mean(Variable))), 
            main = MainTitle,
            las = 2,
            col = color_plot,
            xlab = xlabGLMs,
            ylab = ylabGLMs)
    
    text(x = c(1:length(unique(factor_glm))),
         y = (max(Variable) + mean(Variable)),
         labels = round(P_values, 2))
    
    text(x = c(1:length(unique(factor_glm))),
         y = (max(Variable) + (mean(Variable))/2),
         labels = stars_out)
    
    
    #### Add data points (https://www.r-graph-gallery.com/96-boxplot-with-jitter.html)
    
    data <- data.frame(names = factor_glm, value = Variable)
    
    #### Add data points (https://www.r-graph-gallery.com/96-boxplot-with-jitter.html)
    
    mylevels <- levels(data$names)
    levelProportions <- summary(data$names)/nrow(data)
    
    for(i in 1:length(mylevels)){
      
      thislevel <- mylevels[i]
      thisvalues <- data[data$names==thislevel, "value"]
      
      ##### take the x-axis indices and add a jitter, proportional to the N in each level
      myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
      
      points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
      
    }
    
  }
  
  ## main
  
  input_list <- kAssesmentOutput
  
  factorTest <- c()
  
  for (i in 1:length(unique(factorVector))) {
    
    factorTest[i] <- length(unique(grep(unique(factorVector)[i], rownames(input_list[[1]])) == grep(unique(factorVector)[i], factorVector)))
    
  }
  
  FinalTest <- unique(factorTest) == 1
  
  if (FinalTest == T) {
    
    cat(".....\n")
    cat("input factor vetor coincides with rownames in input matrices")
    cat("\n")
    cat(".....\n")
    
    normality_test <- c()
    
    p_value <- c()
    
    E_size <- c()
    
    MC_vector <- c()
    
    test_mat <- matrix(NA,
                       nrow = dim(input_list[[1]])[1],
                       ncol = length(input_list))
    
    pdf(paste0(PDFname, ".pdf"))
    par(mar = c(5, 10, 3, 2))
    par(mfrow = c(3,1))
    
    for (i in 1:length(input_list)) {
      
      cat("...\n")
      cat(paste0("Running on Entity No. ", i, "\n"))
      cat("...\n")
      
      if (class(input_list[[i]]) == "matrix") {
        
        ## NA positions from failed previous assays
        
        if(is.null(dim(input_list[[i]]))) {
          
          ## Empty mats with no columns
          
          normality_test[i] <- NA
          
          p_value[i] <- NA
          
          E_size[i] <- NA
          
          ## mean comparison (Log2FC)
          
          MC_vector[i] <- NA
          
          test_mat[,i] <- rep(NA, dim(test_mat)[1])
          
        } else if (dim(input_list[[i]])[2] == 0) {
          
          ## Empty mats with no columns
          
          normality_test[i] <- NA
          
          p_value[i] <- NA
          
          E_size[i] <- NA
          
          ## mean comparison (Log2FC)
          
          MC_vector[i] <- NA
          
          test_mat[,i] <- rep(NA, dim(test_mat)[1])
          
        } else {
          
          ## rowmeans
          
          runner <- rowMeans(input_list[[i]])
          
          runner[runner == 0] <- 0.000001
          
          ## Tukey-HSD
          
          TukeyCustomized(variable = runner,
                          factor = as.factor(factorVector),
                          MainTitle = names(input_list)[i],
                          returnObject = "Letters")
          
          ## saving data for GLM test
          
          test_mat[,i] <- runner
          
          ## defining factors
          
          runner_means <- rowMeans(input_list[[i]])
          
          x_ks <- melt(input_list[[i]][grep(unique(factorVector)[1],
                                            names(rowMeans(input_list[[i]]))),])$value
          
          y_ks <- melt(input_list[[i]][grep(unique(factorVector)[2],
                                            names(rowMeans(input_list[[i]]))),])$value
          
          ## mean comparison (Log2FC)
          
          MC_vector[i] <- log2(mean(x_ks) / mean(y_ks))
          
          ## KS test
          
          normality_test[i] <- paste0(suppressWarnings(ks.test(x = x_ks, y = y_ks)[["method"]]),
                                      ": p = ",
                                      suppressWarnings(ks.test(x = x_ks, y = y_ks)[["p.value"]]))
          
          p_value[i] <- suppressWarnings(ks.test(x = x_ks, y = y_ks)[["p.value"]])
          
          
          E_size[i] <- as.numeric(as.character(cohen.d(d = x_ks, f = y_ks))[3])
          
       
          ## plot
          
          y0_limit <- min(density(x_ks)[["y"]], density(y_ks)[["y"]])
          y_limit <- max(density(x_ks)[["y"]], density(y_ks)[["y"]])
          x0_limit <- min(density(x_ks)[["x"]], density(y_ks)[["x"]])
          x_limit <- max(density(x_ks)[["x"]], density(y_ks)[["x"]])
          
          
          plot(density(y_ks), col = "blue",
               xlim = c(x0_limit, x_limit),
               ylim = c(y0_limit, y_limit),
               main = c(names(input_list)[i],
                        paste0("KS-Pvalue: ", p_value[i]),
                        paste0("Cohen's d value: ", E_size[i])))
          lines(density(x_ks), col = "red", lty = 1)
          legend(x = x0_limit, y = y_limit,
                 legend =  c(unique(factorVector)[2],
                             unique(factorVector)[1]),
                 col = c("blue", "red"), lty=1, cex=0.8)
          
        }
        
      } else {
        
        normality_test[i] <- NA
        
        p_value[i] <- NA
        
        E_size[i] <- NA
        
        ## mean comparison
        
        MC_vector[i] <- NA
        
        test_mat[,i] <- rep(NA, dim(test_mat)[1])
        
      }
    }
    
    dev.off()
    
    ## object names
    
    colnames(test_mat) <- names(input_list)
    rownames(test_mat) <- rownames(input_list[[1]])
    
    names(normality_test) <- names(input_list)
    
    ## stat test
    
    Pvalues_test <- OmicsUnivariateStats(class_comparison_mat = test_mat,
                                         Factor1 = factorVector,
                                         Contrast = F,
                                         TukeyReturns = "Letters", 
                                         returnObject = "OmicsTests")
    
    Original_NAs <- which(is.na(as.numeric(colMeans(test_mat))))
    
    if(length(Original_NAs) > 0) {
      
      for (i in 1:length(Original_NAs)) {
        
        Pvalues_tests <- rbind(Pvalues_test[1:(Original_NAs[i]-1),],
                               rep(NA, dim(Pvalues_test)[2]),
                               Pvalues_test[(Original_NAs[i]):(dim(Pvalues_test)[1]),])
        
      }
      
    } else {
      
      Pvalues_tests <- Pvalues_test
      
    }
    
    ## out_object
    
    p_Adjustment <- cbind(normality_test,
                          KS_pAdjusted = p.adjust(p_value, method = "BH"),
                          Cohensd = E_size,
                          Log2FC = MC_vector,
                          t(test_mat),
                          Pvalues_tests)
    
    ## eliminate NA entries for out_object if present from previous step
    
    if(length(which(is.na(p_Adjustment[,1]))) > 0) {
      
      p_Adjustment = p_Adjustment[-c(which(is.na(p_Adjustment[,1]))),]
      
    }
    
    ## returning plots
    
    if(returnGLMplots == T){
      
      pdf(paste0(PDFname, "GLM.pdf"))
      par(mar = c(5, 10, 3, 2))
      
      for (i in 1:dim(p_Adjustment)[1]) {
        
        GLMplot_CustomizedMSI(Variable = as.numeric(p_Adjustment[i,5:(5+(dim(test_mat)[1])-1)]),
                              factor_glm = as.factor(factorVector),
                              P_values =as.numeric(p_Adjustment[i,grep(patternGLMplot, colnames(p_Adjustment))]),
                              MainTitle = rownames(p_Adjustment)[i])
        
      }
      
      dev.off()
      
    }
    
  } else {
    
    stop("ERROR: input factor vetor does not coincide with rownames in input matrices")
    
  }
  
  return(p_Adjustment)
  
}
