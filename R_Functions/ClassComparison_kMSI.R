


ClassComparison_kMSI <- function(ClassDiscoveryList,
                                 factorVector,
                                 colorVector = c("tomato3", "peachpuff3"),
                                 FactorNumber = 2,
                                 ClassComparison = c(FALSE, "GLM", "ANOVA"),
                                 control = NULL,
                                 confoundingFactor = NULL) {
  
  ## functions needed for the main
  
  #### Tukey function
  
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  
  #### handling errors
  
  show_condition <- function(code) {
    tryCatch(code,
             error = function(c) "error",
             message = function(c) "message"
    )
  }
  
  ### list2df
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  

  ### mean comparison
  
  pdf(file = "BoxPlot_with_Pvalues.pdf")
  
  cat("Doing class comparison tests...", "\n")
  
  list_out <- ClassDiscoveryList
  
  P_Values <- list()
  
  count = 0
  
  for (m in 1:length(list_out)) {
    
    count = count + 1
    
    #print(count)
    
    if (is.na(names(list_out)[m])){
      
      #### no mat
      
    } else {
      
      if (is.null(control)) {
        
        #### when not defined, the control is defined to be the cluster with the lowest enrichment
        
        control <- which(min(colMeans(list_out[[m]])) == colMeans(list_out[[m]]))
        
        mat_to_melt <- list_out[[m]]
        
        ##### adding _ to controls to make them first alphabetically
        
        colnames(mat_to_melt)[control] <- paste0("_", colnames(mat_to_melt)[control])
        
      } else if (class(control) == "character") {
        
        ##### adding _ to controls to make them first alphabetically
        
        factorVector[grep(control, factorVector)] <- paste0(rep("_", length(grep(control, factorVector))), grep(control, factorVector, value = T))
        
        mat_to_melt <- list_out[[m]]
        
      } else {
        
        cat("ERROR... control must be a character")
        
      }
      
      #### melting each mat to get the values and fit the glm
      
      cluster_IDs <- colnames(mat_to_melt)
      
      Variable <- melt(mat_to_melt)$value
      
      #### building the factor for the glm based on the user input on the function
      
      factor_glm <- paste0(melt(mat_to_melt)$Var2,
                           rep(factorVector, length(unique(melt(mat_to_melt)$Var2))))
      
      
      if (ClassComparison == FALSE) {
        
        cat("No class comparison made")
        
        #plot(1)
        
        } else if (ClassComparison == "GLM") {
          
          #### getting out the glm P-values for each cluster separately
          
          P_value <- c()
          
          for (i in 1:length(cluster_IDs)) {
            
            glm_indexes <- grep(cluster_IDs[i], factor_glm)
            
            if (is.null(confoundingFactor)) {
              
            } else {
              
              ### testing confounding factor
              
              GLM_simple <- glm(Variable[glm_indexes] ~ factor_glm[glm_indexes])
              
              extra_factor <- paste0(factor_glm[glm_indexes], confoundingFactor)
              
              GLM_extra <- glm(Variable[glm_indexes] ~ extra_factor)
              
              cat("Feature: ", names(list_out)[m], "\n")
              cat("Cluster No. ", cluster_IDs[i], "\n")
              cat("AIC: Simple GLM ",  GLM_simple[["aic"]],
                  " ",
                  "Extra GLM ",  GLM_extra[["aic"]], "\n", "\n")
              
            }
            
            ### final test (not customized yet to minimum AIC)
          
            P_value <- c(P_value, summary(glm(Variable[glm_indexes] ~ factor_glm[glm_indexes]))$coefficients[,"Pr(>|t|)"])
          
          }
        
          P_Values[[m]] <- P_value
        
          color_plot <- rep(colorVector, (length(P_value)/FactorNumber))
        
          #### assigning stars to P-values, ° = +inf:0.1, * = 0.1:0.05, ** 0.05:0.01, *** 0.01:0
        
          stars_out <- c()
        
          for (i in 1:length(P_value)) {
          
            if (P_value[i] < 0.1 & P_value[i] > 0.05) {
              runner = "*"
            } else if (P_value[i] < 0.05 & P_value[i] > 0.01) {
              runner = "**"
            } else if (P_value[i] < 0.01) {
              runner = "***"
            } else {
              runner = "°"
            }
          
            stars_out[i] <- runner 
          
          }
        
          stars_out[grep("Intercept", names(P_value))] <- "C"
          
          P_value[grep("Intercept", names(P_value))] <- NA
        
          #### plotting
        
          par(mar=c(10,2,2,2))
        
          boxplot(Variable ~ factor_glm,
                  ylim = c(0, (max(Variable) + mean(Variable))), 
                  main = names(list_out)[m],
                  las = 2,
                  col = color_plot)
        
          text(x = c(1:length(unique(factor_glm))),
               y = (max(Variable) + mean(Variable)),
               labels = round(P_value, 2))
        
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
        
      } else if (ClassComparison == "ANOVA") {
        
        if (length(unique(melt(mat_to_melt)$Var2)) < 2) {
          
          ## one factor mat
          
          data = cbind.data.frame(treatment = as.factor(factor_glm),
                                  value = as.numeric(Variable))
          
          model=lm(data$value ~ data$treatment)
          ANOVA=aov(model)
          
          P_value = summary(ANOVA)[[1]][["Pr(>F)"]]
          P_Values[[m]] <- as.double(na.omit(P_value))
          
          #### assigning stars to P-values, ° = +inf:0.1, * = 0.1:0.05, ** 0.05:0.01, *** 0.01:0
          
          stars_out <- c()
          
          for (i in 1:length(P_value)) {
            
            if (is.na(P_value[i])) {
              
            } else if (P_value[i] < 0.1 & P_value[i] > 0.05) {
              runner = "*"
            } else if (P_value[i] < 0.05 & P_value[i] > 0.01) {
              runner = "**"
            } else if (P_value[i] < 0.01) {
              runner = "***"
            } else {
              runner = "°"
            }
            
            stars_out[i] <- runner 
            
          }
          
          stars_out[1] <- "C"
          
          #### plotting
          
          color_plot <- rep(colorVector, (length(P_value)/FactorNumber))
          
          par(mar=c(10,2,2,2))
          
          boxplot(Variable ~ factor_glm,
                  ylim = c(0, (max(Variable) + mean(Variable))), 
                  main = names(list_out)[m],
                  las = 2,
                  col = color_plot,
                  xlab = "Factors")
          
          text(x = c(1:length(unique(factor_glm))),
               y = (max(Variable) + mean(Variable)),
               labels = c(NA, round(P_value, 2)))
          
          text(x = c(1:length(unique(factor_glm))),
               y = (max(Variable) + (mean(Variable))/2),
               labels = stars_out)
          
        } else {
          
          ## multiple factor mat
          
          data = cbind.data.frame(treatment = as.factor(factor_glm),
                                  value = as.numeric(Variable))
          
          model=lm(data$value ~ data$treatment)
          ANOVA=aov(model)
          
          # Tukey test to study each pair of treatment
          TUKEY <- TukeyHSD(x=ANOVA, 'data$treatment', conf.level=0.95)
          
          P_value = TUKEY[["data$treatment"]][,"p adj"]
          P_Values[[m]] <- P_value
          
          # Tukey test representation
          par(mar=c(5,2,2,2))
          plot(TUKEY , las=1 , col="brown")
          
          
          ##### I need to group the treatments that are not different from each other together
          generate_label_df <- function(TUKEY, var){
            
            ###### Extract labels and factor levels from Tukey post-hoc 
            Tukey.levels <- TUKEY[[var]][,4]
            Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
            
            ###### I need to put the labels in the same order as in the boxplot
            Tukey.labels$treatment=rownames(Tukey.labels)
            Tukey.labels = Tukey.labels[order(Tukey.labels$treatment),]
            return(Tukey.labels)
          }
          
          ##### Apply the function on my dataset
          LABELS = generate_label_df(TUKEY , "data$treatment")
          
          
          ##### A panel of colors to draw each group with the same color
          my_colors = c(rgb(143,199,74,maxColorValue = 255),
                        rgb(242,104,34,maxColorValue = 255),
                        rgb(111,145,202,maxColorValue = 255),
                        rgb(254,188,18,maxColorValue = 255),
                        rgb(74,132,54,maxColorValue = 255),
                        rgb(236,33,39,maxColorValue = 255),
                        rgb(165,103,40,maxColorValue = 255))
          
          ##### Draw the basic boxplot
          par(mar=c(5,2,2,2))
          a = boxplot(data$value ~ data$treatment,
                      ylim =  c(0, (max(Variable) + mean(Variable))),
                      col=my_colors[as.numeric(LABELS[,1])],
                      main = names(list_out)[m],
                      las = 2,
                      cex.axis = 1) 
          
          ##### I want to write the letter over each box. Over is how high I want to write it.
          over = 0.1*max(a$stats[nrow(a$stats),])
          
          ##### Add the labels
          text(c(1:nlevels(data$treatment)),
               a$stats[nrow(a$stats),] + over ,
               LABELS[,1]  ,
               col=my_colors[as.numeric(LABELS[,1])] )
          
        }
        
      } else {
        
        cat("Error: Input a correct Class Comparison method (ANOVA or GLM)")
        
      }
      
    }
    
  }
  
  dev.off()
  
  names(P_Values) <- names(ClassDiscoveryList)
  
  ### start here the volcano plots part
  
  
  return(P_Values)
  
}
