#' A function to support generalized linear model based class comparison in replicated KineticMSI datasets
#'
#' This function allows KineticMSI class comparison functions to create plots based on generalized linear models.
#' @param Variable Numeric vector where the means are inherited from to perform comparisons according to the factor vector.
#' @param factorVector factor type of vector that needs to define the treatments using the same nomenclature and order in the naming scheme as the names in the input variable.
#' @param Pvalues Vector of Q or P values that are used to return significance on top of the boxplots within the PDF containing the GLM results. Inherits from class comparison functions.
#' @param MainTitle Allows defining the title in the plots
#' @param ylabGLMs defaults to NULL. Character vector that defines the text on y-axes from GLM plots. Inherits from class comparison functions.
#' @param xlabGLMs defaults to NULL. Character vector that defines the text on x-axes from GLM plots. Inherits from class comparison functions.
#' @keywords MSI Replicates Class Comparison Tracer Dynamics
#' @export
#' @examples
#' 
#' ...


GLMplotCustomizedMSI <- function(Variable,
                                 factorVector,
                                 Pvalues,
                                 MainTitle = "",
                                 ylabGLMs = NULL,
                                 xlabGLMs = NULL){
  
  colorVector = rand_color(length(levels(factorVector)))
  
  FactorNumber <- length(levels(factorVector))
  
  
  #### Using the glm P-values
  
  color_plot <- rep(colorVector, (length(Pvalues)/FactorNumber))
  
  #### assigning stars to P-values, ° = +inf:0.1, * = 0.1:0.05, ** 0.05:0.01, *** 0.01:0
  
  stars_out <- c()
  
  for (i in 1:length(Pvalues)) {
    
    if (Pvalues[i] < 0.1 & Pvalues[i] > 0.05) {
      runner = "*"
    } else if (Pvalues[i] < 0.05 & Pvalues[i] > 0.01) {
      runner = "**"
    } else if (Pvalues[i] < 0.01) {
      runner = "***"
    } else {
      runner = "°"
    }
    
    stars_out[i] <- runner 
    
  }
  
  stars_out[1] <- "C"
  
  #### plotting
  
  par(mar=c(10,5,3,2))
  
  boxplot(Variable ~ factorVector,
          ylim = c(0, (max(Variable) + mean(Variable))), 
          main = MainTitle,
          las = 2,
          col = color_plot,
          xlab = xlabGLMs,
          ylab = ylabGLMs)
  
  text(x = c(1:length(unique(factorVector))),
       y = (max(Variable) + mean(Variable)),
       labels = round(Pvalues, 2))
  
  text(x = c(1:length(unique(factorVector))),
       y = (max(Variable) + (mean(Variable))/2),
       labels = stars_out)
  
  
  #### Add data points (https://www.r-graph-gallery.com/96-boxplot-with-jitter.html)
  
  data <- data.frame(names = factorVector, value = Variable)
  
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

