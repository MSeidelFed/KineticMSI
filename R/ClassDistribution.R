#' A function to find out the data distribution from variables
#'
#' This function allows to calculate density shapes and distributions from variable values in a multi-treatment matrix.
#' @param inMat a matrix with multiple treatments and variables from which the density shape and distribution of the variable values needs to be known.
#' @param Treatments a factor vector defining treatments in the input matrix, length must match the dimensions from the matrix.
#' @param plots defaults to TRUE Allows to return the diagnostic plots (ggplot dependency) that specify the distributions.
#' @param returnTable defaults to FALSE. Allows to return the final melted data frame used for density calculations.
#' @param factorVector defaults to NULL. a factor vector that defines the identity of independent factors in the design, in order to color them differently in the plots.
#' @param PlotMain defaults to NULL. When defined as a character vector it will appear as the diagnostic plots main title.
#' @param ScalingFactorXaxis defaults to NULL. numeric that is typically defined in the parent function as a factor to scale the x-axis in density plots with heavy-tailed distributions.
#' @keywords Distribution ggplot stats
#' @export
#' @examples
#'
#' ...


### distribution ggplot function

ClassDistribution <- function(inMat,
                              Treatments,
                              plots = TRUE,
                              returnTable = FALSE,
                              factorVector = NULL,
                              PlotMain = NULL,
                              ScalingFactorXaxis = NULL) {

  mydata4.1 <- as.data.frame(t(inMat))
  colnames(mydata4.1) <- seq(1, dim(mydata4.1)[2])

  mydata4.2 <- cbind(Treatments, mydata4.1)
  colnames(mydata4.2) <- c("X1_1", seq(1,dim(mydata4.1)[2]))

  mydata4.3 <- reshape2::melt(mydata4.2, id.vars= seq(1), measure.vars= seq(2, (dim(mydata4.1)[2]+1)))

  mydata4.4 <- mydata4.3[order(mydata4.3$X1_1),]


  if(is.null(factorVector)){

    out_vec <- rep("black", length(mydata4.4$X1_1))

  } else {

    unique_factors <- unique(factorVector)

    out_vec <- c()

    for (i in 1:length(unique_factors)) {

      runner <- grep(unique(factorVector)[i], mydata4.4$X1_1)

      out_vec <- c(out_vec, rep(unique_factors[i], length(runner)))

    }

  }


  if(plots == T) {

    if(is.null(ScalingFactorXaxis)){

      ScalingFactorXaxis = 1

    }

    suppressMessages(print(ggplot2::ggplot(mydata4.4, ggplot2::aes(x = as.numeric(value), y = X1_1, color = out_vec)) +
                           ggplot2::xlim(summary(as.numeric(mydata4.4$value))["Min."] - summary(as.numeric(mydata4.4$value))["Mean"],
                                         summary(as.numeric(mydata4.4$value))["Mean"] + ScalingFactorXaxis) +
                           ggplot2::ggtitle(PlotMain) +
                           ggridges::geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
                           ggridges::theme_ridges(font_size = 10, grid = TRUE) +
                           ggplot2::theme(axis.title.y = ggplot2::element_blank())))


    suppressMessages(print(ggplot2::ggplot(mydata4.4, ggplot2::aes(x = as.numeric(value), fill = X1_1, color = out_vec)) +
                           ggplot2::geom_density(alpha = 0.2, position = "identity") +
                           ggplot2::xlim(summary(as.numeric(mydata4.4$value))["Min."] - summary(as.numeric(mydata4.4$value))["Mean"],
                                         summary(as.numeric(mydata4.4$value))["Mean"] * ScalingFactorXaxis) +
                           ggplot2::ggtitle(PlotMain)))

  }


  if(returnTable == T) {

    return(mydata4.4)

  }

}
