#' A function that allows class comparison in replicated KineticMSI data subsets
#'
#' This function allows KineticMSI users to use classical class comparison algorithms such as ANOVA + Tukey HSD or generalized linear models (GLMs) + FDR correction. The input for the function comes from the previous kClassDiscoveryMSI.R output. The general assumption of this step is that multi-factorial designs were split at the class discovery level to allow discovery of inner data structures that are factor specific, which are then herein compared. The function outputs graphics embedded in PDF files detailing the results from the class comparison algorithms. Additionally, a table with all the results is returned to the environment when the output is assigned to an object.
#' @param kDiscoveryFactor1 output from the previous KineticMSI function. Namely, kClassDiscoveryMSI.R. The object must be a list of matrices, one matrix for each molecular feature of interest measured across replicates in factor1 (i.e., rows in each matrix).
#' @param kDiscoveryFactor2 output from the previous KineticMSI function. Namely, kClassDiscoveryMSI.R. The object must be a list of matrices, one matrix for each molecular feature of interest measured across replicates in factor2 (i.e., rows in each matrix).
#' @param factor1 character vector that needs to define the treatments using the same nomenclature and order in the naming scheme as the row names in matrices within the input files from kDiscoveryFactor1 (follow the exemplary KineticMSI data for details).
#' @param factor2 character vector that needs to define the treatments using the same nomenclature and order in the naming scheme as the row names in matrices within the input files from kDiscoveryFactor2 (follow the exemplary KineticMSI data for details).
#' @param repNumber1 number of biological replicates in dataset 1
#' @param repNumber2 number of biological replicates in dataset 2
#' @param PDFname defaults to "test". Defines the name of the output PDF files that will be produce as outcome in the current working directory.
#' @param returnGLMplots defaults to TRUE. Defines whether GLM plots are returned.
#' @param patternGLMplot defaults to "Q values". Defines whether Q or P values are used to return significance on top of the boxplots within the PDF containing the GLM results.
#' @param ylabTukey defaults to NULL. Character vector that defines the text on y-axes from Tukey HSD plots.
#' @param xlabTukey defaults to NULL. Character vector that defines the text on x-axes from Tukey HSD plots.
#' @param ylabGLM defaults to NULL. Character vector that defines the text on y-axes from GLM plots.
#' @param xlabGLM defaults to NULL. Character vector that defines the text on x-axes from GLM plots.
#' @keywords MSI Replicates Class Comparison Subsets Tracer Dynamics
#' @export
#' @examples
#'
#' ...


kSubSetClassComparisonMSI <- function(kDiscoveryFactor1,
                                      kDiscoveryFactor2,
                                      factor1,
                                      factor2,
                                      repNumber1,
                                      repNumber2,
                                      PDFname,
                                      returnGLMplots = T,
                                      patternGLMplot = c("Q values", "P values"),
                                      ylabTukey = NULL,
                                      xlabTukey = NULL,
                                      ylabGLM = NULL,
                                      xlabGLM = NULL) {

  ## Main

  ### Data quality

  in_list1 = kDiscoveryFactor1[which(!is.na(names(kDiscoveryFactor1)))]

  in_list2 = kDiscoveryFactor2[which(!is.na(names(kDiscoveryFactor2)))]

  if(length(in_list1) != length(in_list2)) {

    cat("\n")
    cat("Provided lists have different length... subsetting to common entities \n")
    cat("\n")

    minList <- intersect(names(in_list1), names(in_list2))

    in_list1_2 <- in_list1[minList]

    in_list2_2 <- in_list2[minList]

    if (unique(names(in_list1_2) == names(in_list1_2)) == T) {

      cat("Succesfully removed empty and/or missing entities... total No. ",
          length(in_list1_2))

    } else {

      stop("Error: not all entities are in both lists")

    }

  } else {

    in_list1_2 <- in_list1

    in_list2_2 <- in_list2

    cat("\n")
    cat("Provided lists have the same length...\n")
    cat("\n")

  }


  ### loop to access all the cluster comparisons

  iter_object <- names(in_list1_2)

  out_list <- list()

  names <- c()

  pdf(paste0(PDFname,  "TukeyHSD.pdf"))
  par(mar = c(5, 10, 3, 2))
  par(mfrow = c(2,1))

  for (i in 1:length(iter_object)) {

    ### preparing all objects to be ready for the pairwise combinations

    runner_entity <- iter_object[i]

    FactorReps <- c(rep(factor1, length(in_list1_2[[runner_entity]][["Values"]])),
                    rep(factor2, length(in_list2_2[[runner_entity]][["Values"]])))

    if (length(FactorReps) != 0) {

      ### only joined matrices are allowed beyond NULL objects in both input lists

      boundMats <- c(in_list1_2[[runner_entity]][["Values"]], in_list2_2[[runner_entity]][["Values"]])

      row_means_runner <- list()

      Tukey_joined_vec_data <- c()

      Tukey_joined_vec_factor <- c()

      for (j in 1:length(boundMats)) {

        row_means_runner[[j]] <- rowMeans(boundMats[[j]])

        Tukey_joined_vec_data <- c(Tukey_joined_vec_data, rowMeans(boundMats[[j]]))

        Tukey_joined_vec_factor <- c(Tukey_joined_vec_factor,
                                     rep(paste0(FactorReps[j],"_", names(boundMats)[j]),
                                         dim(boundMats[[j]])[1]))

      }

      if (length(unique(Tukey_joined_vec_factor)) > 1) {

        ### only comparisons are allowed beyond one factor

        ### Tukey-HSD

        RandoDiStats::TukeyCustomized(variable = Tukey_joined_vec_data,
                                      factor = as.factor(Tukey_joined_vec_factor),
                                      MainTitle = runner_entity,
                                      returnObject = "MeanComparisons",
                                      ylabTukeys = ylabTukey,
                                      xlabTukeys = xlabTukey)

        ### creating combination matrices for each pairwise comparison

        #### y = 0.5x2 - 0.5x + 6E-14

        combiNo <- 0.5*(length(boundMats)^2) - (0.5*length(boundMats)) + 6E-14

        combiMat <- combn(length(boundMats), m = 2)

        for (j in 1:combiNo) {

          runner_combi <- row_means_runner[combiMat[,j]]

          combi_enrichments <- c()

          names_rows <- c()

          for (m in 1:length(runner_combi)) {

            combi_enrichments <- c(combi_enrichments, runner_combi[[m]])

          }

          combi_enrichments = as.matrix(combi_enrichments)

          names <- c(names, paste0(runner_entity, "_",
                                   names(boundMats)[combiMat[,j]], "_",
                                   FactorReps[combiMat[,j]],
                                   collapse = ";"))

          out_list <- c(out_list, list(combi_enrichments))

        }
      }
    }
  }

  names(out_list) <- names

  dev.off()


  ### stats

  input_list <- out_list

  normality_test <- c()

  p_value <- c()

  E_size <- c()

  MC_vector <- c()

  test_mat <- matrix(NA,
                     nrow = dim(input_list[[1]])[1],
                     ncol = length(input_list))

  factorVectorList <- list()

  pdf(paste0(PDFname, ".pdf"))
  par(mar = c(5, 5, 4, 2))

  for (i in 1:length(out_list)) {

    #### defining factorVector

    factorVector <- as.factor(c(rep(strsplit(names(input_list)[i], ";")[[1]][1],
                                    length(input_list[[i]])/2),
                                rep(strsplit(names(input_list)[i], ";")[[1]][2],
                                    length(input_list[[i]])/2)))

    factorVectorList[[i]] <- factorVector

    rownames(input_list[[i]]) <- factorVector

    #### performing tests

    runner <- rowMeans(input_list[[i]])

    test_mat[,i] <- runner

    ## defining factors

    runner_means <- rowMeans(input_list[[i]])

    x_ks <- reshape2::melt(input_list[[i]][grep(unique(factorVector)[1],
                                                names(rowMeans(input_list[[i]]))),])$value

    y_ks <- reshape2::melt(input_list[[i]][grep(unique(factorVector)[2],
                                                names(rowMeans(input_list[[i]]))),])$value

    ## mean comparison (Log2FC)

    MC_vector[i] <- log2(mean(x_ks) / mean(y_ks))

    ## KS test

    normality_test[i] <- paste0(ks.test(x = x_ks, y = y_ks)[["method"]],
                                ": p = ",
                                ks.test(x = x_ks, y = y_ks)[["p.value"]])

    p_value[i] <- ks.test(x = x_ks, y = y_ks)[["p.value"]]

    ## effect size test

    E_size[i] <- as.numeric(as.character(effsize::cohen.d(d = x_ks, f = y_ks))[3])

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
           legend =  c(as.character(unique(factorVector)[2]),
                       as.character(unique(factorVector)[1])),
           col = c("blue", "red"), lty=1, cex=0.8)

    ## remove to force replacement in next iteration

    rm(x_ks)

    rm(y_ks)

  }

  dev.off()

  ## object names

  colnames(test_mat) <- names(input_list)
  rownames(test_mat) <- NULL

  names(normality_test) <- names(input_list)

  ## stat test (remove NAs from the output objects if present)

  test_mat1 <- test_mat[,!is.na(colMeans(test_mat))]

  if(dim(test_mat1)[2] != dim(test_mat)[2]) {

    test_mat <- test_mat1

    normality_test <- normality_test[!is.na(colMeans(test_mat))]

    p_value <- p_value[!is.na(colMeans(test_mat))]

    E_size <- E_size[!is.na(colMeans(test_mat))]

    MC_vector <- MC_vector[!is.na(colMeans(test_mat))]

    factorVectorList <- factorVectorList[!is.na(colMeans(test_mat))]

  }

  ## use factorVectorList to customize the output of the tests

  Pvalues_test <- RandoDiStats::OmicsUnivariateStats(class_comparison_mat = test_mat,
                                                     Factor1 = as.factor(c(rep("F1", repNumber1),
                                                                           rep("F2_GLM", repNumber2))),
                                                     Contrast = F,
                                                     ReturnTukeyPlots = T,
                                                     TukeyReturns = "MeanComparisons",
                                                     returnObject = "OmicsTests")

  colnames(Pvalues_test)[length(colnames(Pvalues_test))] <- "Factor1F2_TukeyHSD_Padj"

  ## out_object

  p_Adjustment <- cbind(normality_test,
                        KS_pAdjusted = p.adjust(p_value, method = "BH"),
                        Cohensd = E_size,
                        Log2FC = MC_vector,
                        t(test_mat),
                        t(as.data.frame(factorVectorList)),
                        Pvalues_test)

  if(returnGLMplots == T){

    pdf(paste0(PDFname, "GLM.pdf"))

    for (i in 1:dim(p_Adjustment)[1]) {

      GLMplotCustomizedMSI(Variable = as.numeric(p_Adjustment[i,5:(5+(dim(test_mat)[1])-1)]),
                           factorVector = as.factor(p_Adjustment[i,(5+(dim(test_mat)[1])):((5+(dim(test_mat)[1]))+(dim(test_mat)[1])-1)]),
                           Pvalues = as.numeric(p_Adjustment[i,grep(patternGLMplot, colnames(p_Adjustment))]),
						               MainTitle = rownames(p_Adjustment)[i],
						               ylabGLMs = ylabGLM, xlabGLMs = xlabGLM)

    }

    dev.off()

  }

  return(p_Adjustment)

}
