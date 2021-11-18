#' A function that allows summarizing the statistical output from KineticMSI full workflow
#'
#' This function allows KineticMSI users to to summarize KineticMSI output in two complementary steps. First the function offers and option to draw volcano plots using the different previously calculated metrics as axes. For instance, P values, Q values, Tukey HSD Padj, Kolmogorov-Smirnov P values in the y-axis, and Log fold changes or Cohen's D statistic in the x-axis. from parameter two until 25, the function call allows to fine tune the volcano plot that comes out from the calculations, extending the customization for enhanced biological understanding. Second, the function allows to perform a pathway enrichment analyses given the appropriate file format with the pathway information (see exemplary file in ...). Parameters 26 to 41 allow to fine tune the details of the pathway enrichment analysis output, which is performed thorugh a Fisher exact test. The test inherits the full categories from the input files and hence supports its legitimacy on the biological accuracy of the files that are provided by the user. Finally, the function is able to return a list of plots with the desired outcome as well as the Fisher enriched categories, if present, in a .csv file named "FisherResults.csv".
#' @param kComparisonOutput output from "kClassComparisonMSI.R" or "kSubSetClassComparisonMSI.R"
#' @param Abscissa name of the input column in "kComparisonOutput" to be used as the volcano abscissa.
#' @param Ordinate name of the input column in "kComparisonOutput" to be used as the volcano ordinate.
#' @param SigBndrie significance numerical boundary below which data points are labelled and significance vertical lines drawn in the volcano.
#' @param xBndrie numerical boundary above which data points are labelled and significance horizontal lines drawn in the volcano.
#' @param label parameter that allows to use as data point labels either the molecular entity name "EntityName" or the row name from "kComparisonOutput" when set to "Rownames".
#' @param plotTitle character vector inherited as a title for the volcano plots.
#' @param factor1 character that must coincide with factor 1 used at the class comparison functions.
#' @param factor2 character that must coincide with factor 2 used at the class comparison functions.
#' @param returnPlotsPNGs defaults to TRUE, builds .png files for every plot.
#' @param AbscissasName character to name the abscissa axis
#' @param OrdinatesName character to name the ordinate axis
#' @param width defaults to 44.45. graphical parameter for the saved .png volcano plots.
#' @param height defaults to 27.78. graphical parameter for the saved .png volcano plots.
#' @param AxesIndexSize defaults to 10 points.
#' @param AxesTitleSize defaults to 10 points.
#' @param LegendFontSize defaults to 10 points.
#' @param LegendKeySize defaults to 5 points.
#' @param LegendTitleSize defaults to 10 points.
#' @param LabelSize defaults to 5 points.
#' @param ColorAllDots defaults to FALSE. When TRUE colors all the plots with the category color, including those below the "SigBndrie" and "xBndrie".
#' @param DotsSize defaults to 5 points.
#' @param LowerLimitVolcano defaults to NULL. Allows to customize the lower limit of the volcano x-axis.
#' @param UpperLimitVolcano defaults to NULL. Allows to customize the upper limit of the volcano x-axis.
#' @param VolcanoLegendPosition defaults to "right". Allows to customize the position of the legend in the volcano plots.
#' @param FisherTermsDir directory where the file with the functional categories that will be used for the Fisher exact test of pathway enrichment are stored. The file must be a comma separated file (.csv). the first column contains the molecular entities and the subsequent columns the functional categories; each in a separate column, with the first row being the column name. For details see the exemplary file named xx.
#' @param returnSigFisher defaults to FALSE. When TRUE returns only the significant categories.
#' @param returnQvaluesFisher defaults to FALSE. When TRUE significances are only taken from adjusted P values, i.e., Q values.
#' @param h1Fisher defaults to "two.sided". Can be also "greater" or "less" and refers to whether functional categories are judged to be underrepresented and over-represented ("two.sided") as compared to all the significances, or whether only under-representation or over-representation is judged by the Fisher exact test.
#' @param ColorCategory name of the column in "FisherTermsDir" file from where the colors are taken for all the dots and categories in both plots.
#' @param PathwayColorsDir directory where the file defining colors for molecular entities and functional categories is located. The file must be tab separated (.txt). The first column contains the names of the molecular entities, the second column the functional category that will be used to implement colors and the third column R-hexadecimal color identifiers. For details see the exemplary file named xx.
#' @param BarNumbersPlot defaults to "pSigEntInClassToTotEntInClass" (proportion of significant entities in class as compared to total entities in class). Refers to the magnitude displayed in in the Fisher bar plots that show the pathway enrichment analysis results. Can also be "pEntInClassToTotEnt" "SigEntInCat".
#' @param FisherCatNumSize defaults to 3.5 points.
#' @param FisherCatAxisFontSize defaults to 10 points.
#' @param FisherAxesTitleSize defaults to 10 points.
#' @param FisherLegendFontSize defaults to 10 points.
#' @param FisherLegendTitleSize defaults to 10 points.
#' @param FisherLegendKeySize defaults to 10 points.
#' @param LowerLimitFisher defaults to NULL. Allows to customize the lower limit of the volcano x-axis.
#' @param UpperLimitFisher defaults to NULL. Allows to customize the lower limit of the volcano x-axis.
#' @param FisherLegendPosition defaults to "right". Allows to customize the position of the legend in the Fisher plots.
#' @param returnObject defaults to "VolcanoPlots", which returns a list with the ggplot objects. Can also be "FisherPlots" or "matrices".
#' @keywords MSI Replicates Summary Results Subsets Tracer Dynamics
#' @export
#' @examples
#'
#' ...


kSummaryMSI <- function(kComparisonOutput,
                        Abscissa = c("Cohensd", "Log2FC"),
                        Ordinate = c("KS_pAdjusted", "Factor1F2_GLM_P values", "Factor1F2_GLM_Q values", "Factor1F2_TukeyHSD_Padj"),
                        SigBndrie = 0.01,
                        xBndrie = 2,
                        label = c("Rownames", "EntityName"),
                        plotTitle = "test",
                        factor1,
                        factor2,
                        returnPlotsPNGs = TRUE,
                        AbscissasName = "Abscissas",
                        OrdinatesName = c("Ordinates", expression(-Log[10]~italic(P)[adj]~ ~value)),
                        width = 44.45,
                        height = 27.78,
                        AxesIndexSize = 10,
                        AxesTitleSize = 10,
                        LegendFontSize = 10,
                        LegendKeySize = 5,
                        LegendTitleSize = 10,
                        LabelSize = 5,
                        ColorAllDots = FALSE,
                        DotsSize = 5,
                        LowerLimitVolcano = NULL,
                        UpperLimitVolcano = NULL,
                        VolcanoLegendPosition = c("right", "bottom", "none", "top", "left", "right"),
                        FisherTermsDir = paste0(system.file("extdata", package = "KineticMSI"), "/PA_Terms.csv"),
                        returnSigFisher = FALSE,
                        returnQvaluesFisher = FALSE,
                        h1Fisher = c("two.sided", "greater", "less"),
                        ColorCategory,
                        PathwayColorsDir = c(NULL, paste0(system.file("extdata", package = "KineticMSI"), "/PA_TermsColors.txt")),
                        BarNumbersPlot = c("pSigEntInClassToTotEntInClass", "pEntInClassToTotEnt", "SigEntInCat"),
                        FisherCatNumSize = 3.5,
                        FisherCatAxisFontSize = 10,
                        FisherAxesTitleSize = 10,
                        FisherLegendFontSize = 10,
                        FisherLegendTitleSize = 10,
                        FisherLegendKeySize = 10,
                        LowerLimitFisher = NULL,
                        UpperLimitFisher = NULL,
                        FisherLegendPosition = c("right", "bottom", "none", "top", "left", "right"),
                        returnObject = c("VolcanoPlots", "FisherPlots", "matrices")) {


  ## functions needed

  FisherExactMSI <- function(FisherDir,
                             FisherTermsDir,
                             returnSig = F,
                             returnQvalues = F,
                             alternative = c("two.sided", "greater", "less"),
                             BarNumbers = c("pSigEntInClassToTotEntInClass", "pEntInClassToTotEnt", "SigEntInCat"),
                             FisherCatNumSizes = 10,
                             FisherCatAxisFontSizes = 10,
                             FisherAxesTitleSizes = 10,
                             FisherLegendFontSizes = 10,
                             FisherLegendTitleSizes = 10,
                             FisherLegendKeySizes = 10,
                             LowerLimitsFisher = NULL,
                             UpperLimitsFisher = NULL){


    in_Fisher <- FisherDir

    Fisher_terms <- read.table(FisherTermsDir, header = T, sep = ",", row.names = 1)

    Fisher_terms[Fisher_terms == ""] <- NA

    listoflists <- list()

    for (i in 1:dim(Fisher_terms)[2]) {

      Enrichment_Term <- colnames(Fisher_terms)[i]

      Categories <- na.omit(unique(as.character(Fisher_terms[,Enrichment_Term])))

      list_mats <- list()

      names_mats <- c()

      for (j in 1:length(Categories)) {

        cat("...\n")
        cat("working on category: ", Categories[j], "\n")
        cat("...\n")

        # a is distributed as an hypergeometric distribution with
        # a+c draws from a population with
        # a+b successes and
        # c+d failures

        # a  |  b
        # c  |  d

        c <- unique(rownames(Fisher_terms)[grep(Categories[j], Fisher_terms[,Enrichment_Term])]) ## entities in category

        b <- as.character(in_Fisher)  ## entities significantly changed in total

        d <- unique(rownames(Fisher_terms)) ## total entitites

        sign_Categ <- intersect(c, b)

        if (length(sign_Categ) == 0){

          out_a <- 0

          ### print variables
          SigEntInCat <- out_a
          pSigEntInClassToTotEntInClass <- out_a/length(c)
          pEntInClassToTotEnt <- length(c)/length(d)

        } else {

          out_a <- c()

          for (k in 1:length(sign_Categ)) {

            out_a <-c(out_a,
                      grep(sign_Categ[k], d, value = T))   ## entities significantly changed in category

          }

          ### print variables
          SigEntInCat <- length(out_a)
          pSigEntInClassToTotEntInClass <- length(out_a)/length(c)
          pEntInClassToTotEnt <- length(c)/length(d)

        }

        ### contingency mats & proportions

        d = setdiff(d, b)

        if (length(sign_Categ) == 0){

          ContingencyMat <- as.matrix(cbind(c(out_a, length(c)), c(length(b), length(d))))

        } else {

          c = setdiff(c, out_a)

          ContingencyMat <- as.matrix(cbind(c(length(out_a), length(c)), c(length(b), length(d))))

        }


        list_mats[[j]] <- ContingencyMat

        names_mats[j] <- paste0(Enrichment_Term,";", Categories[j], ";",
                                scales::percent(round(pSigEntInClassToTotEntInClass, 2)),
                                ";", scales::percent(round(pEntInClassToTotEnt, 2)),
                                ";", SigEntInCat)

      }

      names(list_mats) <- names_mats

      listoflists <- c(listoflists, list_mats)

    }

    P_values <- c()

    for (i in 1:length(listoflists)) {

      P_values[i] <- fisher.test(listoflists[[i]], alternative = alternative)[["p.value"]]

    }

    Q_values <- p.adjust(P_values, "fdr")

    ## returning only significant or all results

    if (returnSig == T) {

      if (returnQvalues == T & length(unique(Q_values < 0.05)) == 1) {

        cat("...\n")
        cat("No significant Q values found....")
        cat("...\n")
        stop("There are no significant categories (Qvalues) in your dataset, try executing with returnQvalues = FALSE")

      } else if (returnQvalues == T & length(unique(Q_values < 0.05)) > 1) {

        Enriched_categories <- paste0(names(listoflists)[which(Q_values < 0.05)], ";" , Q_values[which(Q_values < 0.05)])

        expressionXaxis <- expression(-Log[10]~italic(P)[adj]~ ~value)

      } else if (returnQvalues == F & length(unique(P_values < 0.05)) == 1) {

        cat("...\n")
        cat("No significant P values found....")
        cat("...\n")
        stop("There are no significant categories (Pvalues) in your dataset, try executing with returnSig = FALSE")

      } else if (returnQvalues == F & length(unique(P_values < 0.05)) > 1) {

        Enriched_categories <- paste0(names(listoflists)[which(P_values < 0.05)], ";" , P_values[which(P_values < 0.05)])

        expressionXaxis <- expression(-Log[10]~italic(P)~value)

      }

    } else if (returnSig == F) {

      Enriched_categories <- paste0(names(listoflists), ";" , P_values)

      expressionXaxis <- expression(-Log[10]~italic(P)~value)

    }

    ## building the nice output; table and barplot

    in_object_test <- t(as.data.frame(strsplit(Enriched_categories, ";")))

    rownames(in_object_test) <- NULL

    colnames(in_object_test) <- c("Column_ID",
                                  "Assigned_function",
                                  "Functional_Categories",
                                  "pSigEntInClassToTotEntInClass",
                                  "pEntInClassToTotEnt",
                                  "SigEntInCat",
                                  "Fisher_P_value")

    in_object_test[in_object_test == "Inf"] <- "NA"

    write.csv(in_object_test, "FisherResults.csv", row.names = F)

    ### then use ggplot to build the barplot output

    in_object <- readr::read_csv("FisherResults.csv")

    ## set the levels in order we want
    ### reordered object

    in_object <- in_object[order(in_object$Functional_Categories),]

    if (BarNumbers == "pSigEntInClassToTotEntInClass"){

      labels <- in_object$pSigEntInClassToTotEntInClass

    } else if (BarNumbers == "pEntInClassToTotEnt") {

      labels <- in_object$pEntInClassToTotEnt

    } else if (BarNumbers == "SigEntInCat") {

      labels <- in_object$SigEntInCat

    } else {

      stop("ERROR: Input wrong Bar Number category....")

    }

    ### automatic axes if not defined by user

    if(is.null(UpperLimitsFisher)) {

      UpperLimitsFisher <- round(max(-log10(in_object$Fisher_P_value)) + 1, 0)

    }

    if(is.null(LowerLimitsFisher)) {

      LowerLimitsFisher <- 0

    }

    ### plot

    simplePlot <- ggplot2::ggplot(in_object, ggplot2::aes(x = -log10(Fisher_P_value),
                                                          y = forcats::fct_reorder((Assigned_function),
                                                                                    as.numeric(as.factor(Functional_Categories))),
                                                           fill = Functional_Categories)) +
                  ggplot2::geom_col(width=0.9) +
                  ggplot2::labs(x = expressionXaxis, y = NULL) +
                  ggplot2::geom_vline(xintercept = -log10(0.05), color='red', linetype = "solid")+
                  ggplot2::geom_text(ggplot2::aes(label = labels), hjust = -0.3, color="black", size = FisherCatNumSizes)+
                  ggplot2::scale_fill_viridis_d(option = "inferno")+
                  ggplot2::theme_classic()+
                  ggplot2::xlim(LowerLimitsFisher, UpperLimitsFisher)+
                  ggplot2::labs(fill = "Functional Categories")+
                  ggplot2::theme(axis.text.y = ggplot2::element_text(size = FisherCatAxisFontSizes),
                                 axis.title = ggplot2::element_text(size= FisherAxesTitleSizes,face="bold"),
                                 legend.text = ggplot2::element_text(size= FisherLegendFontSizes),
                                 legend.title = ggplot2::element_text(size= FisherLegendTitleSizes,face="bold"))+
                  ggplot2::theme(legend.position=FisherLegendPosition)+
                  ggplot2::guides(color = guide_legend(override.aes = list(size= FisherLegendKeySizes)))


    return(simplePlot)

  }


  ## main

  ### fisher terms file for plots

  Fisher_terms <- read.table(FisherTermsDir, header = T, sep = ",", row.names = 1)

  ### main continues

  in_mat <- kComparisonOutput

  in_mat[(as.numeric(in_mat[,Ordinate]) == 0), Ordinate] <- 0.000000000000000000001

  ## ??? why NAs to signifificant stuff --> in_mat[is.na(in_mat[,Ordinate]),Ordinate] <- 0.000000000000000000001

  replaceable <- c(which(as.numeric(in_mat[,Abscissa]) == "-Inf"),
                   which(as.numeric(in_mat[,Abscissa]) == "Inf"),
                   which(is.na(as.numeric(in_mat[,Abscissa]))))

  in_mat[replaceable, Abscissa] <- 0

  ## subsetting the entire dataset into common factorial objects

  Fac1Rows <- grep(factor1, rownames(in_mat))

  Fac2Rows <- grep(factor2, rownames(in_mat))


  if(length(Fac1Rows) == 0 & length(Fac2Rows) == 0) {

    cat("\n", "... Working with the population means..." , "\n")

    info_carried <- 1

    rownames_vec <- paste0(rep("WT-HD", dim(in_mat)[1]),"_", rownames(in_mat))

    rownames(in_mat) <- rownames_vec

    Fac1Rows <- grep(factor1, rownames(in_mat))

    Fac2Rows <- grep(factor2, rownames(in_mat))

  } else {

    info_carried <- 0

  }

  FactorNames <- c(paste0(factor1, "_", factor2), factor1, factor2)

  list_SubSets <- list(F1F2 = in_mat[intersect(Fac1Rows, Fac2Rows),],
                       F1 = in_mat[setdiff(Fac1Rows, Fac2Rows),],
                       F2 = in_mat[setdiff(Fac2Rows, Fac1Rows),])

  names(list_SubSets) <- FactorNames

  for (i in length(list_SubSets):1) {

    if(dim(list_SubSets[[i]])[1] == 0) {

      list_SubSets[[i]] = NULL

    }

  }


  list_HyperTest <- list()

  plot_list <- list()

  for (i in 1:length(list_SubSets)) {

    if(class(list_SubSets[[i]]) == "matrix") {

      in_mat2 <- list_SubSets[[i]]

      if (info_carried == 1) {

        namesVector <- apply(as.data.frame(strsplit(rownames(in_mat2), "_"))[2,], 2, as.character)

      } else {

        namesVector <- apply(as.data.frame(strsplit(rownames(in_mat2), "_"))[1,], 2, as.character)

      }

      rownames(in_mat2) <- namesVector

      Abscissas <- as.numeric(in_mat2[,Abscissa])

      Ordinates <- -log10(as.numeric(in_mat2[,Ordinate]))

      data <- cbind.data.frame(Abscissas, Ordinates)

      ## significant cases according to user input

      PositiveCases <- which(Abscissas > xBndrie & Ordinates > -log10(SigBndrie) | Abscissas < - xBndrie & Ordinates > -log10(SigBndrie))

      ## colors object

      Pathway_colors <- read.delim(PathwayColorsDir, row.names = 1)

      ## objects to fill for the plot

      Labels_corrected <- rep("", length(Abscissas))

      if(length(PositiveCases) > 0) {

        list_HyperTest[[i]] <- rownames(in_mat2)[PositiveCases]

        ## labels

        if (label == "EntityName") {

          #cat("Using entity names as label identifiers \n")

          Labels_corrected[PositiveCases] <- rownames(in_mat2)[PositiveCases]

        } else if (label == "Rownames") {

          #cat("Using entire rownames as label identifiers \n")

          Labels_corrected[PositiveCases] <- rownames(list_SubSets[[i]])[PositiveCases]

        } else {

          stop("ERROR: input either EntityName OR Rownames")

        }

        if (ColorAllDots == T) {

          AllEntities <- rownames(in_mat2)

          ## colors

          Colors_corrected <- as.character(Pathway_colors[AllEntities,2])

          Colors_corrected[is.na(Colors_corrected)] <- "#e5e5e5"

          colors2 <- factor(Colors_corrected)

          #### 20 is small circle and 19 big circle

          Colors_shape_corrected <- rep(20, length(Abscissas))

          Colors_shape_corrected[PositiveCases] <- 19

          #### fallback in case categories are missing is "-"

          Factors_corrected <- as.character(Pathway_colors[AllEntities,1])

          Factors_corrected[is.na(Factors_corrected)] <- "-"

          ### Using the colors generated or inputted from the table

          data[,3] <- as.factor(Factors_corrected)

          ### Manual legend

          Legend_input <- unique(cbind(Factors_corrected,Colors_corrected))

        } else {

          ColorsToCustomize <- as.character(Fisher_terms[Labels_corrected[Labels_corrected != ""], ColorCategory])

          ColorsToCustomize[which(is.na(ColorsToCustomize))] <- "Not-in-input-list"

          ## producing LegendFull object when colors are defined with an input table or randomly produced

          if (is.null(PathwayColorsDir)) {

            LegendFull <- cbind(unique(ColorsToCustomize),
                                rand_color(length(unique(ColorsToCustomize)), luminosity = "bright"))

          } else {

            uniqueCols <-unique(ColorsToCustomize)

            out_color <- c()

            for (j in 1:length(uniqueCols)) {

              test_runner <- as.character(unique(Pathway_colors[grep(strsplit(uniqueCols[j], ";")[[1]][1], Pathway_colors[,1]),2]))

              if (length(test_runner) > 1) {

                stop("ERROR: Redundant colors for the same category, check input table")

              } else if (length(test_runner) == 1) {

                out_color[j] <- test_runner

              } else {

                out_color[j] <- "black"

              }

            }


            LegendFull <- cbind(uniqueCols,
                                out_color)

          }

          ## producing out_colors for the plot from LegendFull

          colors2 <- LegendFull[,2]
          names(colors2) <- LegendFull[,1]

          out_colors <- matrix(data = NA, nrow = length(ColorsToCustomize), ncol = 2)

          for (j in 1:length(ColorsToCustomize)) {

            out_colors[j,] <- LegendFull[which(ColorsToCustomize[j] == LegendFull[,1]),]

          }

          ## colors

          ### colors will be inherited from categories and shapes from PositiveCases

          #### fallback in case categories are not inputed is grey "#e5e5e5"

          Colors_corrected <- rep("#e5e5e5", length(Abscissas))

          Colors_corrected[PositiveCases] <- out_colors[,2]#"black"

          #### 20 is small circle and 19 big circle

          Colors_shape_corrected <- rep(20, length(Abscissas))

          Colors_shape_corrected[PositiveCases] <- 19

          #### fallback in case categories are not inputed is "-"

          Factors_corrected <- rep("-", length(Abscissas))

          Factors_corrected[PositiveCases] <- out_colors[,1]

          ### Using the colors generated or inputed from the table

          data[,3] <- as.factor(Factors_corrected)

          ### Manual legend

          Legend_input <- unique(cbind(Factors_corrected,Colors_corrected))

        }

      } else {

        if (ColorAllDots == T) {

          AllEntities <- rownames(in_mat2)

          ## colors

          Colors_corrected <- as.character(Pathway_colors[AllEntities,2])

          Colors_corrected[is.na(Colors_corrected)] <- "#e5e5e5"

          #### 20 is small circle and 19 big circle

          Colors_shape_corrected <- rep(20, length(Abscissas))

          #### fallback in case categories are missing is "-"

          Factors_corrected <- as.character(Pathway_colors[AllEntities,1])

          Factors_corrected[is.na(Factors_corrected)] <- "-"

          ### Using the colors generated or inputted from the table

          data[,3] <- as.factor(Factors_corrected)

          ### Manual legend

          Legend_input <- unique(cbind(Factors_corrected,Colors_corrected))

        } else {

          Colors_corrected <- rep("#e5e5e5", length(Abscissas))

          Colors_shape_corrected <- rep(20, length(Abscissas))

          Legend_input <- cbind("-", "#e5e5e5")

          data[,3] <- rep("-", dim(data)[1])

        }

      }

      ## Volcano plot

      ### setting abscissa limits

      if(is.null(LowerLimitVolcano)) {

        LowerLimitVolcano <- (abs(min(Abscissas))+(abs(min(Abscissas))*0.1))*sign(min(Abscissas))

      }

      if(is.null(UpperLimitVolcano)) {

        UpperLimitVolcano <- (abs(max(Abscissas))+(abs(max(Abscissas))*0.1))*sign(max(Abscissas))

      }


      plot_list[[i]] <- ggplot2::ggplot(data, ggplot2::aes(x = Abscissas, y = Ordinates, color = factor(V3))) +
                        ggplot2::geom_point(shape = Colors_shape_corrected, size = DotsSize) +
                        ggplot2::scale_color_manual(values = Legend_input[,2], breaks = Legend_input[,1]) + #  , "-" = "#e5e5e5"
                        ggplot2::labs(color = ColorCategory, x = AbscissasName, y = OrdinatesName) +
                        ggplot2::theme_classic(base_size = 12) +
                        ggplot2::theme(axis.text = ggplot2::element_text(size=AxesIndexSize),
                                       axis.title = ggplot2::element_text(size=AxesTitleSize, face="bold"),
                                       legend.text = ggplot2::element_text(size= LegendFontSize),
                                       legend.title = ggplot2::element_text(size=LegendTitleSize, face="bold"),
                                       legend.position = VolcanoLegendPosition) +
                        ggplot2::geom_hline(yintercept = c(-log10(0.05), -log10(SigBndrie)),
                                            linetype=c("dashed", "solid"),
                                            color = c("orange", "red")) +
                        ggplot2::geom_vline(xintercept = c(-xBndrie,-0.2,0.2,xBndrie),
                                           linetype=c("solid","dashed","dashed", "solid"),
                                           color = c("red","orange","orange", "red")) +
                        ggplot2::ggtitle(names(list_SubSets)[i]) +
                        ggrepel::geom_text_repel(label = Labels_corrected, show.legend = F, size = LabelSize) +
                        ggplot2::xlim(LowerLimitVolcano, UpperLimitVolcano) +
                        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=LegendKeySize)))

    } else {

      list_HyperTest[[i]] <- NA

    }

    if (returnPlotsPNGs == T) {

      ggplot2::ggsave(plot_list[[i]], file=paste0("plot_", names(list_SubSets)[i], "_", plotTitle,".png"),
                      width = width, height = height, units = "cm", dpi=300)

    }

  }

  ## Fisher test if anything significant is there

  if (length(list_HyperTest) > 0) {

    names(list_HyperTest) <- names(list_SubSets)

    Fisher_plots <- list()

    for (i in 1:length(list_HyperTest)) {

      cat("\n")
      cat("Running Fisher exact test on: ", names(list_HyperTest)[i])
      cat("\n")

      runner <- FisherExactMSI(FisherDir = list_HyperTest[[i]],
                               FisherTermsDir = FisherTermsDir,
                               returnSig = returnSigFisher,
                               returnQvalues = returnQvaluesFisher,
                               alternative = h1Fisher,
                               BarNumbers = BarNumbersPlot,
                               FisherCatNumSizes = FisherCatNumSize,
                               FisherCatAxisFontSizes = FisherCatAxisFontSize,
                               FisherAxesTitleSizes = FisherAxesTitleSize,
                               FisherLegendFontSizes = FisherLegendFontSize,
                               FisherLegendTitleSizes = FisherLegendTitleSize,
                               FisherLegendKeySizes = FisherLegendKeySize,
                               LowerLimitsFisher = LowerLimitFisher,
                               UpperLimitsFisher = UpperLimitFisher)

      Fisher_plots[[i]] <- runner

      ggsave(runner, file=paste0("FisherBarplot_", names(list_HyperTest)[i], "_", plotTitle,".png"),
             width = 44.45, height = 27.78, units = "cm", dpi=300)

    }

  }

  if (returnObject == "VolcanoPlots") {

    return(plot_list)

  } else if (returnObject == "FisherPlots") {

    return(Fisher_plots)

  } else if (returnObject == "matrices") {

    return(list_SubSets)

  }

}
