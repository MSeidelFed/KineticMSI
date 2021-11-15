#' A function to select an isotope proxy that best reflects the changes seen in steady state pools measured without stable isotope assisted mass spectrometry
#'
#' This function allows to compare the steady state pools from molecular features that were measured through both stable assisted and conventional mass spectrometry. By contrasting the results from different isotope combinations (e.g., a0 to an versus a0 to a1) users can easily define how many isotopes reflect the actual pool changes of the target molecular features in their biological systems. The input files must have the same treatment rows while molecular feature may vary. The algorithm will select the features in common for the calculations.
#' @param LabelledFileDir directory where the enrichment-calculated feature steady state pool input file is stored.
#' @param TreatmentFileDir directory where a treatments csv file defining rows in the both labelled and non-labelled steady state pool files is located. Must contain a "Batch" column defining the batches to be corrected across files.
#' @param NLSteadyStatePoolsDir directory where the non-labelled feature steady state pool input file is stored.
#' @param BatchCorr defaults to TRUE. When TRUE, batch correction is performed.
#' @param TreatmentIntoMeans defaults to TRUE. When true, replicates are turned into means for the final diagnostic heatmap.
#' @param Factor factor vector that defines the treatment of each row in the steady state input files. Follows the same order as the treatments file.
#' @param Duplicate defaults to FALSE. Allows to correct batches of two replicates by duplicating the entries to enhance batch testing.
#' @keywords Isotope Incorporation Proxy Selection
#' @export
#' @examples
#'
#' ...


ProxySelection <- function(LabelledFileDir,
                           TreatmentFileDir,
                           NLSteadyStatePoolsDir,
                           BatchCorr = TRUE,
                           TreatmentIntoMeans = TRUE,
                           Factor,
                           Duplicate = FALSE) {

  ## functions
  ### batch correction function

  BatchCorrection <- function(array_dir,
                              Treatments_dir,
                              Duplicate = Duplicate){

    #### building array

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

    #### preparing the data in the correct format

    array_intensities <- array[4:dim(array)[1],]
    colnames(array_intensities) <- colnames(array)

    Numeric_array <-apply(X = array_intensities, MARGIN = 2, FUN = as.numeric)

    Numeric_array[is.na(Numeric_array)] <- 0
    rownames(Numeric_array) <- rownames(array_intensities)
    colnames(Numeric_array) <- as.character(colnames(array_intensities))

    Numeric_array <- Numeric_array[as.numeric(rowMeans(Numeric_array)) != 0,]

    batch <- array["Batch",]


    #### main

    KineticMSI::ClassDistribution(inMat = Numeric_array,
                                  Treatments = batch)

    batch_necessity <- readline(prompt="Do you want to proceed with Batch correction? (Y/N)")


    if (batch_necessity == "Y") {

      NC_Normalized <- sva::ComBat(dat =  Numeric_array, batch = batch,
                                   mod = NULL,  par.prior = T, prior.plots = T)

      colnames(NC_Normalized) <- colnames(Numeric_array)
      rownames(NC_Normalized) <- rownames(Numeric_array)

      ClassDistribution(inMat = NC_Normalized,
                         Treatments = batch)

      return(NC_Normalized)

    } else {

      stop("Exit")

    }

  }


  ## Main

  if (BatchCorr == T) {

    DF_L <- BatchCorrection(array_dir = LabelledFileDir,
                            Treatments_dir = TreatmentFileDir,
                            Duplicate = Duplicate)

  } else {

    DF_L <- t(read.csv(file = LabelledFileDir, header = T, row.names = 1))

  }

  ### Importing non-labelled pools

  DF_NL <- t(read.csv(file = NLSteadyStatePoolsDir,
                      header = T, row.names = 1))

  TreatmentFile <- read.csv(file = TreatmentFileDir, header = T, row.names = 1)

  ### LvsNL pools
  #### common features between both sets

  common_features <- intersect(rownames(DF_L), rownames(DF_NL))

  DF_NL <- DF_NL[common_features,]

  DF_L <- DF_L[common_features,]

  #### getting ratios

  ratio_L_NL <- DF_L / DF_NL

  ratio_L_NL[is.na(ratio_L_NL)] <- 0
  ratio_L_NL[ratio_L_NL == Inf] <- 0
  ratio_L_NL[ratio_L_NL == -Inf] <- 0

  ### aggregating reps into means

  Treatment = Factor

  if (TreatmentIntoMeans == T) {

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

  ### building comparative heatmap

  ha = ComplexHeatmap::HeatmapAnnotation(df = data.frame(type = type))

  h <- ComplexHeatmap::Heatmap(matrix = ratio_L_NL,
                                name = "Intensities",
                                km = 1,
                                col = circlize::colorRamp2(c(0.5, 0.75, 1, 1.25, 1.5),
                                                           c("yellow", "white", "white", "white", "grey")),
                                show_column_names = F,
                                show_row_names = T,
                                cluster_columns = F,
                                cluster_rows = T,
                                top_annotation = ha)
  print(h)

  return(ratio_L_NL[ComplexHeatmap::row_order(h),])

}
