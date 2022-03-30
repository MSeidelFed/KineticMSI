#' A function to compare proportions within selected tracer enrichment proxies across pixels in replicated KineticMSI datasets
#'
#' This function allows KineticMSI users to compare the proportion of pixels that fall within specific ranges of selected tracer enrichment proxies. The user must define a threshold from which proportions are calculated and comparisons made using the parameter "ProportionLimit". Additionally the function provides the parameter "ProportionOperator" to define whether the comparisons are drawn in pixels "equal" to, "less" or "greater" than the predefined limit. The function returns to the R environment a matrix containing the values of the proportion comparison across molecular features (including statistical test outcomes) and an optional PDF with the graphical HeatMap representation of that matrix.
#' @param path parent directory were all enrichment files are contained, digs within recursive folders.
#' @param PatternEnrichment defaults to "MeanEnrichment". Defines a character vector used to grab input csv enrichment files that can be later subset.
#' @param SubSetRepsIntensities defaults to FALSE. Allows to subset the MSI file list found in path.
#' @param factorVector character vector that needs to define the treatments using the same nomenclature and naming scheme as for the input files (follow the exemplary KineticMSI data for details).
#' @param ProportionOperator allows to define whether the proportions to be compared are "equal", "less" or "greater" than a value defined in the function call parameter "ProportionLimit".
#' @param ProportionLimit defaults to 0. Allows to define the Threshold value against which all comparisons are made.
#' @param kmeans defaults to 5. Allows to define the number of clusters looked for in each feature dataset.
#' @param KmBoot defaults to 10. Inherits from ComplexHeatmap and defines the number of bootstrap iterations used to build the K-mean consensus.
#' @param ClustMethod defaults to "average". Inherits from ComplexHeatmap and supports all clustering methods described there.
#' @param returnProprotionsHeatmap defaults to TRUE. Allows users to decide whether a comparison HeatMap representation is returned as a PDF file named "ProportionsHeatmap.pdf".
#' @keywords MSI Replicates Proportion Tracer Dynamics
#' @export
#' @examples
#'
#' ...


kEnrichmentProportionsMSI <- function(path,
                                      PatternEnrichment = "MeanEnrichment",
                                      SubSetRepsIntensities = FALSE,
                                      factorVector,
                                      ProportionOperator = c("equal", "less", "greater"),
                                      ProportionLimit = 0,
                                      kmeans = 5,
                                      KmBoot = 10,
                                      ClustMethod = "average",
                                      returnProprotionsHeatmap = TRUE) {

  ## Main

  ### Importing all files

  reps <- list.files(path = path, pattern = PatternEnrichment,
                     all.files = FALSE, full.names = TRUE, recursive = TRUE,
                     ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)

  ### Sub-setting files

  if (SubSetRepsIntensities == T) {

    reps = SubSetInputFiles(reps = reps)

  }

  ### Producing files with shared features if necessary

  reps1 = ProduceNewFiles(reps = reps)

  if(!is.null(reps1)) {

    reps = reps1

  }

  ### getting the common features for all datasets to perform analyses

  feature_Nr <- list()

  feature_Nr_vec <- c()

  for (i in 1:length(reps)) {

    file_test <- read.csv(file = reps[i], header = T, row.names = 1)

    feature_Nr[[i]] <- unique(rownames(file_test))

    feature_Nr_vec <- c(feature_Nr_vec, rownames(file_test))

  }

  feature_Nr = Reduce(intersect, feature_Nr)


  ### Enrichment proportion function

  out_mat <- matrix(NA, nrow = length(feature_Nr), ncol = length(reps))

  for (j in 1:length(reps)) {

    test <- read.csv(file = reps[j], header = T, row.names = 1)

    out_vec <- c()

    for (i in 1:dim(test)[1]) {

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

  ## set to 1 NaN (come from no 0 in the matrices - signature region features)

  Padj_values[is.na(Padj_values)] <- 1

  P_values[is.na(P_values)] <- 1

  ## set to null names without sig padj values

  rownames(Padj_values)[which(Padj_values > 0.05)] <- "."

  rownames(P_values)[which(P_values > 0.05)] <- "."

  if (length(unique(rownames(Padj_values))) > 1) {

    P_heatmap = Padj_values

    Name_Heatmap <- "Q values"

  } else {

    P_heatmap = P_values

    Name_Heatmap <- "P values"

  }

  type = gsub("s\\d+_", "", unique(factorVector))

  ha = ComplexHeatmap::HeatmapAnnotation(df = data.frame(type = type))

  if (returnProprotionsHeatmap == T) {

    pdf("ProportionsHeatmap.pdf")

    Scal_Heatmap <- ComplexHeatmap::Heatmap(matrix = as.matrix(t(test_num)),
                                            name = "Intensities",
                                            km = kmeans, ## five K-means
                                            row_km_repeats = KmBoot, ## repeats to get consensus
                                            col = circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
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

      ComplexHeatmap::Heatmap(as.matrix(P_heatmap),
                              col = circlize::colorRamp2(c(0, 0.1, 0.11, 1),
                                                         c("black","grey","white","white")),
                              name = Name_Heatmap,
                              show_row_names = T)

    print(Scal_Heatmap)

    dev.off()

  }


  return_mat <- cbind(t(out_df),
                      P_values = P_values[,1],
                      Padj_values = Padj_values[,1])

  return(return_mat)


}
