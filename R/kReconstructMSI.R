#' A function to produce graphical reconstructions from MSI images using isotope tracer dynamics
#'
#' This function allows to produce the first type of inference on the replicated datasets that KineticMSI offers. The function produces first boxplots that reflect the dispersion of the chosen isotope tracer proxy across MSI pixels and treatments. The boxplots are produced both for separate replicates and as treatment means. Afterwards the function reconstructs the MSI acquired images by extracting the coordinates from the .ibd .imzML input files and using them to build reconstructed images that reflect the isotope tracer dynamics within the samples specimen. The function has two different run possibilities, it may be run "before" or "after" the KineticMSI class discovery function. In the former case the function will build clusters during the spatial reconstruction of images using a K-mean clustering algorithm and attempting at clustering together pixels based on specific similarity or dissimilarity metrics. In the latter case the clusters are inherited from the kClassDiscoveryMSI function output when the parameter returnObject is set to "minDatasetPlusCoords".
#' @param Reconstruct either "Before" or "After". Allows defining if the reconstruction algorithm runs before or after applying the kClassDiscoveryMSI function.
#' @param kClustersMSI only needs to be defined if the previous parameter is set to "After". the input is a kClassDiscoveryMSI output object produced with the parameter returnObject set to "minDatasetPlusCoords".
#' @param Enrpath parent directory were all enrichment files are contained, digs within recursive folders.
#' @param MSIPath parent directory were all MSI files are contained, digs within recursive folders.
#' @param PatternEnrichment defaults to "MeanEnrichment.csv". character vector used to grab input csv enrichment files
#' @param outpath defaults to getwd(). Defines where the output files are written.
#' @param as inherits from Cardinal, type of MSI experimental file.
#' @param PositionClusterLegend defaults to "bottomleft". Allows to define the position of the legend in the cluster plot.
#' @param clustMethod defaults to "average". Inherits from ComplexHeatmap and supports all clustering methods described there.
#' @param clustDistance defaults to "euclidean". Inherits from ComplexHeatmap and supports all clustering methods described there.
#' @param kmeans defaults to 5. Allows to define the number of clusters looked for in each feature dataset.
#' @param KmBoot defaults to 10. Inherits from ComplexHeatmap and defines the number of bootstrap iterations used to build the K-mean consensus.
#' @param RevAbscissas defaults to FALSE. Allows to reverse the x-axis from its default position.
#' @param RevOrdinates defaults to FALSE. Allows to reverse the y-axis from its default position.
#' @param FactorName Character vector used to define the names of the plots x-axis.
#' @param yLabName Character vector used to define the names of the plots x-axis.
#' @param paletteSpatialPlots defines a palette inherited from viridis to color the intensities plots. Can be any of viridis, magma, plasma, inferno, cividis, mako, rocket, turbo.
#' @param ContrastPercentValue defaults to 0.1. Allows to define the minimum value to build contrast in the intensity plots.
#' @param SubSetRepsIntensities defaults to FALSE. Allows to subset the MSI file list found in path.
#' @param SubSetRepsMSI defaults to FALSE. Allows to subset the csv file list found in path.
#' @param returnObject defaults to TRUE. Allows to return the list of matrices to the R enviroment. When FALSE, only the output files are produced
#' @keywords Image Reconstruction Isotope Tracer Dynamics
#' @export
#' @examples
#'
#' ...


kReconstructMSI <- function(Reconstruct = c("After", "Before"),
                            kClustersMSI,
                            EnrPath,
                            MSIPath = system.file("extdata", package = "KineticMSI"),
                            PatternEnrichment = "MeanEnrichment.csv",
                            outpath = getwd(),
                            as = c("MSImageSet","MSImagingExperiment"),
                            PositionClusterLegend = "bottomleft",
                            clustMethod = "average",
                            clustDistance = "euclidean",
                            kmeans = 5,
                            KmBoot = 10,
                            RevAbscissas = FALSE,
                            RevOrdinates = FALSE,
                            FactorName,
                            yLabName = "Enrichment (%)",
                            paletteSpatialPlots = c(viridis, magma,
                                                    plasma, inferno,
                                                    cividis, mako,
                                                    rocket, turbo),
                            ContrastPercentValue = 0.1,
                            SubSetRepsIntensities = FALSE,
                            SubSetRepsMSI = FALSE,
                            returnObject = TRUE) {

  ## Functions needed
  ### plots function

  plots_recontruct_kMSI <- function(Reconstructs = c("After", "Before"),
                                    df_coords,
                                    RevAbscissa = RevAbscissas,
                                    RevOrdinate = RevOrdinates,
                                    plot_coords,
                                    plot_colors,
                                    positionClusterLegend = PositionClusterLegend,
                                    clusterOut,
                                    enrichmentFile,
                                    paletteChosen,
                                    ContrastPercent = ContrastPercentValue) {
    ##### defining coords
    if (RevAbscissa == F) {

      Abscissas = c(min(as.numeric(df_coords["x",])) - 10,
                    max(as.numeric(df_coords["x",])) + 30)

      AbscissasContour = c(min(as.numeric(df_coords["x",])),
                           max(as.numeric(df_coords["x",])))

    } else {

      Abscissas = c(max(as.numeric(df_coords["x",])) + 30,
                    min(as.numeric(df_coords["x",])) - 10)

      AbscissasContour = c(max(as.numeric(df_coords["x",])),
                           min(as.numeric(df_coords["x",])))
    }

    if (RevOrdinate == F) {

      Ordinates = c(min(as.numeric(df_coords["y",])),
                    max(as.numeric(df_coords["y",])))

    } else {

      Ordinates = c(max(as.numeric(df_coords["y",])),
                    min(as.numeric(df_coords["y",])))

    }

    ##### clustered

    ##### Saving legend

    if(Reconstructs == "After") {

      Unique_colors_minus_gray <-  unique(plot_colors)[unique(plot_colors) != "gray80"]

      plot(x = plot_coords[,"x"], y = plot_coords[,"y"],
           col = plot_colors, pch = 18,
           main = paste0("Cluster (#)", " - ", rownames(enrichmentFile)[i]),
           xlim = Abscissas,
           ylim = Ordinates,
           xlab = "Abscissa coordinates", ylab = "Ordinate coordinates")

      legend(positionClusterLegend, col = Unique_colors_minus_gray, pch = 18,
             legend = names(sort(colMeans(na.omit(clusterOut)), decreasing = T)))

      clustPlot <- recordPlot()

    } else if (Reconstructs == "Before") {

      plot(x = plot_coords[,"x"], y = plot_coords[,"y"],
           col = plot_colors, pch = 18,
           main = paste0("Cluster (#)", " - ", rownames(enrichmentFile)[i]),
           xlim = Abscissas,
           ylim = Ordinates,
           xlab = "Abscissa coordinates", ylab = "Ordinate coordinates")

      legend(positionClusterLegend, col = unique(plot_colors), pch = 18,
             legend = names(sort(colMeans(na.omit(clusterOut)), decreasing = T)))

      clustPlot <- recordPlot()

    } else {

      stop("Reconstruct should be either after or before kClassDiscoveryMSI")

    }

    ##### enrichment
    min_contrast_value = (min(as.numeric(enrichmentFile[i,]))*100)-ContrastPercent

    mat_contour <- matrix(min_contrast_value,
                          nrow = max(as.numeric(df_coords["y",])),
                          ncol = max(as.numeric(df_coords["x",])))

    for (ii in 1:dim(df_coords)[2]) {

      x_runner <- df_coords["x",ii]
      y_runner <- df_coords["y",ii]

      if(is.na(round(as.numeric(enrichmentFile[i,])[ii]*100, 2))) {

        stop("Your Enrichment Files Have Different Coord No. compared to your .imzML .ibd files")

      } else if (show_condition(code = mat_contour[y_runner, x_runner] <- round(as.numeric(enrichmentFile[i,])[ii]*100, 2)) == "error") {

        stop("Your Enrichment Files Have Different Coord No. compared to your .imzML .ibd files")

      }

      mat_contour[y_runner, x_runner] <- round(as.numeric(enrichmentFile[i,])[ii]*100, 2)

    }

    in_contour <- t(mat_contour)

    filled.contour(x = 1:nrow(in_contour),
                   y = 1:ncol(in_contour),
                   in_contour,
                   color.palette=paletteChosen,
                   ylim = Ordinates,
                   xlim = AbscissasContour,
                   plot.title = title(main = rownames(enrichmentFile)[i]),
                   key.title = title(main = yLabName, cex.main = 0.65),
                   axes=T,
                   plot.axes = T)

    contourPlot <- recordPlot()

    ## plotting
    clustPlot
    contourPlot

  }

  ### catching errors function
  show_condition <- function(code) {
    tryCatch(code,
             error = function(c) "error",
             warning = function(c) "warning",
             message = function(c) "message"
    )
  }

  #### https://www.rdocumentation.org/packages/qdapTools/versions/1.3.3/topics/list2df
  list2df <- function(x)
  {
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE)
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x)))))
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")
    DF
  }




  ## Main
  ### boxplots generation bit of the function
  #### Importing all files, pairing and subseting them
  ##### Enrichment files

  Enrichment_files <- list.files(path = EnrPath, pattern = PatternEnrichment,
                                 all.files = FALSE, full.names = TRUE, recursive = TRUE,
                                 ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)

  if (SubSetRepsIntensities == T) {

    Enrichment_files = SubSetInputFiles(reps = Enrichment_files)

  }

  ##### MSI files

  MSI_files <- list.files(path = MSIPath, pattern = ".imzML",
                          all.files = FALSE, full.names = TRUE, recursive = TRUE,
                          ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)

  if (SubSetRepsMSI == T) {

    MSI_files = SubSetInputFiles(reps = MSI_files)

  }

  #### testing MSI-enrichment files compatibility

  MSI_file_names <- c()

  out_vector_names <- c()

  for (l in 1:length(MSI_files)) {

    dir_elements <- strsplit(MSI_files[l], split = "/")[[1]]

    dir_file <- strsplit(dir_elements[length(dir_elements)], split = "\\.")[[1]][1]

    MSI_file_names[l] <- dir_file

    if (length(grep(pattern = dir_file, x = Enrichment_files) != l) == 0){

      stop("paired files have different names")

    } else if (grep(pattern = dir_file, x = Enrichment_files) != l) {

      stop("paired files have different names")

    }

    #### getting the min dataset in terms of minimum pixels

    file_runner <- read.csv(file = Enrichment_files[l], header = T, row.names = 1)

    out_vector_names[l] <- dim(file_runner)[2]

  }

  #### Using the min dataset and the entities in common to build boxplots per individual entity
  ##### common entities

  entity_names <- list()

  for (i in 1:length(Enrichment_files)) {

    entity_names[[i]] <- rownames(read.csv(Enrichment_files[i], header = T, row.names = 1))

  }

  entity_common_names <- Reduce(intersect, entity_names)

  ##### min dataset in terms of pixels and entities

  min_dataset <- min(out_vector_names)

  out_mat_minSet <- matrix(NA, nrow = 0, ncol = min_dataset)

  for (i in 1:length(Enrichment_files)) {

    runner_file <- read.csv(Enrichment_files[i], header = T, row.names = 1)

    runner_sampled <- as.matrix(sample(runner_file, min_dataset, replace = F))

    runner_sampled_common <- runner_sampled[entity_common_names,]

    out_mat_minSet <- rbind(out_mat_minSet, runner_sampled_common)

  }

  lipid_Nr <- unique(rownames(out_mat_minSet))

  #### Producing the boxplots

  pdf(file = paste0(outpath, "/", "Boxplots", ".pdf"))
  par(mfrow = c(2,1))

  cat("...\n")
  cat("Building Boxplots\n")
  cat("...\n")

  for (i in 1:length(lipid_Nr)) {

    ##### grabbing all the reps for a single entity

    Single_lipid <- out_mat_minSet[grep(pattern = lipid_Nr[i], x = rownames(out_mat_minSet)),]

    rownames(Single_lipid) <- paste0(rep(lipid_Nr[i],
                                         length(MSI_file_names)), "_", MSI_file_names)

    test_to_melt <- reshape2::melt(Single_lipid)

    ##### plotting

    Factor <- rep(apply(as.data.frame(strsplit(MSI_file_names, "_"))[1,], 2, as.factor),
                  min_dataset)

    names(Factor) <- NULL

    a_reps <- ggplot2::ggplot(test_to_melt, ggplot2::aes(x=Var1, y=value, color= Factor)) +
              ggplot2::geom_boxplot() +
              ggplot2::coord_flip() +
              ggplot2::geom_jitter(shape = 1, position = ggplot2::position_jitter(0.2))+ ggplot2::ylab(yLabName) +
              ggplot2::xlab(FactorName)+
              ggplot2::theme_classic()+
              ggplot2::ggtitle(lipid_Nr[i])

    test_to_melt[,"Var1"] <- Factor

    a_mean <- ggplot2::ggplot(test_to_melt, ggplot2::aes(x=Var1, y=value, color= Factor)) +
              ggplot2::geom_boxplot() +
              ggplot2::coord_flip() +
              ggplot2::geom_jitter(shape = 1, position = ggplot2::position_jitter(0.2))+ ggplot2::ylab(yLabName) +
              ggplot2::xlab(FactorName)+
              ggplot2::theme_classic()+
              ggplot2::ggtitle(lipid_Nr[i])

    gridExtra::grid.arrange(a_reps, a_mean, nrow = 2)

  }

  dev.off()



  ### reconstruct bit of the function

  cluster_list_of_lists <- list()

  for (m in 1:length(Enrichment_files)) {

    #### importing enrichment data

    cat("...\n")
    cat(paste0("Reconstructing kMSI: ", MSI_file_names[m]))
    cat("...\n")

    enrichment_file <- read.csv(file = Enrichment_files[m],
                                header = T,
                                dec = ".",
                                row.names = 1)

    ##### building a dataset with the entities in common

    enrichment_file <- enrichment_file[entity_common_names,]

    #### importing MSI files to get coords

    if (as == "MSImageSet") {

      path_m <- strsplit(strsplit(MSI_files[1], "/")[[1]][length(strsplit(MSI_files[1], "/")[[1]])], "\\.")[[1]][1]

      coords_file <- Cardinal::readImzML(name = path_m,
                                         folder = MSIPath,
                                         as = "MSImageSet")

      coords_ <- Cardinal::coord(coords_file)

      df_coords_ <- as.data.frame(rbind(x = coords_$x, y = coords_$y))
      colnames(df_coords_) <- rownames(coords_)

    } else if (as == "MSImagingExperiment") {

      path_m <- strsplit(strsplit(MSI_files[1], "/")[[1]][length(strsplit(MSI_files[1], "/")[[1]])], "\\.")[[1]][1]

      coords_file <- Cardinal::readImzML(name = path_m,
                                         folder = MSIPath,
                                         as = "MSImagingExperiment")

      File_coords1 <- Cardinal::coord(coords_file)

      coords_ <- as.data.frame(cbind(x = File_coords1$x, y = File_coords1$y),
                               row.names = paste0("x =" , ", ",
                                                  File_coords1$x,
                                                  "y =" , ", ",
                                                  File_coords1$y,
                                                  "z =" , ", ",
                                                  File_coords1$z))

      df_coords_ <- as.data.frame(rbind(x = File_coords1$x, y = File_coords1$y))
      colnames(df_coords_) <- rownames(coords_)


    } else {

      stop("only MSImagingExperiment or MSImageSet classes are supported")

    }

    #### Reconstructing MSI images with enrichment percentage

    pdf(file = paste0(outpath, "/", "sMZ_SubStr", "_" , MSI_file_names[m], ".pdf"), height = 3)

    cluster_list <- vector(mode = "list", length = dim(enrichment_file)[1])

    ##### building single reconstructed plots within each file

    for (i in 1:dim(enrichment_file)[1]) {

      ##### After or Before class_discovery_kMSI ? (After inherits clusters from our algorithms)

      if (Reconstruct == "After") {

        ##### Filtering out null objects

        if (is.null(kClustersMSI[[i]][["Coordinates"]])){

          #### NULL objects by default

        } else if (rownames(enrichment_file)[i] != names(kClustersMSI)[i]) {

          ##### stopping if the order of entities in enrichment and kClustersMSI is different

          stop("ERROR: different entities in kClassDiscoveryMSI object and enrichment tables...")

        } else {

          ##### getting cluster coordinates

          clusters_after <- names(kClustersMSI[[i]][["Coordinates"]])

          ###### counting clusters

          ###### producing a cluster_out object

          max_clust_length <- c()

          for (j in 1:length(clusters_after)) {

            max_clust_length <- c(max_clust_length, length(kClustersMSI[[i]][["Values"]][[clusters_after[j]]][m,]))

          }

          ##### building empty objects to fill in for loop

          coords_clusters_featurei_filem <- c()

          cluster_out_A <- matrix(NA, nrow = max(max_clust_length), ncol = length(clusters_after))

          values_cluster_colors_featurei_filem <- c()

          names_ids <- c()

          for (j in 1:length(clusters_after)) {

            coords_after_filei <- kClustersMSI[[i]][["Coordinates"]][[clusters_after[j]]][m,]

            coords_clusters_featurei_filem <- c(coords_clusters_featurei_filem, coords_after_filei)

            #### exchanging individual values by the mean value and color that represents the cluster..

            cluster_out_A[,j] <- c(kClustersMSI[[i]][["Values"]][[clusters_after[j]]][m,],
                                   rep(NA, max(max_clust_length) - length(kClustersMSI[[i]][["Values"]][[clusters_after[j]]][m,])))

            color_vector <- viridis_pal(option = "D")(length(clusters_after))

            values_cluster_colors_featurei_filem <- c(values_cluster_colors_featurei_filem,
                                                      rep(color_vector[j], length(coords_after_filei)))

            #### grabbing the names for later cluster identification

            names_ids <- c(names_ids, paste0(names(coords_after_filei), "_", clusters_after[j]))

          }

          #### naming the output objects

          colnames(cluster_out_A) <- clusters_after

          names(coords_clusters_featurei_filem) <- names_ids

          names(values_cluster_colors_featurei_filem) <- names_ids

          #### producing the input files for the plotting function

          #### colors & coords: Clust Nr + 1 (pixel background)

          color_background <- "gray80"

          empty_color_vector <- rep(color_background, length(colnames(enrichment_file)))

          for (j in 1:length(coords_clusters_featurei_filem)) {

            indexes_runner <- which(colnames(enrichment_file) == coords_clusters_featurei_filem[j])

            empty_color_vector[indexes_runner] <- values_cluster_colors_featurei_filem[j]

          }

          #### Final objects named for plot

          plot_colors_ordered <- empty_color_vector ## plot_colors_ordered ## 1 color per coordinate

          plot_coords_ordered <- t(df_coords_) ##  plot_coords_ordered ## coords in X Y Z format

          cluster_out <- cluster_out_A

          cluster_list[[i]] <- cluster_out

          #### plots

          plots_recontruct_kMSI(Reconstructs = Reconstruct,
                                df_coords = df_coords_,
                                RevAbscissa = RevAbscissas,
                                RevOrdinate = RevOrdinates,
                                plot_coords = plot_coords_ordered,
                                plot_colors = plot_colors_ordered,
                                positionClusterLegend = PositionClusterLegend,
                                clusterOut = cluster_out,
                                enrichmentFile = enrichment_file,
                                paletteChosen = paletteSpatialPlots,
                                ContrastPercent = ContrastPercentValue)

        }

      } else if (Reconstruct == "Before") {

        #### Before class_discovery_kMSI, i.e., conventional K-means clustering

        ##### setting a boundary for almost empty matrices

        if (show_condition(code = column_order(Heatmap(as.numeric(enrichment_file[i,]),
                                                       km = kmeans)) == "error") == "error" |
            mean(apply(enrichment_file[i,], 2, as.numeric)) == 0) {

          ###### null plots

          plot(1,1, main = rownames(enrichment_file[i,]))
          plot(1,1, main = rownames(enrichment_file[i,]))

          ###### out object

          cluster_out <- as.data.frame(matrix(NA, nrow = 1, ncol = kmeans))

        } else {

          HM <- Heatmap(as.numeric(enrichment_file[i,]),
                        km = kmeans,
                        row_km_repeats = KmBoot,
                        cluster_rows = F,
                        clustering_method_rows = clustMethod,
                        clustering_distance_rows = clustDistance,
                        column_title = rownames(enrichment_file[i,]))

          cluster_number <- row_order(HM)

          #### loop to take out a matrix with cluster membership

          out <- list()

          for (j in 1:length(cluster_number)) {

            runner <- as.matrix(enrichment_file)[i, cluster_number[[j]]]

            out[[j]] <- runner

            names(out[[j]]) <- NULL

          }

          cluster_out <- list2df(out[1:length(out)])
          colnames(cluster_out) <- names(cluster_number)

          #### loop to assemble the plot elements

          plot_colors_ordered <- c()

          plot_coords_ordered <- matrix(NA, nrow = 0, ncol = dim(coords_)[2])

          color_vector <- viridis_pal(option = "D")(length(cluster_number))

          for (k in 1:length(cluster_number)) {

            ## sorting clusters according to their mean magnitudes

            plot_color_runner <- rep(color_vector[k],
                                     length(cluster_number[names(sort(colMeans(na.omit(cluster_out)),
                                                                      decreasing = T))[k]][[1]]))

            plot_colors_ordered <- c(plot_colors_ordered, plot_color_runner)

            plot_coords_runner <- coords_[cluster_number[names(sort(colMeans(na.omit(cluster_out)),
                                                                    decreasing = T))[k]][[1]],]

            plot_coords_ordered <- rbind(plot_coords_ordered, plot_coords_runner)

          }

          #### plots

          plots_recontruct_kMSI(Reconstructs = Reconstruct,
                                df_coords = df_coords_,
                                RevAbscissa = RevAbscissas,
                                RevOrdinate = RevOrdinates,
                                plot_coords = plot_coords_ordered,
                                plot_colors = plot_colors_ordered,
                                positionClusterLegend = PositionClusterLegend,
                                clusterOut = cluster_out,
                                enrichmentFile = enrichment_file,
                                paletteChosen = paletteSpatialPlots,
                                ContrastPercent = ContrastPercentValue)

        }

        #### return object unfiltered list

        cluster_list[[i]] <- cluster_out


      } else {

        stop("User needs to define if the reconstruction is Before or After ClassDiscovery kMSI...")

      }

    }

    dev.off()

    #### return object filtered list featuring only clustered entities

    names(cluster_list) <- rownames(enrichment_file)

    out_mat <- matrix(NA, nrow = length(cluster_list), ncol = kmeans)

    names_vector <- c()

    for (n in 1:length(cluster_list)) {

      names_vector[n] <- paste0(names(cluster_list)[n], "_", MSI_file_names[m])

      if (is.null(dim(cluster_list[[n]]))) {
        #### filters species for which the k means algorithm fails

        out_mat[n,] <- rep(NA, kmeans)

      } else if (dim(cluster_list[[n]])[2] != kmeans){
        #### filters species could only be partitioned into less than the predefined n k means

        out_mat[n,] <- rep(NA, kmeans)

      } else {

        out_mat[n,] <- colMeans(cluster_list[[n]], na.rm = T)

      }

    }

    rownames(out_mat) <- names_vector

    ### this object is a matrix that can be used as a merged format of the filtered list if necessary

    cluster_list_of_lists[[m]] <- cluster_list

  }

  if(returnObject == T) {

    return(cluster_list_of_lists)

  } else {

    cat("Finished building the output plots")

  }

}
