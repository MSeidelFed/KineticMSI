











reconstruct_kMSI <- function(path = "Data/", 
                             as = "MSImageSet" ,
                             position_intensities_legend = "topright",
                             position_cluster_legend = "bottomleft",
                             clust_method = "average",
                             clust_distance = "euclidean",
                             Km_bootstrap = 10,
                             k_means = 5)
  
{
  
  ## functions needed
  
  show_condition <- function(code) {
    tryCatch(code,
             error = function(c) "error",
             warning = function(c) "warning",
             message = function(c) "message"
    )
  }
  
  make_colour_gradient = function(x, brewer_palette = "Spectral") {
    min_x = min(x)
    max_x = max(x)
    range_x = max_x - min_x
    x_scaled = (x - min_x) / range_x
    colours = scales::brewer_pal("seq", brewer_palette)(5)[2:5]
    colour_vals = scales::colour_ramp(colours)(x_scaled)
    colour_vals
  }
  
  #### https://www.rdocumentation.org/packages/qdapTools/versions/1.3.3/topics/list2df
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  } 
  
  
  ## Importing all files
  
  Enrichment_files <- list.files(path = path, pattern = "MeanEnrichment.csv",
                                 full.names = T, recursive = T)
  
  MSI_files <- list.files(path = path, pattern = ".imzML", all.files = FALSE,
                          full.names = T, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  ### testing MSI-enrichment files compatibility
  
  for (l in 1:length(MSI_files)) {
    
    if (grep(pattern = strsplit(strsplit(MSI_files[l],
                                         split = "/")[[1]][2], split = "\\.")[[1]][1],
             x = Enrichment_files) == l) {
      
      runner = T} 
    
    else {
      
      stop("paired files have different names")
      
    }
    
  }
  
  
  ## main function
  
  cluster_lists <- matrix(NA, nrow = 0, ncol = k_means)
  
  for (m in 1:length(Enrichment_files)) {
    
    ### importing data
    
    MSI_file_m <- strsplit(MSI_files[m], split = "\\.")[[1]][1]
    
    enrichment_file <- read.csv(file = Enrichment_files[m],
                                header = T,
                                dec = ".", 
                                row.names = 1)
    
    coords_file <- readImzML(name = MSI_file_m, as = as)
    
    coords_ <- coord(coords_file)
    
    df_coords_ <- as.data.frame(rbind(xAxis = coords_$x, yAxis = coords_$y))
    colnames(df_coords_) <- rownames(coords_)
    
    ### procedure
    
    pdf(file = paste0(path,
                      "sMZ_SubStr",
                      "_" ,
                      strsplit(x = MSI_file_m,
                               split = "/")[[1]][length(strsplit(x = MSI_file_m,
                                                                 split = "/")[[1]])], ".pdf"))
    par(mfrow = c(2,1))
    
    cluster_list <- list()
    
    for (i in 1:dim(enrichment_file)[1]) {
      
      #### setting a boundary for almost empty matrices
      
      if (show_condition(code = column_order(Heatmap(as.numeric(enrichment_file[i,]),
                                                     km = k_means)) == "error") == "error") {
        
        plot(1,1, main = rownames(enrichment_file[i,]))
        plot(1,1, main = rownames(enrichment_file[i,]))}
      
      else {
        
        HM <- Heatmap(as.numeric(enrichment_file[i,]),
                      km = k_means,
                      row_km_repeats = Km_bootstrap,
                      cluster_rows = F,
                      clustering_method_rows = clust_method,
                      clustering_distance_rows = clust_distance,
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
        
        color_vector <- c("tomato", "sandybrown", "lemonchiffon2", "gray70", "gray20")
        
        for (k in 1:length(cluster_number)) {
          
          plot_color_runner <- rep(color_vector[k],
                                   length(cluster_number[names(sort(colMeans(na.omit(cluster_out)),
                                                                    decreasing = T))[k]][[1]]))
          
          plot_colors_ordered <- c(plot_colors_ordered, plot_color_runner)
          
          plot_coords_runner <- coords_[cluster_number[names(sort(colMeans(na.omit(cluster_out)),
                                                                  decreasing = T))[k]][[1]],]
          
          plot_coords_ordered <- rbind(plot_coords_ordered, plot_coords_runner)
          
        }
        
        #### plots
        ##### clustered
        
        plot(x = plot_coords_ordered[,"x"], y = plot_coords_ordered[,"y"],
             col = plot_colors_ordered, pch = 18,
             main = paste0("Cluster (#)", " - ", rownames(enrichment_file)[i]),
             xlim = c(min(as.numeric(df_coords_[1,])) - 10,
                      max(as.numeric(df_coords_[1,])) + 30),
             xlab = "Abscissa coordinates", ylab = "Ordinate coordinates") 
        
        legend(position_cluster_legend, col = unique(plot_colors_ordered), pch = 18,
               legend = names(sort(colMeans(na.omit(cluster_out)), decreasing = T)))
        
        ##### enrichment
        
        plot(x = as.numeric(df_coords_[1,]),
             y = as.numeric(df_coords_[2,]), 
             col = make_colour_gradient(x = as.numeric(enrichment_file[i,])),
             bty = "n",
             type = 'p',
             pch = 15,
             xlim = c(min(as.numeric(df_coords_[1,])) - 10,
                      max(as.numeric(df_coords_[1,])) + 30),
             main = paste0("Enrichment (%)", " - ", rownames(enrichment_file)[i]),
             xlab = "Abscissa coordinates", ylab = "Ordinate coordinates")
        
        legend(position_intensities_legend, pch = 18, col = "white", fill = "white", bty = "o",
               legend = c(paste0(c("Min:",
                                   "1st-Qu:",
                                   "Median:",
                                   "Mean:",
                                   "3rd-Qu:",
                                   "Max:"),
                                 " ",
                                 round(x = as.numeric(summary(as.numeric(enrichment_file[i,]))),
                                       digits = 2))))
        
        #### return object part 1
        
        cluster_list[[i]] <- cluster_out
        
      }
      
    }
    
    dev.off()
    
    
    #### return object part 2
    
    names(cluster_list) <- rownames(enrichment_file)
    
    out_mat <- matrix(NA, nrow = length(cluster_list), ncol = k_means)
    
    names_vector <- c()
    
    for (n in 1:length(cluster_list)) {
      
      names_vector[n] <- paste0(names(cluster_list)[n], "_",
                                strsplit(x = MSI_file_m,
                                         split = "/")[[1]][length(strsplit(x = MSI_file_m,
                                                                           split = "/")[[1]])])
      
      if (is.null(dim(cluster_list[[n]]))) {
        ## filters species for which the k means algorithm fails
        
        out_mat[n,] <- rep(NA, k_means)
        
      } else if (dim(cluster_list[[n]])[2] != k_means){
        ## filters species could only be partitioned into less than the predefined n k means
        
        out_mat[n,] <- rep(NA, k_means)
        
      } else {
        
        out_mat[n,] <- colMeans(cluster_list[[n]], na.rm = T) 
        
      }
      
    }
    
    rownames(out_mat) <- names_vector
    
    cluster_lists <- rbind(cluster_lists, out_mat)
    
  }
  
  return(cluster_lists)
  
}

