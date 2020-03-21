segmentation_initial_steps <- function(file_name_WO_extension = "85wt_hippocampus_female_rep2", type = "MSImageSet", r=1, k=5, s=3, cluster_Nr = 0) {
  
if(!require(Cardinal)) {BiocManager::install("Cardinal"); require(Cardinal)}
library(Cardinal)
  
msset <- readImzML(name = file_name_WO_extension, attach.only = FALSE, as = type)
msset2 <- normalize(msset, method="tic")
msset3 <- smoothSignal(msset2, method="gaussian", window=9)
msset4 <- reduceBaseline(msset3, method="median", blocks=50)
msset5 <- msset4[,!duplicated(coord(msset4))]
msset6 <- peakPick(msset5, method="adaptive")
msset7 <- peakAlign(msset6, method="diff")
msset8 <- peakFilter(msset7, method="freq")

ssc <- spatialShrunkenCentroids(msset8, r=r, k=k, s=s, method="adaptive")

## segmented image output
image(ssc, col=c("black", "red", "orange", "grey", "blue"), key=TRUE, superpose = T, strip = F)


## matrix of one cluster output
if (cluster_Nr > 0) {

  File <-  ifelse(as.data.frame(ssc$classes) == cluster_Nr, yes = TRUE, no = FALSE)

File_spectra <- spectra(msset8)
File_colour_Spectra <- File_spectra[,File[,1]]

File_coords <- coord(msset8)
File_colour_coords <- File_coords[File[,1],]

File_mzs <- mz(msset8)

File_Mat <- as.matrix(File_colour_Spectra)
colnames(File_Mat) <- rownames(File_colour_coords)
rownames(File_Mat) <- File_mzs

return(File_Mat)
}

}
