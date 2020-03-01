tracer_dynamics <- function(file,
                        MSI_type = MSImageSet,
                        Enrichment_File)
  
  {
    
if(!require(Cardinal)) {BiocManager::install("Cardinal"); require(Cardinal)}
library(Cardinal)

getting_mat <- function(msset) {
msset2 <- normalize(msset, method="tic")
msset3 <- smoothSignal(msset2, method="gaussian", window=9)
msset4 <- reduceBaseline(msset3, method="median", blocks=50)
msset5 <- msset4[,!duplicated(coord(msset4))]
msset6 <- peakPick(msset5, method="adaptive")
msset7 <- peakAlign(msset6, method="diff")
msset8 <- peakFilter(msset7, method="freq")
}

make_colour_gradient = function(x, brewer_palette = "Dark2") {
    min_x = min(x)
    max_x = max(x)
    range_x = max_x - min_x
    x_scaled = (x - min_x) / range_x
    colours = scales::brewer_pal("seq", brewer_palette)(5)[2:5]
    colour_vals = scales::colour_ramp(colours)(x_scaled)
    colour_vals
}

msset <- readImzML(name = file , as = MSI_type)
mat <- getting_mat(msset = msset)
File_spectra <- spectra(mat)
File_coords <- coord(mat)
File_mzs <- mz(mat)

File_Mat <- as.matrix(File_spectra)
colnames(File_Mat) <- rownames(File_coords)
rownames(File_Mat) <- File_mzs

common_mzs <- as.numeric(rownames(File_Mat))
mycol <- gradient.colors(100, start="white", end="black") 


pdf(file = "Best_mzs.pdf")

for (i in 2:length(common_mzs)) {
  image(msset, mz= common_mzs[i], col.regions=mycol, smooth.image="gaussian")
}

dev.off()

coords <- File_coords
Enrichment <- read.csv(file = Enrichment_File, header = T, dec = ".", row.names = 1)[,c(-dim(read.csv(file = Enrichment_File, header = T, dec = ".", row.names = 1))[2])]
df <- as.data.frame(rbind(xAxis = File_coords$x, yAxis = File_coords$y))
colnames(df) <- rownames(File_coords)

pdf(file = "Enrichment_Reconstructions.pdf")

for (i in 1:dim(Enrichment)[1]) {
  
  plot(x = as.numeric(df[1,]),
     y = as.numeric(df[2,]), 
     col = make_colour_gradient(x = as.numeric(Enrichment[i,])),
     bty = "n",
     type = 'p',
     pch = 15,
     xlim = c(-10,250),
     main = rownames(Enrichment)[i])

    }

dev.off()


return(File_Mat)

}
