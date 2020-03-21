# kinteticMSI
functions to interpret stable isotope assisted mass spec imaging experiments

## Distribution of the tracer across the tissue

The function depends mainly on Cardinal R package. The function uses percentages of enrichments of specific mass features to reconstruct the tissue slide gaining insigths on which areas of the tissue have incorporated more tracer.

``` 
example_kMSI <- Kinetic_MSI(file = "Imaging_File_Directory", MSI_type = "MSImageSet", Enrichment_File = Enrichment_File_Directory)
``` 


## drawing the SCC image and selecting an specific cluster as output

``` 
example_SSC <- segmentation_initial_steps(file_name_WO_extension = "Imaging_File_Directory", type = "MSImageSet", r = 1, k = 5, s = 3, cluster_Nr = 2)
``` 

