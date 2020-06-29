# kinteticMSI
Functions to interpret stable isotope assisted mass spec imaging experiments

## Usage Instructions
kineticMSI has been divided in several steps:

1.  The first step is a calculation of enrichment percentages from the molecular species of interest.

1.  The second step is meant to reconstruct MSI images based on the derived proxies of isotope enrichment.

1.  The third step classifies subclasses inside images based on the enrichment percentages and compares them to anatomical regions of interest obtained thorugh unsupervised statsitcal methods.

1.  The fourth step entails an integrated user-assisted relative quantitation and comparison analyses of the enrichment dynamics of the labelled metabolic targets.

## Step 1 - Enrichment percentages

Enrichment percentages are calculated using the algorithms described in 

  *Heinrich, P., Kohler, C., Ellmann, L., Kuerner, P., Spang, R., Oefner, P. J., and Dettmer, K. (2018). Correcting for natural isotope abundance and tracer impurity in MS-, MS/MS- and high-resolution-multiple-tracer-data from stable isotope labeling experiments with IsoCorrectoR. Sci. Rep. 8.*
  
  *Millard, P., Delépine, B., Guionnet, M., Heuillet, M., Bellvert, F., Létisse, F., and Wren, J. (2019). IsoCor: Isotope correction for high-resolution MS labeling experiments. Bioinformatics 35:4484–4487.*

which correct the endogenous metabolite or peptide pools for the natural isotopic abundance (NIA) according to the chemical formula before calculating via a simple A0 to An division. The IsoCorrectoR is used in R to obtain percentages of enrichment of molecular species with less than 100 hydrogens, and the IsoCor is used in python to obtain percentages of enrichment of molecular species with 100 or more hydrogens (this is due to an intrinsic software limitation within the IsoCorrectoR)

### IsoCorrectoR workflow


### IsoCor workflow


## Step 2 - Spatial dynamics of the tracer

The function depends mainly on Cardinal R package. The function uses percentages of enrichments of specific mass features to reconstruct the tissue slide gaining insigths on which areas of the tissue have incorporated more tracer.

``` 
example_kMSI <- Kinetic_MSI(file = "Imaging_File_Directory", MSI_type = "MSImageSet", Enrichment_File = Enrichment_File_Directory)
``` 

## Step 3 - Clustering regions of differential enrichment 

### drawing the SCC image and selecting an specific cluster as output

There is one mandatory parameter in the function, file_name_WO_extension, which refers to the directory and name of your .ibd and .imzML files. Furthermore there are four optional parameters; type, refers to the kind of Cardinal object that the workflow will use (options are "MSImageSet" and "MSImagingExperiment"). Then you can tune the mathematical parameters of your partitions r, k and s. And the last parameter refers to the cluster_Nr that you want to obtain a data matrix of. Specially usefull to isolate a matrix of your experimental segment

``` 
example_SSC <- segmentation_initial_steps(file_name_WO_extension = "Imaging_File_Directory", type = "MSImageSet", r = 1, k = 5, s = 3, cluster_Nr = 2)
``` 

The output is a segmented image of your file and a data matrix containing the information of your selected cluster.

``` 
example_SSC[1:5,1:3]

                             x = 49, y = 2, z = 1         x = 50, y = 2, z = 1          x = 40, y = 3, z = 1
286.984753095237             7.712556                     8.307419                      7.192391
289.058499875001             5.743113                     5.350798                      5.144504
290.100985837298             0.000000                     0.000000                      3.810359
315.017782775377             0.000000                     0.000000                      0.000000
325.13129680233              0.000000                     0.000000                      0.000000
``` 

## Step 4 - Relative quantitation



