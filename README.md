# kineticMSI
Functions to interpret stable isotope assisted mass spec imaging experiments

## Usage Instructions
kineticMSI has been divided in several steps:

1.   Step 1 – data preprocessing

 * Filtering step 1 – All pixels with M0 = 0 are deleted (replaced with NA).
 
 * Filtering step 2 – All pixels with all isotopologues = 0 are deleted (replaced with NA).


# unoptimized
1. The first step is a calculation of enrichment percentages from the molecular species of interest.

1.  The second step is meant to reconstruct MSI images based on the derived proxies of isotope enrichment.

1.  The third step classifies the coordinates from individual molecular species in subclasses based on the enrichment percentages and compares them to anatomical regions of interest obtained thorugh unsupervised statistical methods.

1.  The fourth step entails an integrated user-assisted relative quantitation and comparison analyses of the enrichment dynamics of the labelled metabolic targets.

## Step 1 - Enrichment calculations

Enrichment percentages are calculated using the algorithms described in 

  *Heinrich, P., Kohler, C., Ellmann, L., Kuerner, P., Spang, R., Oefner, P. J., and Dettmer, K. (2018). Correcting for natural isotope abundance and tracer impurity in MS-, MS/MS- and high-resolution-multiple-tracer-data from stable isotope labeling experiments with IsoCorrectoR. Sci. Rep. 8.*
  
  *Millard, P., Delépine, B., Guionnet, M., Heuillet, M., Bellvert, F., Létisse, F., and Wren, J. (2019). IsoCor: Isotope correction for high-resolution MS labeling experiments. Bioinformatics 35:4484–4487.*

which correct the endogenous metabolite or peptide pools for the natural isotopic abundance (NIA) according to the chemical formula before calculating via a simple A0 to An division. The IsoCorrectoR is used in R to obtain percentages of enrichment of molecular species with less than 100 hydrogens, and the IsoCor is used in python to obtain percentages of enrichment of molecular species with 100 or more hydrogens (this is due to an intrinsic software limitation within the IsoCorrectoR)

### IsoCorrectoR workflow
IsoCorrectoR has been installed and used according to the instructions provided upon releasing of the package in BioConductor (https://www.bioconductor.org/packages/release/bioc/html/IsoCorrectoR.html) as follows:

```
Enrichment <- IsoCorrectoR::IsoCorrection(MeasurementFile = "Data/MeasurementFile_HD_D8.csv",
                                     ElementFile = "Data/ElementFile.csv",
                                     MoleculeFile = "Data/MoleculeFile_WT_D8.csv",
                                     CorrectTracerImpurity = T,
                                     CorrectTracerElementCore = T,
                                     CalculateMeanEnrichment = T,
                                     UltraHighRes = F,
                                     FileOut = "trial.csv",
                                     FileOutFormat = "csv",
                                     ReturnResultsObject = T,
                                     CorrectAlsoMonoisotopic = T)
```

For MSI, each column in the input table "MeasurementFile.csv" belongs to a single coordinate on the original image where the isotopologues could be measured and mined out. The files were built manually from our own universal input csv format file (Data/Universal_Isotopologue_File.csv).

### IsoCor workflow

Python version: 
```
module add devel/Python-3.8.0
```

IsoCor input tables have the following format:

 sample              | metabolite | derivative | isotopologue  | area           | resolution  |
 ------------------- | -----------|------------| --------------| ---------------| ------------|
  [1,] "Sample_X"    | "PIP2492"  | ""         | "0"           |  "52335.21982" | "70000"     |
  [2,] "Sample_X"    | "PIP2492"  | ""         | "1"           |  "75684.458"   | "70000"     |
  [3,] "Sample_X"    | "PIP2492"  | ""         | "2"           |  "0"           | "70000"     |
  [4,] "Sample_X"    | "PIP2492"  | ""         | "3"           |  "0"           | "70000"     |
  [5,] "Sample_X"    | "PIP2492"  | ""         | "4"           |  "0"           | "70000"     |
  [6,] "Sample_X"    | "PIP2492"  | ""         | "5"           |  "0"           | "70000"     |

and can be obtained from our universal input csv format file (Data/Universal_Isotopologue_File.csv) using the function ProduceIsoCorTables as follows:

```
WT1 <- ProduceIsoCorTables(PathToCSV_file = "Replicate 1/Pyr layer replicate1  29WT.csv")
```
Subsequently, the output table can be directly used by IsoCor following the published instructions (https://github.com/MetaSys-LISBP/IsoCor).

## Step 2 - Spatial dynamics of the tracer

Following the enrichment calculation procedure, and aiming at taking advantage of the gained spatial dimensions provided by MSI, we explored the kinetics of the tracer in the tissue. To do that we built a fuction that uses percentages of enrichment from specific mass features to reconstruct the tissue slide gaining insigths on which tissue areas, if at all, have incorporated more tracer. The function has a dependency to the Cardinal R package (https://www.bioconductor.org/packages/release/bioc/html/Cardinal.html).

``` 
if(!require(Cardinal)) {BiocManager::install("Cardinal"); require(Cardinal)}
library(Cardinal)

example_kMSI <- Kinetic_MSI(file = "Imaging_File_Directory", MSI_type = "MSImageSet", Enrichment_File = Enrichment_File_Directory)
``` 

The outcome from this visual evaluation is a prioritization of peptides or metabolites for the subsequent steps.

## Step 3 - Clustering active regions with differential tracer incorporation 

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



