# KineticMSI
**Functions to interpret stable isotope assisted mass spec imaging experiments**

This repository only holds the R package itself and the exemplary data to execute it in default mode. For questions, threads and assistance please use our sister repository [KineticMSI_HelpThreads](https://github.com/MSeidelFed/KineticMSI_HelpThreads).

## Introduction
KineticMSI is a collection of scripts to assist in accurate data preparation and analysis of stable isotope assisted (kinetic) Mass Spectrometry Imaging experiments in order to derive functional biological interpretations. Additionally, the functionality available in KineticMSI is compatibile with stable isotope assisted liquid-chromatography mass spectrometry data. See below for a link to the KineticMSI to Kinetic LCMS (KineticMSI_2_kLCMS) repository and installation. The procedures described here are detailed in the [KineticMSI preprint](https://www.biorxiv.org/content/10.1101/2022.08.31.505954v1.full.pdf).

The repo follows this file structure:

1. [Usage Instructions](https://github.com/MSeidelFed/KineticMSI/blob/master/USAGE.md): _detailed and recommended usage of R script code to run the analysis step-by-step._
1. [Data](https://github.com/MSeidelFed/KineticMSI/tree/master/inst/extdata): _sample data used in the original project from which the usage examples are based. Use this to reproduce our results._
1. [R_Functions](https://github.com/MSeidelFed/KineticMSI/tree/master/R): _collection of R scripts to carry out various steps of the analysis._

1. [Images](https://github.com/MSeidelFed/KineticMSI/tree/master/images): _some figures relevant to the repo_

Below is an illustration of the workflow.

**Workflow - Module I and II**

![Workflow - Module I and II](images/Fig2_Modified_workflow_GitHub.png)

## Installation

### First step

It is necessary to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/history.html) before starting. Make sure to install the appropriate version depending on your R version and add the Rtools directory to the system PATH (instructions for this can be found [here](https://datag.org/resources/documents/spring-2018/37-de-barros-installing-r-on-windows/file))

### Installing KineticMSI

```
### Get the latest installation of RandoDiStats (KineticMSI depends on it)

devtools::install_github("MSeidelFed/RandodiStats_package")
library(RandoDiStats)

### Additional packages required for KineticMSI (depending on the R version these may or may not be automatically installed with your KineticMSI installation, better to install them manually beforehand)

install.packages(“BiocManager”)
BiocManager::install("pcaMethods")
BiocManager::install("ComplexHeatmap")
BiocManager::install(“Cardinal”)
BiocManager::install("sva")


### Install KineticMSI

devtools::install_github("MSeidelFed/KineticMSI")
library(KineticMSI)

```

## Getting the exemplary datasets directory after installation

```
system.file("extdata", package = "KineticMSI")
```

## Link to KineticMSI to kLCMS repository and installation

https://github.com/MSeidelFed/KineticMSI_2_kLCMS

