# <ins>D</ins>enoising <ins>u</ins>sing <ins>Re</ins>plicate <ins>S</ins>pectra (DuReS)

## Introduction

This package provides easy-to-implement functions to denoise tandem mass spectrometry data. It requires a set of mzML files and a txt file containing feature information (from standard untargeted metabolomics software), such as the precursor mz and RT as input. It outputs a set of mzML files with the same number of samples but containing denoised MS/MS spectra. 

## Installation

First, you need to install two dependencies from Bioconductor: S4Vectors and Spectra
```r
install.packages("BiocManager")
BiocManager::install(c("Spectra", "S4Vectors", "mzR"))
```
After this you can proceed with the installation of the development version of the package DuReS as follows:

```r
install.packages("devtools")
devtools::install_github("BiosystemEngineeringLab-IITB/dures")
```

## Documentation
- A detailed description of all the functions included within the package is available [here](https://biosystemengineeringlab-iitb.github.io/dures/reference/index.html). 
- A vignette outlining how to run DuReS using several test datasets is available [here](https://biosystemengineeringlab-iitb.github.io/dures/articles/dures-vignette.html).
- Another vignette describing how to derive the optimal recurrence frequency cutoff is available [here](https://biosystemengineeringlab-iitb.github.io/dures/articles/dures-vignette-tuning.html).  

## Test Datasets
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13778168.svg)](https://doi.org/10.5281/zenodo.13778168)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15132924.svg)](https://doi.org/10.5281/zenodo.15132924).



## Quick Start
```r
#This step reads in the mzML files, prepares the stats.txt file in a format that extracts MS2 spectra and returns a list
folder_path = "~/metabolomics/test_1/" #folder path containing mzml/ and Stats.txt in required format
l1 = preprocess(folder_path = folder_path, tol_mz = 5, tol_rt = 0.1) #reads mzml files, prepares Stats file, extracts spectra and concatenates spectra

#This step extracts the top 80% TIC spectra and groups fragments within a given mass tolerance
l2 = extract_raw_spectra(folder_path = folder_path, l1, 0.05, 0.8) #extract top x% (where x = 0.8) TIC spectra, groups fragments within a given tolerance (0.05 Da)

#This step aggregates the top 80% TIC spectra from step2 and calculates the fragment frequencies
l3 = call_aggregate(l2$sps_top_tic_2, 0.05, folder_path) 

#This step labels the individual spectrum with frequencies learnt from Step 3
l4 = label_individual_spectrum(l3, folder_path, 0.05)

#This step removed fragments with frequencies below the given threshold (denoising step)
l5 = generate_denoised_spectra(l4, folder_path, ion_mode = "pos") 
```

## Tuning the optimal recurrence frequency parameter

DuReS provides users with the ability to determine the optimal frequency cutoff for denoising MS/MS spectra. A detailed walkthrough is available through an example analysis using an open-source experimental tandem mass spectrometry dataset [here](https://biosystemengineeringlab-iitb.github.io/dures/articles/dures-vignette-tuning.html).

## Metabolite annotation

**DuReS** includes a curated reference library built from *11 publicly available MS/MS datasets*, comprising:

- **1,259,372 unique spectra**  
- **356,330 unique compounds**  
- Coverage across both **positive** and **negative** ionization modes  
- Redundant entries removed using **SPLASH keys** (see Supporting Information, Section 3 of the main manuscript)

The processed positive and the negative library spectra is available here. Users may choose to install GIT LFS for using the files in an automated manner (see vigentte titled [Parameter Tuning](https://biosystemengineeringlab-iitb.github.io/dures/articles/dures-vignette-tuning.html) or may directly upload the files from this [link](https://drive.google.com/drive/folders/1CIlmggvodtsPSUyxY1mL_M1UazSOHvGN?usp=sharing) to the package directory /inst/extdata/.

Experimental and denoised spectra are matched against this reference library using the following scoring function:

```
Matching Score = 5 × [Modified Dot Product × 100 + 20 × log₂(max(NMF, 1))]
```

Where:  
- **Modified Dot Product** measures spectral similarity  
- **NMF** = Number of Matching Fragments

In addition to the Matching Score, DuReS reports:
- **Forward and Reverse Dot Products**
- **Fragment Matching Ratio (FMR)**
- Number of matched fragments
- Reference compound identifiers

To ensure high-confidence annotations:
- Matches with fewer than **2 fragments** or **dot products < 0.25** are filtered out  
- Annotations must occur in **≥30%** of replicate spectra

All parameters are fully **tunable**, and the spectral matching algorithm is **integrated into DuReS**.

## Package Workflow

![Workflow Diagram](https://github.com/BiosystemEngineeringLab-IITB/dures/blob/main/images/Denoising_workflow_tuning_testing_combined.drawio-1.png)
=======

## Citation
Banerjee, Shayantan, Prajval Nakrani, Aviral Singh, and Pramod Wangikar. "DuReS: An R package for denoising experimental tandem mass spectrometry-based metabolomics data." [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.16.613198v1) (2024): 2024-09.
