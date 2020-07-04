# GRAPPLE


This is the R package for GRAPPLE (Genome-wide  mR  Analysis  under  Pervasive  PLEiotropy), a comprehensive framework for Mendelian Randomization. 

![](http://jingshuw.org/uploads/1/2/2/1/122138403/grapple.png)
For details of the GRAPPLE framework, please refer to our [manuscript](https://www.biorxiv.org/content/10.1101/2020.05.06.077982v1).


## Installation
The GRAPPLE package can be installed in R from this Github repository:

```
library(devtools)
install_github("jingshuw/grapple")
```

One may also need to install [PLINK](https://www.cog-genomics.org/plink/) for LD clumping, and prepare a suitable LD clumping reference dataset for PLINK, if he/she requires GRAPPLE to perform data preprocessing. One resource of the reference dataset of the European population for clumping, which is also used in the GRAPPLE paper, is the 1000 genome European reference panel downloaded here http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz, from the [MRCIEU cite](https://github.com/MRCIEU/gwas2vcf).

## Input Datasets


The analysis with GRAPPLE starts with preprocessing the raw GWAS summary statistics, which includes selection of genetic instruments, harmonizing the datasets to share the same effect and reference allele, calculating a correlation matrix of the noise in estimated marginal associations which is shared among SNPs, and selecting a list of candidate SNPs for mode marker detection.

In order to do this, GRAPPLE needs a list of raw GWAS summary statistics files. Each GWAS file is a ".csv" or ".txt" file containing a data frame of at least 6 columns with these column names: 'SNP', 'effect allele', 'other allele', 'beta', 'se', 'pval'. The SNP column contains rsID for each SNP. Both the effect allele and other allele column need to have capital letters. The beta column contains the estimated effect size for continuous trait and log odds ratio for binary trait, and the se column is the standard deviation of the corresponding 'beta'. 

All GWAS summary statistics files in the above format that are used in the GRAPPLE paper can be download from [datasets](https://www.dropbox.com/sh/vv6pz09cknyz9ca/AAAV_WWLsJmI2LZwL1da45q0a?dl=0). For the original publicly available source of each file, please refer to Table S1 and S2 of our [manuscript](https://www.biorxiv.org/content/10.1101/2020.05.06.077982v1).                             


## Data Preprocessing

  Deta preprocessing of GRAPPLE needs one file for the disease (outcome file), and two files for each risk factor, where one is for SNP selection (selection file) and the other is for estimating the marginal effects of the SNPs on this risk factor (exposure file). We allow overlapping individuals for the GWAS data in the outcome and exposure files, which there should be as few overlapping individuals as possible between the selection files an other files to avoid SNP selection bias. 



 As an example, one can download PLINK and save the executable file as "./plink", and also download the reference panel and unzip the ".tgz" file in the current folder.


```
library(grapple)
sel.file <- "bmi_ukbb_full.csv"
exp.file <- "bmi_giant_full.csv"
out.file <- "t2d_morris2012_impute.csv"
plink.exe <- "./plink"
plink_refdat <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"
data <- getInput(sel.file, exp.file, out.file, plink.exe, plink_refdat)
```

## Basic Usage

If the number of risk factors is 1, GRAPPLE can detect the number of pleiotropic pathways by finding the number of modes in the robustified profile likelihood, given a p-value selection threshold.
```
diagnosis <- findModes(data, p.thres = 1e-4)
diagnosis$modes
diagnosis$p
```

If there are multiple modes, then we can check which SNP markers contribute to the mode. To correct for multiple pleiotropic pathways, one can either remove outliers or add more risk factors. GRAPPLE can also estimate the causal effect of the risk factors assuming a random effect model on SNPs' pleiotropic effects.
```
result <- grappleRobustEst(data, p.thres = 1e-4)
```
