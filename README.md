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

In order to do this, GRAPPLE needs a list of raw GWAS summary statistics files. Each GWAS file is a ".csv" or ".txt" file containing a data frame of at least 6 columns with these column names: `SNP`, `effect allele`, `other allele`, `beta`, `se`, `pval`. The `SNP` column contains rsID for each SNP. Both the `effect allele` and `other allele` column need to have capital letters. The `beta` column contains the estimated effect size for continuous trait and log odds ratio for binary trait, and the `se` column is the standard deviation of the corresponding `beta`. 

All GWAS summary statistics files in the above format that are used in the GRAPPLE paper can be download from [datasets](https://www.dropbox.com/sh/vv6pz09cknyz9ca/AAAV_WWLsJmI2LZwL1da45q0a?dl=0). For the original publicly available source of each file, please refer to Section 3 of the manuscript's [Supplementary Information](https://www.biorxiv.org/content/biorxiv/early/2020/05/08/2020.05.06.077982/DC1/embed/media-1.pdf).                             


## Data Preprocessing

  Deta preprocessing of GRAPPLE needs one file for the disease (outcome file), and two files for each risk factor, where one is for SNP selection (selection file) and the other is for estimating the marginal effects of the SNPs on this risk factor (exposure file). We allow overlapping individuals for the GWAS data in the outcome and exposure files, which there should be as few overlapping individuals as possible between the selection files an other files to avoid SNP selection bias. Sometimes it can be hard to find such GWAS datasets for SNP selection. An option is to use GWAS data from other ancestries.


 As an example, one can download the summary statistics files for BMI and T2D from [datasets](https://www.dropbox.com/sh/vv6pz09cknyz9ca/AAAV_WWLsJmI2LZwL1da45q0a?dl=0), and 
 also download the reference panel and unzip the ".tgz" file in the current folder. 
 Then, we can start by selecting SNPs from these files as genetic instruments.


```
library(GRAPPLE)
sel.file <- "BMI-ukb.csv"
exp.file <- "BMI-giant17eu.csv"
out.file <- "T2D-diagram12-M.csv"
plink_refdat <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"
```

Extracting the independent SNPs can take a few minutes depending on the size of the GWAS summary statistics files.



```
data.list <- getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres = 0.01, plink_exe = 'plink')
```

This function extracts all independent SNPs whose selection p-values do not exceed 0.01. It returns three elements. One is 'data' which is a data frame of the summary statistics of the selected SNPs, the other is 'marker.data' which is a data frame for all candidate marker SNPs. The difference between 'data' and 'marker.data' is that we use a more stringent r2 (0.001) in LD clumping for the selection of genetic instruments to guarantee independence than in selecting candidate marker SNPs (r2 = 0.05). The third element is the estimated correlation matrix for the GWAS cohorts. 

We use PLINK to perform LD clumping. The default name/path of the plink command is 'plink' (passed through the argument \code{plink_exe}). For users with Linux / Windows systems, this command would not work and one need to specify the path of the exe file, like "./plink" depending on where they install plink. If R fails to run the plink command, there will be an error stating that the clumped file is not found.

## Basic Usage

If the number of risk factors is 1, GRAPPLE can detect the number of pleiotropic pathways by finding the number of modes in the robustified profile likelihood, given a p-value selection threshold. 

### Estimation when no pleiotropic pathways are detected

We can try the above example. 
```
## Here we take the p-value threshold be 1e-4, but one can try a series of p-values and see if the result is consistent
diagnosis <- findModes(data.list$data, p.thres = 1e-4)

diagnosis$modes
diagnosis$p
```
There is only one mode in the likelihood, so we do not find any evidence of multiple pleiotropic pathways. It gives us some 
guarantee to use MR-RAPs to estimate the causal effect of the risk factor. 
```
result <- grappleRobustEst(data.list$data, p.thres = 1e-4)
```
This function returns our estimate results.

### Estimation when multiple pleiotropic pathways are detected

Here is another example where we can detect multiple modes

```
library(GRAPPLE)
sel.file <- "CRP-Prins17.csv"
exp.file <- "CRP-Dehghan11.csv"
out.file <- "CAD-Nelson17.csv"
data.list <- getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres = 0.01)
diagnosis <- findModes(data.list$data, p.thres = 1e-4, marker.data = data.list$marker.data)
diagnosis$p
```

We find three modes in the profile likelihood. Checking the marker SNPs, marker genes and mapped GWAS trait, we can find that 
LDL-C can be a confounding trait. We can then run GRAPPLE adjusting for the LDL-C effects


```
library(GRAPPLE)
sel.file <-c(sel.file, "LDL-gera18.csv")
exp.file <- c(exp.file, "LDL-glgc13.csv")
out.file <- "CAD-Nelson17.csv"
data.list <- getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres = 0.01)
result <- grappleRobustEst(data.list$data, p.thres = 1e-4)
```

After adjusting for LDL-C, we can not find any evidence that there is a causal effect of CRP on CAD.



