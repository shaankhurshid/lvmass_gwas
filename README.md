# lvmass_gwas
These scripts were used to perform the analyses described in the LV Mass GWAS manuscript by Khurshid et al. (https://www.medrxiv.org/content/10.1101/2022.01.09.22268962v1)

All scripts written by Shaan Khurshid 2022, excepting the GWAS bash scripts which were written by James Pirruccello (https://github.com/carbocation), see header information.

# Overview
* "prep_gwas" - contains scripts used to prepare phenotype files for GWAS using BOLT-LMM
* "prs" - scripts used to process GWAS results and create polygenic risk scores using pruning-and-thresholding in plink v1.90b
* "phenotype_outcomes" - scripts used to perform longitudinal analyses using LVM and LVM PRS in UK Biobank and MGB datasets
* "plots" - visualizations
* "baseline characteristics" - scripts used to summarize baseline characteristics of the samples
* "gwas" - scripts to perform GWAS using BOLT-LMM on GCP using dsub
