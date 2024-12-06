# A-TWAS

## Overview
A-TWAS is an omnibus tool integrating multiple imputation models with different structural assumptions on the genotype structure. It assuming individual level eQTL data and summary level GWAS data.
## Model Implementation and Data Preparation
The command for running the model is:
```R
source("my_path/A_TWAS.R")
model <- A_TWAS(genotype, expression, GWAS_summary, CHR, LD, option = c("BLasso", "Horseshoe", "Horseshoe+"), extra_eQTL_weights = NA)
```
The detailed description of the data that need to be prepared for implementing the model is stated as below:
1. genotype (genotype information): the genotype information for the imputation model in stage I TWAS should be a dataframe where each row represents a sample and each column represents a SNP; the names of the columns should be the basic position for the corresponding SNPs. Remind that the SNPs of the dataframe should be in the same chromosome.
2. expression (gene expression): the gene expression information should be saved as a vector; the order of sample in the vector (i.e. the order of the element) should be aligned with the order of sample in the dataframe of genotype information (i.e. the order of the row).
3. GWAS_summary (GWAS summary data): a dataframe mandatorily containing at least the column of zscore of SNPs, basic position of the SNPs and the chromosome the SNPs belonging to; the name of these three columns should be c("zscore", "BP", "CHR"). Reminder that the major and minor allele for SNPs of GWAS summary dataset should be consistent with the one for the genotype information.
4. CHR: the chromosome that the SNPs lies on.
5. LD (reference LD panel): a dataframe measures the association between the SNPs. The column names of the dataframe should be the basic position of the corresponding SNPs.
6. option: This parameter is used for choosing which imputation models are used. The choice includes "BLasso", "Horseshoe" and "Horseshoe_plus". By default, all three methods are used in model.
7. extra_eQTL_weights: this term is optional. If needs to combine the outcome from other imputation models, input a dataframe storing the basic position and the weights of the SNPs. The first column of the dataframe indicating the basic position for the SNPs that the weights corresponding to; the remaining column of the dataframe indicating the weights from each of the external model. The names for the first column should be "BP" and for the remaining column should be the names of the external models.

## Dependency
A-TWAS is built based on R package "bayesreg" and "ACAT", we sincerely appreciate the authors of these packages.
