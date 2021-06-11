# R2D2 Project
# Created by Becca Gould
# Updated May 2021

#### Phenotypes ####
# Chr 2 (WSB): GSH, GSSG, Total GSH (omitted), Eh, NADP, NADP/NADPH
# Chr 14 (NOD): NADPH, NADP/NADPH

library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (pander)
library (writexl)
library (RSQLite)


#Load in R2D2-Liver-GSH-NAD-RankZ.Rdata
#Run 1-R2D2-Setup.R prior to this script




########################################################################################################
##
##
##
########################################################################################################

  # pull genotype frequencies 
  pull <- pull_genoprobpos(probs, map = R01_GSH_DO_QTLdata$gmap, 2, 47)
  head(pull)






########################################################################################################
## Testing regression matrix - is it ill-conditioned? 
## why the covariates affected the scan
########################################################################################################

  # instructions from Gary Churchill
  
  # kappa() computes the condition number of a matrix
  # Construct a matrix with the genoprobs at the maker locus.  It should have 8 columns A-H.  Compute kappa.
  # Repeat this but drop column H and replace it with a column of all ones – kappa should be the same.
  # Now add a column with 1 for females and 0 for males…compute kappa.
  # Next add several column with dummy variables for cohort…compute kappa.
  # Just to get a feel for things, swap in genotypes from loci where there is no problem.
  # You can compute the Sex and Cohort columns using the model.matrix() function,
  # e.g. model.matrix( ~Sex, data=my.dataframe) or model.matrix(~Sex+Cohort, data=my.dataframe)
  
  kappa()






########################################################################################################
## 
## Fit1 model
##
########################################################################################################

  #sex
  #BLUPs
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sex, model = "normal", se = TRUE, blup = TRUE, cores = 2)
  #default
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sex, model = "normal", se = TRUE, blup = FALSE, cores = 2)
  #sexgen  
  #BLUPs
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, model = "normal", se = TRUE, blup = TRUE, cores = 2)
  #default
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, model = "normal", se = TRUE, blup = FALSE, cores = 2)




