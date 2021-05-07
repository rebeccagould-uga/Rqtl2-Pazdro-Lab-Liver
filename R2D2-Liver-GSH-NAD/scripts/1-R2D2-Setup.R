# R2D2 Project
# Created by Becca Gould
# Updated May 2021

#### Phenotypes ####
# Chr 2 (WSB): GSH, GSSG, Total GSH, Eh, NADP, NADP/NADPH
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

#control file
  R01_GSH_DO_QTLdata <- read_cross2(file = "~/Rqtl2-Glutathione-Genetics/data/R01_GSH_DO_control.json")

#genoprobs/alleleprobs
  probs <- readRDS("~/Rqtl2-Glutathione-Genetics/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")

#variant files
  query_variants <- create_variant_query_func("~/Rqtl2-Glutathione-Genetics/data/cc_variants.sqlite")
  query_genes_mgi <- create_gene_query_func("~/Rqtl2-Glutathione-Genetics/data/mouse_genes_mgi.sqlite")
  query_genes <- create_gene_query_func("~/Rqtl2-Glutathione-Genetics/data/mouse_genes.sqlite")

#kinship calculation
  kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE, cores = 2)

#pheno file
  pheno <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/data/R01_GSH_DO_pheno_covar.csv", header = TRUE)
#make row names the ID of each sample
  rownames(pheno) <- pheno$id
#change sexes to numeric variables
  pheno$sex[pheno$sex == "M"] <- 1
  pheno$sex[pheno$sex == "F"] <- 0
#both added covariates must be numeric, not characters 
  pheno$sex <- as.numeric(pheno$sex)
  pheno$generation <- as.numeric(pheno$generation)  
#check pheno file
  pheno[1:10,]
  str(pheno)

#rankZ transformation function
  rankZ <- function(x) {x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))}

#transform phenotypes of interest
  pheno$zLiverGSH = rankZ(pheno$Liver_GSH)
  pheno$zLiverGSSG = rankZ(pheno$Liver_Adj_GSSG)
  pheno$zLiverTotalGSH = rankZ(pheno$Liver_Adj_Total_GSH)
  pheno$zLiverEh = rankZ(pheno$Liver_Adj_Redox_Potential_GSSG_2GSH)
  pheno$zLiverNADP = rankZ(pheno$Liver_NADP)
  pheno$zLiverNADPH = rankZ(pheno$Liver_NADPH)
  pheno$zLiverNADP_NADPHRatio = rankZ(pheno$Liver_NADP_NADPH_Ratio)

#plot transformations
  boxplot(pheno$zLiverGSH, main = "Rank Z Liver GSH Box Plot")
  boxplot(pheno$zLiverGSH~pheno$generation, main = "Rank Z Liver GSH Box Plot - by generation")

#normality plots
  qqnorm(pheno$zLiverGSH, main = "Normal QQ Plot - Rank Z Liver GSH")
  qqnorm(pheno$zLiverGSSG, main = "Normal QQ Plot - Rank Z Liver GSSG")
  qqnorm(pheno$zLiverTotalGSH, main = "Normal QQ Plot - Rank Z Liver Total GSH")
  qqnorm(pheno$zLiverEh, main = "Normal QQ Plot - Rank Z Liver Eh")
  qqnorm(pheno$zLiverNADP, main = "Normal QQ Plot - Rank Z Liver NADP")
  qqnorm(pheno$zLiverNADPH, main = "Normal QQ Plot - Rank Z Liver NADPH")
  qqnorm(pheno$zLiverNADP_NADPHRatio, main = "Normal QQ Plot - Rank Z Liver NADP/NADPH")

#covariates
  sexgen = model.matrix(~ sex + generation, data = pheno)[,-1]
  sex = model.matrix(~ sex, data = pheno)[,-1] 
  



