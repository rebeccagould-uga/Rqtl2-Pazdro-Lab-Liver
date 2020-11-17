# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - LOG TRANSFORMATION AND DATA PREP

#Make a folder under users (for me, my user is "becca") and title it based on your project. For mine, it's R01_GSH_DO_mapping. Then make a data, results, scripts, and docs folder.

setwd("/users/becca/R01_GSH_DO_mapping_Liver/data")

#I can use "~" instead of "/users/becca/" everytime as it represents my home base

#load the command line tools - see https://github.com/Rdatatable/data.table/wiki/Installation for more information - must do every time you open up the Rproject!
library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(RSQLite)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (RSQLite)
library (writexl)

#For all plots, use these parameters to keep everything consistent. You can adjust them accordingly.
#par(mar=c(4.1, 4.1, 2.6, 2.6))

####################################################
## Helpful Sites
####################################################

#https://kbroman.org/qtl2/assets/vignettes/user_guide.html#QTL_analysis_in_Diversity_Outbred_mice
#https://smcclatchy.github.io/mapping/ 


####################################################
## Read in the control file (gm.json)
####################################################

#Load in the control file to tell R/qtl2 to read the data and title it according to your project. For mine, it's R01_GSH_DO_QTLdata. 
#For this to work, all of the files in the control file need to be in the folder. Ex: if calling for the "genoprobs" file, it needs to actually be in the folder.
  R01_GSH_DO_QTLdata <- read_cross2(file = "~/R01_GSH_DO_mapping_Liver/data/R01_GSH_DO_control.json")

####################################################
## Genotype probabilities and allele probabilities - provided by Belinda and Vivek
####################################################

#read in the genoprobs file that is sorted by chromosomes in numerical order - the 8state.rds is the allele probabilities, the 32state.rds is the genotype probabilities
#^this is actually the ALlELE probabilities, but for simplicity, we will call it "probs"
  probs <- readRDS("~/R01_GSH_DO_mapping_Liver/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")

  nrow(data.frame(R01_GSH_DO_QTLdata$gmap[1]))
  #should be 10415
  dim(probs[[1]])
  #should be 347 individuals, 8 alleles, 10415 markers


####################################################
## Variant files
####################################################

#Will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice. Make sure they are the most updated versions!
  query_variants <- create_variant_query_func("~/R01_GSH_DO_mapping_Liver/data/cc_variants.sqlite")
  query_genes_mgi <- create_gene_query_func("~/R01_GSH_DO_mapping_Liver/data/mouse_genes_mgi.sqlite")
  query_genes <- create_gene_query_func("~/R01_GSH_DO_mapping_Liver/data/mouse_genes.sqlite")

####################################################
## Calculating kinship
####################################################

#calculate the kinship loco
#you can increase the cores amount if you have more cores in your computer. For mine, I have 18 cores available so to speed it up, I'll use 10 of them.
  kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE, cores = 10)

#Create the r plot of the kinship matrix
#pdf('output/rplot_kinship.pdf')
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  image(1:nrow(kinship_loco[[1]]), 1:ncol(kinship_loco[[1]]), kinship_loco[[1]][,ncol(kinship_loco[[1]]):1], xlab = "Samples", 
        ylab = "Samples", yaxt = "n", main = "Kinship between samples", 
        breaks = 0:100/100, col = heat.colors(length(0:100) - 1))


####################################################
## Editing the phenotype file to make it R/qtl2-friendly
####################################################
#to have the phenotype file for reference - can be used when plotting the data to see if it needs to be transformed
  pheno <- read.csv(file = "~/R01_GSH_DO_mapping_Liver/data/R01_GSH_DO_pheno_covar.csv", header = TRUE)

#make row names the ID of each sample
  rownames(pheno) <- pheno$id
#checking pheno file
  pheno[1:10,]

#change sexes to numeric variables
  pheno$sex[pheno$sex == "M"] <- 1
  pheno$sex[pheno$sex == "F"] <- 0

#both added covariates must be numeric, not characters 
  pheno$sex <- as.numeric(pheno$sex)
  pheno$generation <- as.numeric(pheno$generation)  

  #when I ran this, it still gives a faulty covariate file
  #sex <- (pheno$sex == "1")*1
  #names(sex) <- rownames(pheno$sex)
  #sex
  #gen <- (pheno$generation)
  #names(gen) <- rownames(pheno$generation)
  #gen

#check pheno file
  pheno[1:10,]
  str(pheno)


####################################################
## checking if data should be transformed
####################################################

#gives you the names of each phenotype
  R01_GSH_DO_QTLdata[["pheno"]]


#Log transformations of each phenotype
pheno$logLiverGSH = log(pheno$Liver_GSH)
pheno$logLiverGSSG = log(pheno$Liver_Adj_GSSG)
pheno$logLiverTotalGSH = log(pheno$Liver_Adj_Total_GSH)
pheno$logLiverGSH_GSSGRatio = log(pheno$Liver_Adj_GSH_GSSG_Ratio)
pheno$logLiverNADH = log(pheno$Liver_NADH)
pheno$logLiverNADP = log(pheno$Liver_NADP)
pheno$logLiverNADPH = log(pheno$Liver_NADPH)
pheno$logLiverNADP_NADPHRatio = log(pheno$Liver_NADP_NADPH_Ratio)
pheno$logAST = log(pheno$AST)
pheno$logALT = log(pheno$ALT)
pheno$logBUN = log(pheno$BUN)

#####Plot the transformations  

#setting the parameters for the plots
par(mar=c(4.1, 4.1, 2.6, 2.6))

#For Liver GSH
  boxplot(pheno$Liver_GSH, main = "Liver GSH Box Plot")
  boxplot(pheno$Liver_GSH~pheno$generation, main = "Liver GSH Box Plot - by generation")
  boxplot(pheno$logLiverGSH, main = "Log Liver GSH Box Plot")
  boxplot(pheno$logLiverGSH~pheno$generation, main = "Log Liver GSH Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverGSH, main = "Normal QQ Plot - Log Liver GSH")

#For Liver GSSG
  boxplot(pheno$Liver_Adj_GSSG, main = "Liver GSSG Box Plot")
  boxplot(pheno$Liver_Adj_GSSG~pheno$generation, main = "Liver GSSG Box Plot - by generation")
  boxplot(pheno$logLiverGSSG, main = "Log Liver GSSG Box Plot")
  boxplot(pheno$logLiverGSSG~pheno$generation, main = "Log Liver GSSG Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverGSSG, main = "Normal QQ Plot - Log Liver GSSG")

#For Liver Total GSH
  boxplot(pheno$Liver_Adj_Total_GSH, main = "Liver Total GSH Box Plot")
  boxplot(pheno$Liver_Adj_Total_GSH~pheno$generation, main = "Liver Total GSH Box Plot - by generation")
  boxplot(pheno$logLiverTotalGSH, main = "Log Liver Total GSH Box Plot")
  boxplot(pheno$logLiverTotalGSH~pheno$generation, main = "Log Liver Total GSH Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverTotalGSH, main = "Normal QQ Plot - Log Liver Total GSH")

#For Liver GSH/GSSG Ratio
  boxplot(pheno$Liver_Adj_GSH_GSSG_Ratio, main = "Liver GSH/GSSG Box Plot")
  boxplot(pheno$Liver_Adj_GSH_GSSG_Ratio~pheno$generation, main = "Liver GSH/GSSG Box Plot - by generation")
  boxplot(pheno$logLiverGSH_GSSGRatio, main = "Log Liver GSH/GSSG Box Plot")
  boxplot(pheno$logLiverGSH_GSSGRatio~pheno$generation, main = "Log Liver GSH/GSSG Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverGSH_GSSGRatio, main = "Normal QQ Plot - Log Liver GSH/GSSG")

#For Liver NADH
  boxplot(pheno$Liver_NADH, main = "Liver NADH Box Plot")
  boxplot(pheno$Liver_NADH~pheno$generation, main = "Liver NADH Box Plot - by generation")
  boxplot(pheno$logLiverNADH, main = "Log Liver NADH Box Plot")
  boxplot(pheno$logLiverNADH~pheno$generation, main = "Log Liver NADH Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverNADH, main = "Normal QQ Plot - Log Liver NADH")

#For Liver NADP
  boxplot(pheno$Liver_NADP, main = "Liver NADP Box Plot")
  boxplot(pheno$Liver_NADP~pheno$generation, main = "Liver NADP Box Plot - by generation")
  boxplot(pheno$logLiverNADP, main = "Log Liver NADP Box Plot")
  boxplot(pheno$logLiverNADP~pheno$generation, main = "Log Liver NADP Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverNADP, main = "Normal QQ Plot - Log Liver NADP")

#For Liver NADPH
  boxplot(pheno$Liver_NADPH, main = "Liver NADPH Box Plot")
  boxplot(pheno$Liver_NADPH~pheno$generation, main = "Liver NADPH Box Plot - by generation")
  boxplot(pheno$logLiverNADPH, main = "Log Liver NADPH Box Plot")
  boxplot(pheno$logLiverNADPH~pheno$generation, main = "Log Liver NADPH Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverNADPH, main = "Normal QQ Plot - Log Liver NADPH")

#For Liver NADP/NADPH Ratio
  boxplot(pheno$Liver_NADP_NADPH_Ratio, main = "Liver NADP/NADPH Box Plot")
  boxplot(pheno$Liver_NADP_NADPH_Ratio~pheno$generation, main = "Liver NADP/NADPH Box Plot - by generation")
  boxplot(pheno$logLiverNADP_NADPHRatio, main = "Log Liver NADP/NADPH Box Plot")
  boxplot(pheno$logLiverNADP_NADPHRatio~pheno$generation, main = "Log Liver NADP/NADPH Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logLiverNADP_NADPH, main = "Normal QQ Plot - Log Liver NADP/NADPH")

#For AST
  boxplot(pheno$AST, main = "AST Box Plot")
  boxplot(pheno$AST~pheno$generation, main = "AST Box Plot - by generation")
  boxplot(pheno$logAST, main = "Log AST Box Plot")
  boxplot(pheno$logAST~pheno$generation, main = "Log AST Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logAST, main = "Normal QQ Plot - Log AST")

#For ALT
  boxplot(pheno$ALT, main = "ALT Box Plot")
  boxplot(pheno$ALT~pheno$generation, main = "ALT Box Plot - by generation")
  boxplot(pheno$logALT, main = "Log ALT Box Plot")
  boxplot(pheno$logALT~pheno$generation, main = "Log ALT Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logALT, main = "Normal QQ Plot - Log ALT")

#For BUN
  boxplot(pheno$BUN, main = "BUN Box Plot")
  boxplot(pheno$BUN~pheno$generation, main = "BUN Box Plot - by generation")
  boxplot(pheno$logBUN, main = "Log BUN Box Plot")
  boxplot(pheno$logBUN~pheno$generation, main = "Log BUN Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$logBUN, main = "Normal QQ Plot - Log BUN")


####################################################
## add covariates
####################################################

#because of how large my cohorts were, it is important to account for generation + sex

#adding sex and generation as covariates
sexgen = model.matrix(~ sex + generation, data = pheno)[,-1]

  #doing it this way does not transfer the row names correctly
  #sexgen = model.matrix(~ sex + gen)[,-1]


#For heritability calculation, need a linear mixed model
#####make kinship function using linear mixed model, not loco
#####default type of kinship is "overall" aka "linear mixed model" -- did not need to specify a type
  
kinship_lmm <- calc_kinship(probs = probs, use_allele_probs = TRUE, cores = 10)

#adding sex as covariate to compare to sexgen
sex = model.matrix(~ sex, data = pheno)[,-1] 

####################################################
## exporting data
####################################################

#set working directory
write_xlsx(list("GSH chr14" = LiverGSH_Genes_MGI_chr14, 
                "GSH chr18" = LiverGSH_Genes_MGI_chr18, 
                "GSSG chr1" = LiverGSSG_Genes_MGI_chr1, 
                "GSSG chr2" = LiverGSSG_Genes_MGI_chr2, 
                "GSSG chr18" = LiverGSSG_Genes_MGI_chr18, 
                "Total GSH chr14" = LiverTotalGSH_Genes_MGI_chr14, 
                "Total GSH chr18" = LiverTotalGSH_Genes_MGI_chr18, 
                "GSH GSSG Ratio chr16" = LiverGSH_GSSGRatio_Genes_MGI_chr16, 
                "GSH GSSG Ratio chr11" = LiverGSH_GSSGRatio_Genes_MGI_chr11),
                "GlutathioneGenesMGI.xlsx")

#set working directory
write_xlsx(list("NADH chr14" = LiverNADH_Genes_MGI_chr14,
                "NADP chr3" = LiverNADP_Genes_MGI_chr3,
                "NADP chr8" = LiverNADP_Genes_MGI_chr8,
                "NADPH chr12" = LiverNADPH_Genes_MGI_chr12,
                "NADPH chr3" = LiverNADPH_Genes_MGI_chr3,
                "NADP NADPH Ratio chr12a" = LiverNADP_NADPHRatio_Genes_MGI_chr12a,
                "NADP NADPH Ratio chr12b" = LiverNADP_NADPHRatio_Genes_MGI_chr12b,
                "NADP NADPH Ratio chr6" = LiverNADP_NADPHRatio_Genes_MGI_chr6,
                "NADP NADPH Ratio chr3" = LiverNADP_NADPHRatio_Genes_MGI_chr3),
                "NADSystemsGenesMGI.xlsx")

#set working directory
write_xlsx(list("AST chr16" = AST_Genes_MGI_chr16),
           "BloodValuesGenesMGI - log sexgen.xlsx")
