# R01 GSH DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - RankZ TRANSFORMATION AND DATA PREP

#setwd

#load the command line tools - see https://github.com/Rdatatable/data.table/wiki/Installation for more information - must do every time you open up the Rproject!
library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (writexl)
library (RSQLite)

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
  R01_GSH_DO_QTLdata <- read_cross2(file = "~/Rqtl2-Glutathione-Genetics/data/R01_GSH_DO_control.json")

####################################################
## Genotype probabilities and allele probabilities - provided by Belinda and Vivek
####################################################

#read in the genoprobs file that is sorted by chromosomes in numerical order - the 8state.rds is the allele probabilities, the 32state.rds is the genotype probabilities
#^this is actually the ALlELE probabilities, but for simplicity, we will call it "probs"
  probs <- readRDS("~/Rqtl2-Glutathione-Genetics/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")

  nrow(data.frame(R01_GSH_DO_QTLdata$gmap[1]))
  #should be 10415
  dim(probs[[1]])
  #should be 347 individuals, 8 alleles, 10415 markers


####################################################
## Variant files
####################################################

#Will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice. Make sure they are the most updated versions!
  query_variants <- create_variant_query_func("~/Rqtl2-Glutathione-Genetics/data/cc_variants.sqlite")
  query_genes_mgi <- create_gene_query_func("~/Rqtl2-Glutathione-Genetics/data/mouse_genes_mgi.sqlite")
  query_genes <- create_gene_query_func("~/Rqtl2-Glutathione-Genetics/data/mouse_genes.sqlite")

####################################################
## Calculating kinship
####################################################

#calculate the kinship loco
#you can increase the cores amount if you have more cores in your computer. For mine, I have 18 cores available so to speed it up, I'll use 10 of them.
  kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE, cores = 2)

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
  pheno <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/data/R01_GSH_DO_pheno_covar.csv", header = TRUE)

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

#The rankZ function is a nonparametric method that replaces the values of a variable with their rank in 
#ascending order (e.g. values 0.2, 1.3, and 4.3 are replaced by 1, 2 and 3). This effectively forces the data into a normal distribution and eliminates outliers.
  rankZ <- function(x) {x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))}
  
#Rank Z transformations of each phenotype
  pheno$zAST = rankZ(pheno$AST)
  pheno$zALT = rankZ(pheno$ALT)
  pheno$zASTALTRatio = rankZ(pheno$AST_ALT_Ratio)
  pheno$zSteatosis = rankZ(pheno$Liver_Steatosis)
  pheno$zBallooning = rankZ(pheno$Liver_Ballooning)
  pheno$zFibrosis = rankZ(pheno$Liver_Fibrosis)
  pheno$zInflammation = rankZ(pheno$Liver_Inflammation)
  
#####Plot the transformations  

  
#set working directory
  
  pdf(file = "Box Plots and QQ Norm Plots - RankZ and Raw.pdf")
  ##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##
  
#For Steatosis
  boxplot(pheno$Steatosis, main = "Liver Steatosis Box Plot")
  boxplot(pheno$Steatosis~pheno$generation, main = "Liver Steatosis Box Plot - by generation")
  boxplot(pheno$zSteatosis, main = "RankZ Liver Steatosis Box Plot")
  boxplot(pheno$zSteatosis~pheno$generation, main = "RankZ Liver Steatosis Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zSteatosis, main = "Normal QQ Plot - RankZ Liver Steatosis") 
  
#For Ballooning
  boxplot(pheno$Ballooning, main = "Liver Ballooning Box Plot")
  boxplot(pheno$Ballooning~pheno$generation, main = "Liver Ballooning Box Plot - by generation")
  boxplot(pheno$zBallooning, main = "RankZ Liver Ballooning Box Plot")
  boxplot(pheno$zBallooning~pheno$generation, main = "RankZ Liver Ballooning - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zBallooning, main = "Normal QQ Plot - RankZ Liver Ballooning") 
  
#For AST
  boxplot(pheno$AST, main = "AST Box Plot")
  boxplot(pheno$AST~pheno$generation, main = "AST Box Plot - by generation")
  boxplot(pheno$zAST, main = "RankZ AST Box Plot")
  boxplot(pheno$zAST~pheno$generation, main = "RankZ AST - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zAST, main = "Normal QQ Plot - RankZ AST") 
  
#For ALT
  boxplot(pheno$ALT, main = "ALT Box Plot")
  boxplot(pheno$ALT~pheno$generation, main = "ALT Box Plot - by generation")
  boxplot(pheno$zALT, main = "RankZ ALT Box Plot")
  boxplot(pheno$zALT~pheno$generation, main = "RankZ ALT - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zALT, main = "Normal QQ Plot - RankZ ALT") 
  
#For AST/ALT Ratio
  boxplot(pheno$AST_ALT_Ratio, main = "AST/ALT Box Plot")
  boxplot(pheno$AST_ALT_Ratio~pheno$generation, main = "AST/ALT Box Plot - by generation")
  boxplot(pheno$zASTALTRatio, main = "RankZ AST/ALT Box Plot")
  boxplot(pheno$zASTALTRatio~pheno$generation, main = "RankZ AST/ALT - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zASTALTRatio, main = "Normal QQ Plot - RankZ AST/ALT") 
  
#For Fibrosis
  boxplot(pheno$Fibrosis, main = "Liver Fibrosis Box Plot")
  boxplot(pheno$Fibrosis~pheno$generation, main = "Liver Fibrosis Box Plot - by generation")
  boxplot(pheno$zFibrosis, main = "RankZ Liver Fibrosis Box Plot")
  boxplot(pheno$zFibrosis~pheno$generation, main = "RankZ Liver ReFibrosis - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zFibrosis, main = "Normal QQ Plot - RankZ Liver Fibrosis") 
  
#For Inflammation
  boxplot(pheno$Inflammation, main = "Liver Inflammation Box Plot")
  boxplot(pheno$Inflammation~pheno$generation, main = "Liver Inflammation Box Plot - by generation")
  boxplot(pheno$zInflammation, main = "RankZ Liver Inflammation Box Plot")
  boxplot(pheno$zInflammation~pheno$generation, main = "RankZ Liver Inflammation - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zInflammation, main = "Normal QQ Plot - RankZ Liver Inflammation") 
  
dev.off()
  
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
  
kinship_lmm <- calc_kinship(probs = probs, use_allele_probs = TRUE, cores = 2)

#adding sex as covariate to compare to sexgen
sex = model.matrix(~ sex, data = pheno)[,-1] 



