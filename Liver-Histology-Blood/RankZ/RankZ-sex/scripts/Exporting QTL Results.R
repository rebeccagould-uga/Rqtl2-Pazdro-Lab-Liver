# R01 GSH DO Mapping Code 
# Updated November 2020
# Becca Gould 


#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Exporting QTL Results

#Load in Liver-Histology-Blood-RankZ-Sex.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


#load the command line tools 
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


#set working directory to store the plots

#tells you all of the qtlscans that you have
ls(pattern = "qtl")


####################################################
## Grab the QTL from all of the qtl scans 
####################################################
  
## I just set threshold to 6 (tells you all of the important qtl peaks with a LOD score > 6)
## map is the qtl2 map you want to use (gmap or pmap)
  qtl_gmap <- find_peaks(scans, map = R01_GSH_DO_QTLdata$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #shows you the first 6 results
  head(qtl_gmap)
  #shows you all of the results - check that these match your records
  qtl_gmap
  

####################################################
## exporting data
####################################################
  
  #set working directory
  write_xlsx(list("AST chr16" = AST_Genes_MGI_chr16),
             "BloodValuesGenesMGI - RankZ sex.xlsx")
  
  
  

