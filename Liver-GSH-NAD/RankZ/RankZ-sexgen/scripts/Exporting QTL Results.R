# R01 GSH DO Mapping Code 
# Updated September 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - Exporting QTL Results

#Load in Liver-GSH-NAD-RankZ-SexGen.Rdata


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
write_xlsx(list(  "GSH chr14" = LiverGSH_Genes_MGI_chr14, 
                  "Total GSH chr14" = LiverTotalGSH_Genes_MGI_chr14, 
                  "GSH GSSG Ratio chr16" = LiverGSH_GSSGRatio_Genes_MGI_chr16), 
              "GlutathioneGenesMGI - rankZ sexgen.xlsx")
  
#set working directory
write_xlsx(list(  "NADP chr3" = LiverNADP_Genes_MGI_chr3,
                  "NADP chr8" = LiverNADP_Genes_MGI_chr8,
                  "NADPH chr12" = LiverNADPH_Genes_MGI_chr12,
                  "NADP NADPH Ratio chr12" = LiverNADP_NADPHRatio_Genes_MGI_chr12),
             "NADSystemsGenesMGI - rankZ sexgen.xlsx")
  
  