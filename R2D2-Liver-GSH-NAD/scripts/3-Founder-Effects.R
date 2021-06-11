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

pdf(file = "Founder-Effects.pdf")

par(mar=c(4.1, 4.1, 2.6, 2.6))


########################################################################################################
## 
## Liver GSH (WSB - Chromosome 2)
##
########################################################################################################

# founder allele effects 
  chr = 2
  
  #sex
    #BLUPs
    coef_blup_LiverGSH_chr2_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sex, main = "Liver GSH BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sex, main = "Liver GSH BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverGSH_chr2_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver GSH default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver GSH default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
  #sexgen
    #BLUPs
    coef_blup_LiverGSH_chr2_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverGSH_chr2_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  

########################################################################################################
## 
## Liver GSSG (WSB - Chromosome 2)
##
########################################################################################################
    
    # founder allele effects 
    chr = 2
    
    #sex
    #BLUPs
    coef_blup_LiverGSSG_chr2_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverGSSG_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG_sex, main = "Liver GSSG BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverGSSG_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG_sex, main = "Liver GSSG BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverGSSG_chr2_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverGSSG_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver GSSG default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverGSSG_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver GSSG default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
    #sexgen
    #BLUPs
    coef_blup_LiverGSSG_chr2_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverGSSG_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG_sexgen, main = "Liver GSSG BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverGSSG_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG_sexgen, main = "Liver GSSG BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverGSSG_chr2_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverGSSG_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG_sexgen, main = "Liver GSSG default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverGSSG_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG_sexgen, main = "Liver GSSG default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    
########################################################################################################
## 
## Liver Total GSH (WSB - Chromosome 2)
##
########################################################################################################
    
    # founder allele effects 
    chr = 2
    
    #sex
    #BLUPs
    coef_blup_LiverTotalGSH_chr2_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverTotalGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH_sex, main = "Liver TotalGSH BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverTotalGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH_sex, main = "Liver TotalGSH BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverTotalGSH_chr2_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverTotalGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver TotalGSH default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverTotalGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver TotalGSH default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
    #sexgen
    #BLUPs
    coef_blup_LiverTotalGSH_chr2_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverTotalGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH_sexgen, main = "Liver TotalGSH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverTotalGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH_sexgen, main = "Liver TotalGSH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverTotalGSH_chr2_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverTotalGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH_sexgen, main = "Liver TotalGSH default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverTotalGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH_sexgen, main = "Liver TotalGSH default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    

########################################################################################################
## 
## Liver Eh (WSB - Chromosome 2)
##
########################################################################################################
    
    # founder allele effects 
    chr = 2
    
    #sex
    #BLUPs
    coef_blup_LiverEh_chr2_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverEh"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverEh_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverEh_sex, main = "Liver Eh BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverEh_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverEh_sex, main = "Liver Eh BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverEh_chr2_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverEh"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverEh_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver Eh default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverEh_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver Eh default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
    #sexgen
    #BLUPs
    coef_blup_LiverEh_chr2_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverEh"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverEh_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverEh_sexgen, main = "Liver Eh BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverEh_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverEh_sexgen, main = "Liver Eh BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverEh_chr2_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverEh"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverEh_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverEh_sexgen, main = "Liver Eh default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverEh_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverEh_sexgen, main = "Liver Eh default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    
########################################################################################################
## 
## Liver NADP (WSB - Chromosome 2)
##
########################################################################################################
    
    # founder allele effects 
    chr = 2
    
    #sex
    #BLUPs
    coef_blup_LiverNADP_chr2_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADP_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sex, main = "Liver NADP BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverNADP_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sex, main = "Liver NADP BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADP_chr2_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverNADP_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADP default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverNADP_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADP default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
    #sexgen
    #BLUPs
    coef_blup_LiverNADP_chr2_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADP_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sexgen, main = "Liver NADP BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverNADP_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sexgen, main = "Liver NADP BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADP_chr2_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverNADP_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sexgen, main = "Liver NADP default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverNADP_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sexgen, main = "Liver NADP default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    

########################################################################################################
## 
## Liver NADP/NADPH (WSB - Chromosome 2)
##
########################################################################################################
    
    # founder allele effects 
    chr = 2
    
    #sex
    #BLUPs
    coef_blup_LiverNADP_NADPHRatio_chr2_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, main = "Liver NADP_NADPHRatio BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, main = "Liver NADP_NADPHRatio BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADP_NADPHRatio_chr2_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADP_NADPHRatio default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADP_NADPHRatio default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
    #sexgen
    #BLUPs
    coef_blup_LiverNADP_NADPHRatio_chr2_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADP_NADPHRatio_chr2_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    

########################################################################################################
## 
## Liver NADPH (WSB - Chromosome 14)
##
########################################################################################################
    
    # founder allele effects 
    chr = 14
    
    #sex
    #BLUPs
    coef_blup_LiverNADPH_chr14_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADPH_chr14_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sex, main = "Liver NADPH BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_blup_LiverNADPH_chr14_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sex, main = "Liver NADPH BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADPH_chr14_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverNADPH_chr14_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADPH default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_def_LiverNADPH_chr14_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADPH default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
    #sexgen
    #BLUPs
    coef_blup_LiverNADPH_chr14_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADPH_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sexgen, main = "Liver NADPH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_blup_LiverNADPH_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sexgen, main = "Liver NADPH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADPH_chr14_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverNADPH_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sexgen, main = "Liver NADPH default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_def_LiverNADPH_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sexgen, main = "Liver NADPH default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

    
########################################################################################################
## 
## Liver NADP/NADPH (WSB - Chromosome 14)
##
########################################################################################################
    
    # founder allele effects 
    chr = 14
    
    #sex
    #BLUPs
    coef_blup_LiverNADP_NADPHRatio_chr14_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, main = "Liver NADP_NADPHRatio BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, main = "Liver NADP_NADPHRatio BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADP_NADPHRatio_chr14_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr14_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADP_NADPHRatio default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr14_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver NADP_NADPHRatio default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
    #sexgen
    #BLUPs
    coef_blup_LiverNADP_NADPHRatio_chr14_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverNADP_NADPHRatio_chr14_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(35,50)
    plot_coefCC(x = coef_def_LiverNADP_NADPHRatio_chr14_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, main = "Liver NADP_NADPHRatio default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    
    
    
dev.off()
    
    
  
  

