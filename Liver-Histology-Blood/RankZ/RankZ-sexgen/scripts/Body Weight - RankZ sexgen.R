# R01 BodyWeight DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Body Weight

#Load in Liver-Histology-Blood-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_BodyWeight_DO_mapping_Liver/data")

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


####################################################
## Plot Genome Scans with Permutation Tests
####################################################

qtlscan_BodyWeight <- scan1(genoprobs = probs, pheno = pheno["zBodyWeight"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_BodyWeight <- scan1perm(genoprobs = probs, pheno = pheno["zBodyWeight"], addcovar = sexgen, n_perm = 1000, cores=10)


#set working directory
pdf(file = "Body Weight QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_BodyWeight = summary(perm_BodyWeight, alpha = c(0.2, 0.1, 0.05))

plot_scan1(x = qtlscan_BodyWeight, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Body Weight", ylim = c(0,11))
abline(h = threshold_BodyWeight, col = c("purple", "red", "blue"), lwd = 2)
plot_scan1(x = qtlscan_BodyWeight, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Body Weight", ylim = c(0,11))
abline(h = threshold_BodyWeight, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")

#using gmap (cM)
find_peaks(scan1_output = qtlscan_BodyWeight, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_BodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
gmap_peaksBodyWeight <- find_peaks(scan1_output = qtlscan_BodyWeight, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_BodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(gmap_peaksBodyWeight)

#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_BodyWeight, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_BodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
peaksBodyWeight <- find_peaks(scan1_output = qtlscan_BodyWeight, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_BodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(peaksBodyWeight)


write_xlsx(list("BodyWeight gmap (cM)" = gmap_peaksBodyWeight,
                "BodyWeight pmap (Mbp)" = peaksBodyWeight),
           "BodyWeight Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Body Weight GWAS - RankZ sexgen.pdf")
out_gwas_BodyWeight <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zBodyWeight"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_BodyWeight$lod, out_gwas_BodyWeight$snpinfo, altcol="green4", gap=0, main = "BodyWeight GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_BodyWeight_sex <- est_herit(pheno["zBodyWeight"], kinship_lmm, sex, cores = 2)
herit_BodyWeight_sexgen <- est_herit(pheno["zBodyWeight"], kinship_lmm, sexgen, cores = 2)

