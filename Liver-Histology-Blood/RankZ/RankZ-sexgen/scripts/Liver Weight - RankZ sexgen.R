# R01 LiverWeight DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - LiverWeight

#Load in Liver-Histology-Blood-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_LiverWeight_DO_mapping_Liver/data")

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

qtlscan_LiverWeight <- scan1(genoprobs = probs, pheno = pheno["zLiverWeight"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_LiverWeight <- scan1perm(genoprobs = probs, pheno = pheno["zLiverWeight"], addcovar = sexgen, n_perm = 1000, cores=10)

qtlscan_BodyWeight <- scan1(genoprobs = probs, pheno = pheno["zBodyWeight"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_BodyWeight <- scan1perm(genoprobs = probs, pheno = pheno["zBodyWeight"], addcovar = sexgen, n_perm = 1000, cores=10)

qtlscan_LiverWeightBodyWeight <- scan1(genoprobs = probs, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_LiverWeightBodyWeight <- scan1perm(genoprobs = probs, pheno = pheno["zLiverWeightBodyWeight"], addcovar = sexgen, n_perm = 1000, cores=10)


#set working directory
pdf(file = "Liver Weight QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverWeight = summary(perm_LiverWeight, alpha = c(0.2, 0.1, 0.05))

plot_scan1(x = qtlscan_LiverWeight, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Weight", ylim = c(0,11))
abline(h = threshold_LiverWeight, col = c("purple", "red", "blue"), lwd = 2)
plot_scan1(x = qtlscan_LiverWeight, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Weight", ylim = c(0,11))
abline(h = threshold_LiverWeight, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")

#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverWeight, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
gmap_peaksLiverWeight <- find_peaks(scan1_output = qtlscan_LiverWeight, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(gmap_peaksLiverWeight)

#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverWeight, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
peaksLiverWeight <- find_peaks(scan1_output = qtlscan_LiverWeight, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(peaksLiverWeight)


write_xlsx(list("LiverWeight gmap (cM)" = gmap_peaksLiverWeight,
                "LiverWeight pmap (Mbp)" = peaksLiverWeight),
           "LiverWeight Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver Weight --- Chromosome 11
par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
chr = 11
coef_blup_LiverWeight_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverWeight_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverWeight, main = "LiverWeight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(15,35)
plot_coefCC(x = coef_blup_LiverWeight_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverWeight, main = "LiverWeight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 11
#could use ci_lo and ci_hi, but in this case I want a specific chr 2 position
#start = peaksLiverWeight[peaksLiverWeight$chr ==  chr,"ci_lo"]
#end = peaksLiverWeight[peaksLiverWeight$chr == chr, "ci_hi"] 

pander(peaksLiverWeight)
#based on peaksAST, peak of interest is ~40 Mbp
variants_LiverWeight_chr11 <- query_variants(chr, 16, 18)
out_snps_LiverWeight_chr11 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                      chr = chr, start = 16, end = 18, keep_all_snps = TRUE)
plot_snpasso(out_snps_LiverWeight_chr11$lod, out_snps_LiverWeight_chr11$snpinfo, main = "Liver Weight SNPs")

LiverWeight_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = 16, end = 18)
plot(out_snps_LiverWeight_chr11$lod, out_snps_LiverWeight_chr11$snpinfo, drop_hilit=1.5, genes = LiverWeight_Genes_MGI_chr11, main = "Liver Weight Genes MGI")
plot_genes(LiverWeight_Genes_MGI_chr11, main = "Liver Weight Genes MGI")


dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Liver Weight GWAS - RankZ sexgen.pdf")
out_gwas_LiverWeight <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeight"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverWeight$lod, out_gwas_LiverWeight$snpinfo, altcol="green4", gap=0, main = "LiverWeight GWAS", ylim = c(0,6))
dev.off()

out_gwas_BodyWeight <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zBodyWeight"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)

out_gwas_LiverWeightBodyWeight <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)


####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverWeight_sex <- est_herit(pheno["zLiverWeight"], kinship_lmm, sex, cores = 2)
herit_LiverWeight_sexgen <- est_herit(pheno["zLiverWeight"], kinship_lmm, sexgen, cores = 2)


herit_BodyWeight_sex <- est_herit(pheno["zBodyWeight"], kinship_lmm, sex, cores = 2)
herit_BodyWeight_sexgen <- est_herit(pheno["zBodyWeight"], kinship_lmm, sexgen, cores = 2)

herit_LiverWeightBodyWeight_sex <- est_herit(pheno["zLiverWeightBodyWeight"], kinship_lmm, sex, cores = 2)
herit_LiverWeightBodyWeight_sexgen <- est_herit(pheno["zLiverWeightBodyWeight"], kinship_lmm, sexgen, cores = 2)
