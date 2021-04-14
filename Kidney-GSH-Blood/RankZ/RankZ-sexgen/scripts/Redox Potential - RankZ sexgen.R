# R01 GSH DO Mapping Code 
# Updated December 2020
# Becca Gould 

#KIDNEY GLUTATHIONE + BLOOD (BUN) MAPPING - Redox Potential GSSG/2GSH

#Load in Kidney-QTL-Mapping-RankZ-sexgen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


#setwd

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

#for Kidney Redox Potential GSSG 2GSH
qtlscan_KidneyRedoxPotentialGSSG2GSH<- scan1(genoprobs = probs, pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_KidneyRedoxPotentialGSSG2GSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Redox Potential GSSG 2GSH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyRedoxPotentialGSSG2GSH = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney Redox Potential GSSG/2GSH", ylim = c(0,11))
  abline(h = threshold_KidneyRedoxPotentialGSSG2GSH, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney Redox Potential GSSG/2GSH", ylim = c(0,11))
  abline(h = threshold_KidneyRedoxPotentialGSSG2GSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksKidneyRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksKidneyRedoxPotentialGSSG2GSH)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksKidneyRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksKidneyRedoxPotentialGSSG2GSH)
  
  
  write_xlsx(list("KidneyRedoxPotential gmap (cM)" = gmap_peaksKidneyRedoxPotentialGSSG2GSH,
                  "KidneyRedoxPotential pmap (Mbp)" = peaksKidneyRedoxPotentialGSSG2GSH),
             "Redox Potential GSSG 2GSH Peaks - RankZ sexgen.xlsx")


  #to compare with Liver GSH
  #saveRDS(qtlscan_KidneyRedoxPotentialGSSG2GSH, file = "~/Rqtl2-Glutathione-Genetics/Kidney-GSH-Blood/QTL-Eh-Kidney.rds")
  
  
####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Kidney Eh --- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "14"
  coef_blup_KidneyRedoxPotential_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyRedoxPotential_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Eh BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,25)
  plot_coefCC(x = coef_blup_KidneyRedoxPotential_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Eh BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "14"
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 14 peak
  #start = peaksKidneyRedoxPotentialGSSG2GSH[peaksKidneyRedoxPotentialGSSG2GSH$chr ==  chr,"ci_lo"]
  #end = peaksKidneyRedoxPotentialGSSG2GSH[peaksKidneyRedoxPotentialGSSG2GSH$chr == chr, "ci_hi"] 
  
  pander(peaksKidneyRedoxPotentialGSSG2GSH)
  #based on peaksGSH, peak of interest is ~100 Mbp
  variants_KidneyRedoxPotential_chr14 <- query_variants(chr, 21, 25)
  out_snps_KidneyRedoxPotential_chr14 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                        chr = chr, start = 21, end = 25, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyRedoxPotential_chr14$lod, out_snps_KidneyRedoxPotential_chr14$snpinfo, main = "Kidney Eh SNPs")
  
  KidneyEh_Genes_MGI_chr14 <- query_genes_mgi(chr = chr, start = 21, end = 25)
  plot(out_snps_KidneyRedoxPotential_chr14$lod, out_snps_KidneyRedoxPotential_chr14$snpinfo, drop_hilit=1.5, genes = KidneyEh_Genes_MGI_chr14, main = "Kidney Eh Genes MGI")
  
  plot_genes(KidneyEh_Genes_MGI_chr14, main = "Kidney Eh Genes MGI")
  
dev.off()

  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Redox-Potential-GSSG-2GSH-GWAS-RankZ-sexgen.pdf")
out_gwas_KidneyRedoxPotentialGSSG2GSH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_KidneyRedoxPotentialGSSG2GSH$lod, out_gwas_KidneyRedoxPotentialGSSG2GSH$snpinfo, altcol="green4", gap=0, main = "Kidney Redox Potential GSSG/2GSH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_KidneyRedoxPotentialGSSG2GSH_sex <- est_herit(pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship_lmm, sex, cores = 2)
herit_KidneyRedoxPotentialGSSG2GSH_sexgen <- est_herit(pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship_lmm, sexgen, cores = 2)



##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "Redox Potential Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gpx1 Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gclc Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_KidneyRedoxPotentialGSSG2GSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gclm Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_KidneyRedoxPotentialGSSG2GSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gss Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_KidneyRedoxPotentialGSSG2GSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gsr Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################



