#Fibrosis Statistical Analysis
#Kruskal Wallis tests + summary statistics
#created by Becca Gould
#updated December 2020

#helpful link: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r


data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data-Liver-Histology-QTL-paper-narowsdeleted.csv", na.strings = "NA")
head(data, 6)

data$group <- ordered(data$Fibrosis,
                         levels = c("0", "1", "2"))


library("dplyr")
library("ggpubr")

# Visualize: Specify the comparisons you want
my_comparisons <- list( c("0", "1"), c("0", "2"),
                        c("1", "2"))


#install.packages("ggpubr")
## Install
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")


#summary stats - GSH by Fibrosis grade
#can do for any phenotype of interest
GSHstats <- group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(GSH, na.rm = TRUE),
    sd = sd(GSH, na.rm = TRUE),
    median = median(GSH, na.rm = TRUE),
    IQR = IQR(GSH, na.rm = TRUE)
  )

FibrosisStats <- group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(Fibrosis, na.rm = TRUE),
    sd = sd(Fibrosis, na.rm = TRUE),
    median = median(Fibrosis, na.rm = TRUE),
    IQR = IQR(Fibrosis, na.rm = TRUE)
  )

#by sex
FibrosisStatsSex <- group_by(data, group, Sex) %>%
  summarise(
    count = n(),
    mean = mean(Fibrosis, na.rm = TRUE),
    sd = sd(Fibrosis, na.rm = TRUE),
    median = median(Fibrosis, na.rm = TRUE),
    IQR = IQR(Fibrosis, na.rm = TRUE)
  )
# Box plots
# ++++++++++++++++++++
# Plot variable by Fibrosis grade and color by grade

  p1 <- ggboxplot(data, x = "Fibrosis", y = "LiverGSH",
            color = "Fibrosis", palette = "jco",
            ylab = "GSH (nmol/mg)", xlab = "Fibrosis Grade")#+ 
  p1 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  # ggboxplot(data, x = "Fibrosis", y = "GSH", color = "Fibrosis", 
  #           add = "jitter", legend = "none") +
  #   rotate_x_text(angle = 45)+
  #   #geom_hline(yintercept = mean(data$GSH), linetype = 2)+ # Add horizontal line at base mean
  #   stat_compare_means(method = "kruskal.test")+        # Add global annova p-value
  #   stat_compare_means(label = "p.signif", method = "t.test",
  #                      ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all
  
  p2 <- ggboxplot(data, x = "Fibrosis", y = "LiverGSSG",
                  color = "Fibrosis", palette = "jco",
                  ylab = "GSSG (nmol/mg)", xlab = "Fibrosis Grade")#+ 
  p2 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p3 <- ggboxplot(data, x = "Fibrosis", y = "LiverTotalGSH",
                  color = "Fibrosis", palette = "jco",
                  ylab = "Total Glutathione (nmol/mg)", xlab = "Fibrosis Grade")#+ 
  p3 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p4 <- ggboxplot(data, x = "Fibrosis", y = "LiverGSHGSSGRatio",
                  color = "Fibrosis", palette = "jco",
                  ylab = "GSH/GSSG", xlab = "Fibrosis Grade")#+ 
  p4 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p5 <- ggboxplot(data, x = "Fibrosis", y = "LiverEh",
                  color = "Fibrosis", palette = "jco",
                  ylab = "Eh", xlab = "Fibrosis Grade")#+ 
  p5 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p6 <- ggboxplot(data, x = "Fibrosis", y = "NADH",
                  color = "Fibrosis", palette = "jco",
                  ylab = "NADH (pmol/ug)", xlab = "Fibrosis Grade")#+ 
  p6 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p7 <- ggboxplot(data, x = "Fibrosis", y = "NADP",
                  color = "Fibrosis", palette = "jco",
                  ylab = "NADP (pmol/ug)", xlab = "Fibrosis Grade")#+ 
  p7 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p8 <- ggboxplot(data, x = "Fibrosis", y = "NADPH",
                  color = "Fibrosis", palette = "jco",
                  ylab = "NADPH (pmol/ug)", xlab = "Fibrosis Grade")#+ 
  p8 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p9 <- ggboxplot(data, x = "Fibrosis", y = "NADPNADPHRatio",
                  color = "Fibrosis", palette = "jco",
                  ylab = "NADP/NADPH", xlab = "Fibrosis Grade")#+ 
  p9 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p10 <- ggboxplot(data, x = "Fibrosis", y = "AST",
                  color = "Fibrosis", palette = "jco",
                  ylab = "AST (U/L)", xlab = "Fibrosis Grade")#+ 
  p10 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p11 <- ggboxplot(data, x = "Fibrosis", y = "ALT",
                   color = "Fibrosis", palette = "jco",
                   ylab = "ALT (U/L)", xlab = "Fibrosis Grade")#+ 
  p11 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
    
  p12 <- ggboxplot(data, x = "Fibrosis", y = "BUN",
                   color = "Fibrosis", palette = "jco",
                   ylab = "BUN (mg/dL)", xlab = "Fibrosis Grade")#+ 
  p12 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p13 <- ggboxplot(data, x = "Fibrosis", y = "Glucose",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Glucose (mg/dL)", xlab = "Fibrosis Grade")#+ 
  p13 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p14 <- ggboxplot(data, x = "Fibrosis", y = "LiverWeight",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Liver Weight (g)", xlab = "Fibrosis Grade")#+ 
  p14 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
    #stat_compare_means(method = "kruskal.test")+ # Add global p-value
    #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p15 <- ggboxplot(data, x = "Fibrosis", y = "Ballooning",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Hydropic Degeneration", xlab = "Fibrosis Grade")#+ 
  p15 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p16 <- ggboxplot(data, x = "Fibrosis", y = "ASTALTRatio",
                   color = "Fibrosis", palette = "jco",
                   ylab = "AST/ALT", xlab = "Fibrosis Grade")#+ 
  p16 + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p17 <- ggboxplot(data, x = "Fibrosis", y = "KidneyGSH",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Kidney GSH (nmol/mg)", xlab = "Fibrosis Grade")+ 
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format") #, symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p18 <- ggboxplot(data, x = "Fibrosis", y = "KidneyGSSG",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Kidney GSSG (nmol/mg)", xlab = "Fibrosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format") #, symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p19 <- ggboxplot(data, x = "Fibrosis", y = "KidneyTotalGSH",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Kidney Total Glutathione (nmol/mg)", xlab = "Fibrosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format") #, symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p20 <- ggboxplot(data, x = "Fibrosis", y = "KidneyGSHGSSGRatio",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Kidney GSH/GSSG", xlab = "Fibrosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format") #, symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p21 <- ggboxplot(data, x = "Fibrosis", y = "KidneyEh",
                   color = "Fibrosis", palette = "jco",
                   ylab = "Kidney Eh (mV)", xlab = "Fibrosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format") #, symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  
  
  # # Mean plots
  # # ++++++++++++++++++++
  # # Plot GSH by Fibrosis grade
  # # Add error bars: mean_se
  # # (other values include: mean_sd, mean_ci, median_iqr, ....)
  # p2 <- ggline(data, x = "Fibrosis", y = "GSH", 
  #        color = "Black",
  #        #color = "Fibrosis", palette = "npg",
  #        add = c("mean_se", "jitter"),
  #        ylab = "GSH (nmol/mg)", xlab = "Fibrosis Grade")
  
  pdf(file = "Fibrosis-plots-pvalues.pdf")
  ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, ncol = 1, nrow = 1)
  dev.off()
  
  
  
  
  
############################################
## Spearman's Correlations
############################################

  #use Spearman's correlations, not Kruskal-Wallis
  #see "CorrelationMatrix.R"
  

############################################
## Met with Dr. Kim Love (statistician) and she said Spearman's is preferred compared to Kruskal Wallis
## included Kruskal-Wallis tests below for records, but use Spearman's
############################################

  #Kruskal Wallis Tests
  KW_GSH <- compare_means(formula = GSH ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSSG <- compare_means(formula = GSSG ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_TotalGSH <- compare_means(formula = TotalGSH ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSHGSSGRatio <- compare_means(formula = GSHGSSGRatio ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_EhGSSG2GSH <- compare_means(formula = EhGSSG2GSH ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADH <- compare_means(formula = NADH ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADP <- compare_means(formula = NADP ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPH <- compare_means(formula = NADPH ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPNADPHRatio <- compare_means(formula = NADPNADPHRatio ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_AST <- compare_means(formula = AST ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_ALT <- compare_means(formula = ALT ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_BUN <- compare_means(formula = BUN ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_Glucose <- compare_means(formula = Glucose ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_LiverWeight <- compare_means(formula = LiverWeight ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_Ballooning <- compare_means(formula = Ballooning ~ Fibrosis, data = data, method = "kruskal.test", paired = FALSE)
  
  
  library(writexl)
  write_xlsx(list("GSH" = KW_GSH,
                  "GSSG" = KW_GSSG,
                  "TotalGSH" = KW_TotalGSH,
                  "GSHGSSGRatio" = KW_GSHGSSGRatio,
                  "EhGSSG2GSH" = KW_EhGSSG2GSH,
                  "NADH" = KW_NADH,
                  "NADP"=KW_NADP,
                  "NADPH"=KW_NADPH,
                  "NADPNADPHRatio"=KW_NADPNADPHRatio,
                  "AST"=KW_AST,
                  "ALT"=KW_ALT,
                  "BUN"=KW_BUN,
                  "Glucose"=KW_Glucose,
                  "LiverWeight"=KW_LiverWeight,
                  "Ballooning"=KW_Ballooning),
             "KWtests-Fibrosis.xlsx")
  
  
  library("ggplot2")
  library("ggsignif")
  

# #individual KW tests for reference
# KW_GSH <- kruskal.test(GSH ~ Fibrosis, data = data)
# KW_GSSG <- kruskal.test(GSSG ~ Fibrosis, data = data)
# KW_TotalGSH <- kruskal.test(TotalGSH ~ Fibrosis, data = data)
# KW_GSHGSSGRatio <- kruskal.test(GSHGSSGRatio ~ Fibrosis, data = data)
# KW_Eh <- kruskal.test(EhGSSG2GSH ~ Fibrosis, data = data)
# KW_NADH <- kruskal.test(NADH ~ Fibrosis, data = data)
# KW_NADP <- kruskal.test(NADP ~ Fibrosis, data = data)
# KW_NADPH <- kruskal.test(NADPH ~ Fibrosis, data = data)
# KW_NADPNADPHRatio <- kruskal.test(NADPNADPHRatio ~ Fibrosis, data = data)
# KW_AST <- kruskal.test(AST ~ Fibrosis, data = data)
# KW_ALT <- kruskal.test(ALT ~ Fibrosis, data = data)
# KW_BUN <- kruskal.test(BUN ~ Fibrosis, data = data)
# KW_Glucose <- kruskal.test(Glucose ~ Fibrosis, data = data)
# KW_LiverWeight <- kruskal.test(LiverWeight ~ Fibrosis, data = data)
# KW_Ballooning <- kruskal.test(Ballooning ~ Fibrosis, data = data)
# 
# #to get them all in one print out
# KWtests = list(KW_GSH, KW_GSSG, KW_TotalGSH, KW_GSHGSSGRatio, KW_Eh, KW_NADH, KW_NADP, KW_NADPH, KW_NADPNADPHRatio, KW_ALT, KW_AST, KW_BUN, KW_Glucose, KW_LiverWeight, KW_Ballooning)
# KWtests
# KWtests_corrected <- data.frame(unlist(KWtests))
# #setwd
# write.table(KWtests_corrected, "~/Rqtl2-Glutathione-Genetics/correlations/KWtests_Fibrosis.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# 
# 
# 
# 
