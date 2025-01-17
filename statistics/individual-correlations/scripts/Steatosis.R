#Steatosis Statistical Analysis
#Kruskal Wallis tests + summary statistics
#created by Becca Gould
#updated December 2020

#helpful link: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r


data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data-Liver-Histology-QTL-paper-narowsdeleted.csv", na.strings = "NA")
head(data, 6)

data$group <- ordered(data$Steatosis,
                         levels = c("0", "1", "2", "3", "4", "5"))


library("dplyr")
library("ggpubr")
library ("ggsignif")

# Visualize: Specify the comparisons you want
my_comparisons <- list( c("0", "1"), c("0", "2"), c("0", "3"), c("0", "5"),
                        c("1", "2"), c("1", "3"), c("1", "5"),
                        c("2", "3"), c("2", "5"),
                        c("3", "5"))


#install.packages("ggpubr")
## Install
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")


#summary stats - GSH by Steatosis grade
#can do for any phenotype of interest
GSHstats <- group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(GSH, na.rm = TRUE),
    sd = sd(GSH, na.rm = TRUE),
    median = median(GSH, na.rm = TRUE),
    IQR = IQR(GSH, na.rm = TRUE)
  )

SteatosisStats <- group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(Steatosis, na.rm = TRUE),
    sd = sd(Steatosis, na.rm = TRUE),
    median = median(Steatosis, na.rm = TRUE),
    IQR = IQR(Steatosis, na.rm = TRUE)
  )

#by sex
SteatosisStatsSex <- group_by(data, group, Sex) %>%
  summarise(
    count = n(),
    mean = mean(Steatosis, na.rm = TRUE),
    sd = sd(Steatosis, na.rm = TRUE),
    median = median(Steatosis, na.rm = TRUE),
    IQR = IQR(Steatosis, na.rm = TRUE)
  )
    
# Box plots
# ++++++++++++++++++++
# Plot variable by Steatosis grade and color by grade

#reference to get all p-values
#ggboxplot(data, x = "Steatosis", y = "GSH",
 #         color = "Steatosis", palette = "jco",
  #        ylab = "GSH (nmol/mg)", xlab = "Steatosis Grade")+
  #stat_compare_means(comparisons = my_comparisons, label = "p.format")

  p1 <- ggboxplot(data, x = "Steatosis", y = "LiverGSH",
            color = "Steatosis", palette = "jco",
            ylab = "Liver GSH (nmol/mg)", xlab = "Steatosis Grade")+ 
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  # ggboxplot(data, x = "Steatosis", y = "GSH", color = "Steatosis", 
  #           add = "jitter", legend = "none") +
  #   rotate_x_text(angle = 45)+
  #   #geom_hline(yintercept = mean(data$GSH), linetype = 2)+ # Add horizontal line at base mean
  #   stat_compare_means(method = "kruskal.test")+        # Add global annova p-value
  #   stat_compare_means(label = "p.signif", method = "t.test",
  #                      ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all
  
  p2 <- ggboxplot(data, x = "Steatosis", y = "LiverGSSG",
                  color = "Steatosis", palette = "jco",
                  ylab = "Liver GSSG (nmol/mg)", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p3 <- ggboxplot(data, x = "Steatosis", y = "LiverTotalGSH",
                  color = "Steatosis", palette = "jco",
                  ylab = "Liver Total Glutathione (nmol/mg)", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p4 <- ggboxplot(data, x = "Steatosis", y = "LiverGSHGSSGRatio",
                  color = "Steatosis", palette = "jco",
                  ylab = "Liver GSH/GSSG", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p5 <- ggboxplot(data, x = "Steatosis", y = "LiverEh",
                  color = "Steatosis", palette = "jco",
                  ylab = "Liver Eh (mV)", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p6 <- ggboxplot(data, x = "Steatosis", y = "NADH",
                  color = "Steatosis", palette = "jco",
                  ylab = "NADH (pmol/ug)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p7 <- ggboxplot(data, x = "Steatosis", y = "NADP",
                  color = "Steatosis", palette = "jco",
                  ylab = "NADP (pmol/ug)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p8 <- ggboxplot(data, x = "Steatosis", y = "NADPH",
                  color = "Steatosis", palette = "jco",
                  ylab = "NADPH (pmol/ug)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p9 <- ggboxplot(data, x = "Steatosis", y = "NADPNADPHRatio",
                  color = "Steatosis", palette = "jco",
                  ylab = "NADP/NADPH", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p10 <- ggboxplot(data, x = "Steatosis", y = "AST",
                  color = "Steatosis", palette = "jco",
                  ylab = "AST (U/L)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p11 <- ggboxplot(data, x = "Steatosis", y = "ALT",
                   color = "Steatosis", palette = "jco",
                   ylab = "ALT (U/L)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
    
  p12 <- ggboxplot(data, x = "Steatosis", y = "LiverWeightBodyWeight",
                   color = "Steatosis", palette = "jco",
                   ylab = "Liver Weight/Body Weight (%)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p13 <- ggboxplot(data, x = "Steatosis", y = "Glucose",
                   color = "Steatosis", palette = "jco",
                   ylab = "Glucose (mg/dL)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p14 <- ggboxplot(data, x = "Steatosis", y = "LiverWeight",
                   color = "Steatosis", palette = "jco",
                   ylab = "Liver Weight (g)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
    #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p15 <- ggboxplot(data, x = "Steatosis", y = "HydropicDegeneration",
                   color = "Steatosis", palette = "jco",
                   ylab = "Hydropic Degeneration", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p16 <- ggboxplot(data, x = "Steatosis", y = "ASTALTRatio",
                   color = "Steatosis", palette = "jco",
                   ylab = "AST/ALT", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p17 <- ggboxplot(data, x = "Steatosis", y = "BodyWeight",
                   color = "Steatosis", palette = "jco",
                   ylab = "Body Weight (g)", xlab = "Steatosis Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p18 <- ggboxplot(data, x = "Steatosis", y = "KidneyGSH",
                  color = "Steatosis", palette = "jco",
                  ylab = "Kidney GSH (nmol/mg)", xlab = "Steatosis Grade")+ 
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  

  p19 <- ggboxplot(data, x = "Steatosis", y = "KidneyGSSG",
                  color = "Steatosis", palette = "jco",
                  ylab = "Kidney GSSG (nmol/mg)", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p20 <- ggboxplot(data, x = "Steatosis", y = "KidneyTotalGSH",
                  color = "Steatosis", palette = "jco",
                  ylab = "Kidney Total Glutathione (nmol/mg)", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p21 <- ggboxplot(data, x = "Steatosis", y = "KidneyGSHGSSGRatio",
                  color = "Steatosis", palette = "jco",
                  ylab = "Kidney GSH/GSSG", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p22 <- ggboxplot(data, x = "Steatosis", y = "KidneyEh",
                  color = "Steatosis", palette = "jco",
                  ylab = "Kidney Eh (mV)", xlab = "Steatosis Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  
  
  # # Mean plots
  # # ++++++++++++++++++++
  # # Plot GSH by Steatosis grade
  # # Add error bars: mean_se
  # # (other values include: mean_sd, mean_ci, median_iqr, ....)
  # p2 <- ggline(data, x = "Steatosis", y = "GSH", 
  #        color = "Black",
  #        #color = "Steatosis", palette = "npg",
  #        add = c("mean_se", "jitter"),
  #        ylab = "GSH (nmol/mg)", xlab = "Steatosis Grade")
  
  pdf(file = "Steatosis-plots-pvalues-symbol.pdf")
  ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, ncol = 1, nrow = 1)
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
  KW_GSH <- compare_means(formula = GSH ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSSG <- compare_means(formula = GSSG ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_TotalGSH <- compare_means(formula = TotalGSH ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSHGSSGRatio <- compare_means(formula = GSHGSSGRatio ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_EhGSSG2GSH <- compare_means(formula = EhGSSG2GSH ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADH <- compare_means(formula = NADH ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADP <- compare_means(formula = NADP ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPH <- compare_means(formula = NADPH ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPNADPHRatio <- compare_means(formula = NADPNADPHRatio ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_AST <- compare_means(formula = AST ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_ALT <- compare_means(formula = ALT ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_BUN <- compare_means(formula = BUN ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_Glucose <- compare_means(formula = Glucose ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_LiverWeight <- compare_means(formula = LiverWeight ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  KW_Ballooning <- compare_means(formula = Ballooning ~ Steatosis, data = data, method = "kruskal.test", paired = FALSE)
  
  
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
             "KWtests-steatosis.xlsx")
  
  
  library("ggplot2")
  library("ggsignif")
  

# #individual KW tests for reference
# KW_GSH <- kruskal.test(GSH ~ Steatosis, data = data)
# KW_GSSG <- kruskal.test(GSSG ~ Steatosis, data = data)
# KW_TotalGSH <- kruskal.test(TotalGSH ~ Steatosis, data = data)
# KW_GSHGSSGRatio <- kruskal.test(GSHGSSGRatio ~ Steatosis, data = data)
# KW_Eh <- kruskal.test(EhGSSG2GSH ~ Steatosis, data = data)
# KW_NADH <- kruskal.test(NADH ~ Steatosis, data = data)
# KW_NADP <- kruskal.test(NADP ~ Steatosis, data = data)
# KW_NADPH <- kruskal.test(NADPH ~ Steatosis, data = data)
# KW_NADPNADPHRatio <- kruskal.test(NADPNADPHRatio ~ Steatosis, data = data)
# KW_AST <- kruskal.test(AST ~ Steatosis, data = data)
# KW_ALT <- kruskal.test(ALT ~ Steatosis, data = data)
# KW_BUN <- kruskal.test(BUN ~ Steatosis, data = data)
# KW_Glucose <- kruskal.test(Glucose ~ Steatosis, data = data)
# KW_LiverWeight <- kruskal.test(LiverWeight ~ Steatosis, data = data)
# KW_Ballooning <- kruskal.test(Ballooning ~ Steatosis, data = data)
# 
# #to get them all in one print out
# KWtests = list(KW_GSH, KW_GSSG, KW_TotalGSH, KW_GSHGSSGRatio, KW_Eh, KW_NADH, KW_NADP, KW_NADPH, KW_NADPNADPHRatio, KW_ALT, KW_AST, KW_BUN, KW_Glucose, KW_LiverWeight, KW_Ballooning)
# KWtests
# KWtests_corrected <- data.frame(unlist(KWtests))
# #setwd
# write.table(KWtests_corrected, "~/Rqtl2-Glutathione-Genetics/correlations/KWtests_Steatosis.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# 
# 
# 
# 
