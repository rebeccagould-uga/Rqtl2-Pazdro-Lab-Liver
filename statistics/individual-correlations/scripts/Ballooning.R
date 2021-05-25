#HydropicDegeneration Statistical Analysis
#Kruskal Wallis tests + summary statistics
#created by Becca Gould
#updated December 2020

#helpful link: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r


data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data_na-rows_deleted.csv")
head(data, 6)

data$group <- ordered(data$HydropicDegeneration,
                      levels = c("0", "1", "2", "3", "4", "5"))


library("dplyr")
library("ggpubr")
library ("ggsignif")

#install.packages("ggpubr")
## Install
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

my_comparisons <- list( c("0", "1"), c("0", "2"), c("0", "3"), c("0", "4"), c("0", "5"),
                        c("1", "2"), c("1", "3"), c("1", "4"),c("1", "5"),
                        c("2", "3"), c("2", "4"),c("2", "5"),
                        c("3", "4"),c("3", "5"),
                        c("4", "5"))

#summary stats - GSH by HydropicDegeneration grade
#can do for any phenotype of interest
GSHstats <- group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(GSH, na.rm = TRUE),
    sd = sd(GSH, na.rm = TRUE),
    median = median(GSH, na.rm = TRUE),
    IQR = IQR(GSH, na.rm = TRUE)
  )


# Box plots
# ++++++++++++++++++++
# Plot variable by HydropicDegeneration grade and color by grade

#reference to get all p-values
#ggboxplot(data, x = "HydropicDegeneration", y = "GSH",
#         color = "HydropicDegeneration", palette = "jco",
#        ylab = "GSH (nmol/mg)", xlab = "Hydropic Degeneration Grade")+
#stat_compare_means(comparisons = my_comparisons, label = "p.format")

  p1 <- ggboxplot(data, x = "HydropicDegeneration", y = "GSH",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "GSH (nmol/mg)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p2 <- ggboxplot(data, x = "HydropicDegeneration", y = "GSSG",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "GSSG (nmol/mg)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p3 <- ggboxplot(data, x = "HydropicDegeneration", y = "TotalGSH",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "Total Glutathione (nmol/mg)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  
  p4 <- ggboxplot(data, x = "HydropicDegeneration", y = "GSHGSSGRatio",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "GSH/GSSG", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p5 <- ggboxplot(data, x = "HydropicDegeneration", y = "Eh",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "Eh (mV)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p6 <- ggboxplot(data, x = "HydropicDegeneration", y = "NADH",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "NADH (pmol/ug)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p7 <- ggboxplot(data, x = "HydropicDegeneration", y = "NADP",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "NADP (pmol/ug)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p8 <- ggboxplot(data, x = "HydropicDegeneration", y = "NADPH",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "NADPH (pmol/ug)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p9 <- ggboxplot(data, x = "HydropicDegeneration", y = "NADPNADPHRatio",
                  color = "HydropicDegeneration", palette = "jco",
                  ylab = "NADP/NADPH", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p10 <- ggboxplot(data, x = "HydropicDegeneration", y = "AST",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "AST (U/L)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p11 <- ggboxplot(data, x = "HydropicDegeneration", y = "ALT",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "ALT (U/L)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p12 <- ggboxplot(data, x = "HydropicDegeneration", y = "LiverWeightBodyWeight",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "Liver Weight/Body Weight (%)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p13 <- ggboxplot(data, x = "HydropicDegeneration", y = "Glucose",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "Glucose (mg/dL)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p14 <- ggboxplot(data, x = "HydropicDegeneration", y = "LiverWeight",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "Liver Weight (g)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p15 <- ggboxplot(data, x = "HydropicDegeneration", y = "Steatosis",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "Steatosis Grade", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p16 <- ggboxplot(data, x = "HydropicDegeneration", y = "ASTALTRatio",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "AST/ALT", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")
  
  p17 <- ggboxplot(data, x = "HydropicDegeneration", y = "BodyWeight",
                   color = "HydropicDegeneration", palette = "jco",
                   ylab = "Body Weight (g)", xlab = "Hydropic Degeneration Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test", label.y = 100)+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0")


  
  
    
  # # Mean plots
  # # ++++++++++++++++++++
  # # Plot GSH by HydropicDegeneration grade
  # # Add error bars: mean_se
  # # (other values include: mean_sd, mean_ci, median_iqr, ....)
  # p2 <- ggline(data, x = "HydropicDegeneration", y = "GSH", 
  #        color = "Black",
  #        #color = "HydropicDegeneration", palette = "npg",
  #        add = c("mean_se", "jitter"),
  #        ylab = "GSH (nmol/mg)", xlab = "Hydropic Degeneration Grade")
  
  pdf(file = "HydropicDegeneration-plots-pvalues-symbol.pdf")
  ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, ncol = 1, nrow = 1)
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
  KW_GSH <- compare_means(formula = GSH ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSSG <- compare_means(formula = GSSG ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_TotalGSH <- compare_means(formula = TotalGSH ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSHGSSGRatio <- compare_means(formula = GSHGSSGRatio ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_EhGSSG2GSH <- compare_means(formula = EhGSSG2GSH ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADH <- compare_means(formula = NADH ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADP <- compare_means(formula = NADP ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPH <- compare_means(formula = NADPH ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPNADPHRatio <- compare_means(formula = NADPNADPHRatio ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_AST <- compare_means(formula = AST ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_ALT <- compare_means(formula = ALT ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_BUN <- compare_means(formula = BUN ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_Glucose <- compare_means(formula = Glucose ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_LiverWeight <- compare_means(formula = LiverWeight ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  KW_Steatosis <- compare_means(formula = Steatosis ~ HydropicDegeneration, data = data, method = "kruskal.test", paired = FALSE)
  
  
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
                  "Steatosis"=KW_Steatosis),
             "KWtests-HydropicDegeneration.xlsx")
  
  
  library("ggplot2")
  library("ggsignif")
  
  
  # Visualize: Specify the comparisons you want
  my_comparisons <- list( c("0", "1"), c("0", "2"), c("0", "3"), c("0", "5"),
                          c("1", "2"), c("1", "3"), c("1", "5"),
                          c("2", "3"), c("2", "5"),
                          c("3", "5"))
  
  


# #individual KW tests for reference
# KW_GSH <- kruskal.test(GSH ~ HydropicDegeneration, data = data)
# KW_GSSG <- kruskal.test(GSSG ~ HydropicDegeneration, data = data)
# KW_TotalGSH <- kruskal.test(TotalGSH ~ HydropicDegeneration, data = data)
# KW_GSHGSSGRatio <- kruskal.test(GSHGSSGRatio ~ HydropicDegeneration, data = data)
# KW_Eh <- kruskal.test(EhGSSG2GSH ~ HydropicDegeneration, data = data)
# KW_NADH <- kruskal.test(NADH ~ HydropicDegeneration, data = data)
# KW_NADP <- kruskal.test(NADP ~ HydropicDegeneration, data = data)
# KW_NADPH <- kruskal.test(NADPH ~ HydropicDegeneration, data = data)
# KW_NADPNADPHRatio <- kruskal.test(NADPNADPHRatio ~ HydropicDegeneration, data = data)
# KW_AST <- kruskal.test(AST ~ HydropicDegeneration, data = data)
# KW_ALT <- kruskal.test(ALT ~ HydropicDegeneration, data = data)
# KW_BUN <- kruskal.test(BUN ~ HydropicDegeneration, data = data)
# KW_Glucose <- kruskal.test(Glucose ~ HydropicDegeneration, data = data)
# KW_LiverWeight <- kruskal.test(LiverWeight ~ HydropicDegeneration, data = data)
# KW_HydropicDegeneration <- kruskal.test(HydropicDegeneration ~ HydropicDegeneration, data = data)
# 
# #to get them all in one print out
# KWtests = list(KW_GSH, KW_GSSG, KW_TotalGSH, KW_GSHGSSGRatio, KW_Eh, KW_NADH, KW_NADP, KW_NADPH, KW_NADPNADPHRatio, KW_ALT, KW_AST, KW_BUN, KW_Glucose, KW_LiverWeight, KW_HydropicDegeneration)
# KWtests
# KWtests_corrected <- data.frame(unlist(KWtests))
# #setwd
# write.table(KWtests_corrected, "~/Rqtl2-Glutathione-Genetics/correlations/KWtests_HydropicDegeneration.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# 
# 
# 
# 
