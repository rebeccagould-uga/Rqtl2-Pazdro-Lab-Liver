#Inflammation Statistical Analysis
#Kruskal Wallis tests + summary statistics
#created by Becca Gould
#updated May 2021

#helpful link: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r


data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data-Liver-Histology-QTL-paper-narowsdeleted.csv", na.strings = "NA")
head(data, 6)

data$group <- ordered(data$Inflammation,
                         levels = c("0", "1", "2", "3", "4", "5"))


library("dplyr")
library("ggpubr")
library ("ggsignif")

# Visualize: Specify the comparisons you want
my_comparisons <- list( c("0", "1"), c("0", "2"), c("0", "3"), c("0", "4"),
                        c("1", "2"), c("1", "3"), c("1", "4"),
                        c("2", "3"), c("2", "4"),
                        c("3", "4") )

#install.packages("ggpubr")
## Install
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")


#summary stats - GSH by Inflammation grade
#can do for any phenotype of interest
GSHstats <- group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(GSH, na.rm = TRUE),
    sd = sd(GSH, na.rm = TRUE),
    median = median(GSH, na.rm = TRUE),
    IQR = IQR(GSH, na.rm = TRUE)
  )

InflammationStats <- group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(Inflammation, na.rm = TRUE),
    sd = sd(Inflammation, na.rm = TRUE),
    median = median(Inflammation, na.rm = TRUE),
    IQR = IQR(Inflammation, na.rm = TRUE)
  )

#by sex
InflammationStatsSex <- group_by(data, group, Sex) %>%
  summarise(
    count = n(),
    mean = mean(Inflammation, na.rm = TRUE),
    sd = sd(Inflammation, na.rm = TRUE),
    median = median(Inflammation, na.rm = TRUE),
    IQR = IQR(Inflammation, na.rm = TRUE)
  )


# Box plots
# ++++++++++++++++++++
# Plot variable by Inflammation grade and color by grade

#reference to get all p-values
#ggboxplot(data, x = "Inflammation", y = "GSH",
 #         color = "Inflammation", palette = "jco",
  #        ylab = "GSH (nmol/mg)", xlab = "Inflammation Grade")+
  #stat_compare_means(comparisons = my_comparisons, label = "p.format")

  p1 <- ggboxplot(data, x = "Inflammation", y = "GSH",
            color = "Inflammation", palette = "jco",
            ylab = "GSH (nmol/mg)", xlab = "Inflammation Grade")+ 
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  # ggboxplot(data, x = "Inflammation", y = "GSH", color = "Inflammation", 
  #           add = "jitter", legend = "none") +
  #   rotate_x_text(angle = 45)+
  #   #geom_hline(yintercept = mean(data$GSH), linetype = 2)+ # Add horizontal line at base mean
  #   stat_compare_means(method = "kruskal.test")+        # Add global annova p-value
  #   stat_compare_means(label = "p.signif", method = "t.test",
  #                      ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all
  
  p2 <- ggboxplot(data, x = "Inflammation", y = "GSSG",
                  color = "Inflammation", palette = "jco",
                  ylab = "GSSG (nmol/mg)", xlab = "Inflammation Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p3 <- ggboxplot(data, x = "Inflammation", y = "TotalGSH",
                  color = "Inflammation", palette = "jco",
                  ylab = "Total Glutathione (nmol/mg)", xlab = "Inflammation Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  
  p4 <- ggboxplot(data, x = "Inflammation", y = "GSHGSSGRatio",
                  color = "Inflammation", palette = "jco",
                  ylab = "GSH/GSSG", xlab = "Inflammation Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p5 <- ggboxplot(data, x = "Inflammation", y = "Eh",
                  color = "Inflammation", palette = "jco",
                  ylab = "Eh (mV)", xlab = "Inflammation Grade")+
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p6 <- ggboxplot(data, x = "Inflammation", y = "NADH",
                  color = "Inflammation", palette = "jco",
                  ylab = "NADH (pmol/ug)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p7 <- ggboxplot(data, x = "Inflammation", y = "NADP",
                  color = "Inflammation", palette = "jco",
                  ylab = "NADP (pmol/ug)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p8 <- ggboxplot(data, x = "Inflammation", y = "NADPH",
                  color = "Inflammation", palette = "jco",
                  ylab = "NADPH (pmol/ug)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p9 <- ggboxplot(data, x = "Inflammation", y = "NADPNADPHRatio",
                  color = "Inflammation", palette = "jco",
                  ylab = "NADP/NADPH", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p10 <- ggboxplot(data, x = "Inflammation", y = "AST",
                  color = "Inflammation", palette = "jco",
                  ylab = "AST (U/L)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p11 <- ggboxplot(data, x = "Inflammation", y = "ALT",
                   color = "Inflammation", palette = "jco",
                   ylab = "ALT (U/L)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
    
  p12 <- ggboxplot(data, x = "Inflammation", y = "LiverWeightBodyWeight",
                   color = "Inflammation", palette = "jco",
                   ylab = "Liver Weight/Body Weight (%)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p13 <- ggboxplot(data, x = "Inflammation", y = "Glucose",
                   color = "Inflammation", palette = "jco",
                   ylab = "Glucose (mg/dL)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p14 <- ggboxplot(data, x = "Inflammation", y = "LiverWeight",
                   color = "Inflammation", palette = "jco",
                   ylab = "Liver Weight (g)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
    #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p15 <- ggboxplot(data, x = "Inflammation", y = "HydropicDegeneration",
                   color = "Inflammation", palette = "jco",
                   ylab = "Hydropic Degeneration Grade", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p16 <- ggboxplot(data, x = "Inflammation", y = "ASTALTRatio",
                   color = "Inflammation", palette = "jco",
                   ylab = "AST/ALT", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p17 <- ggboxplot(data, x = "Inflammation", y = "BodyWeight",
                   color = "Inflammation", palette = "jco",
                   ylab = "Body Weight (g)", xlab = "Inflammation Grade")+ 
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  p18 <- ggboxplot(data, x = "Inflammation", y = "Steatosis",
                   color = "Inflammation", palette = "jco",
                   ylab = "Steatosis Grade", xlab = "Inflammation Grade")+ 
    stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.format", symnum.args = list(cutpoints = c(0, 0.001, 0.05, 1), symbols = c("**", "*", "ns")))
  #stat_compare_means(method = "kruskal.test")+ # Add global p-value
  #stat_compare_means(label = "p.signif", ref.group = "0") 
  
  

  
  # # Mean plots
  # # ++++++++++++++++++++
  # # Plot GSH by Inflammation grade
  # # Add error bars: mean_se
  # # (other values include: mean_sd, mean_ci, median_iqr, ....)
  # p2 <- ggline(data, x = "Inflammation", y = "GSH", 
  #        color = "Black",
  #        #color = "Inflammation", palette = "npg",
  #        add = c("mean_se", "jitter"),
  #        ylab = "GSH (nmol/mg)", xlab = "Inflammation Grade")
  
  pdf(file = "Inflammation-plots-pvalues-symbol.pdf")
  ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, ncol = 1, nrow = 1)
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
  KW_GSH <- compare_means(formula = GSH ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSSG <- compare_means(formula = GSSG ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_TotalGSH <- compare_means(formula = TotalGSH ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_GSHGSSGRatio <- compare_means(formula = GSHGSSGRatio ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_Eh <- compare_means(formula = Eh ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADH <- compare_means(formula = NADH ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADP <- compare_means(formula = NADP ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPH <- compare_means(formula = NADPH ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_NADPNADPHRatio <- compare_means(formula = NADPNADPHRatio ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_AST <- compare_means(formula = AST ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_ALT <- compare_means(formula = ALT ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_BUN <- compare_means(formula = BUN ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_Glucose <- compare_means(formula = Glucose ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_LiverWeight <- compare_means(formula = LiverWeight ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  KW_HydropicDegeneration <- compare_means(formula = HydropicDegeneration ~ Inflammation, data = data, method = "kruskal.test", paired = FALSE)
  
  
  library(writexl)
  write_xlsx(list("GSH" = KW_GSH,
                  "GSSG" = KW_GSSG,
                  "TotalGSH" = KW_TotalGSH,
                  "GSHGSSGRatio" = KW_GSHGSSGRatio,
                  "Eh" = KW_Eh,
                  "NADH" = KW_NADH,
                  "NADP"=KW_NADP,
                  "NADPH"=KW_NADPH,
                  "NADPNADPHRatio"=KW_NADPNADPHRatio,
                  "AST"=KW_AST,
                  "ALT"=KW_ALT,
                  "BUN"=KW_BUN,
                  "Glucose"=KW_Glucose,
                  "LiverWeight"=KW_LiverWeight,
                  "HydropicDegeneration"=KW_HydropicDegeneration),
             "KWtests-Inflammation.xlsx")
  
  
  library("ggplot2")
  library("ggsignif")
  

# #individual KW tests for reference
# KW_GSH <- kruskal.test(GSH ~ Inflammation, data = data)
# KW_GSSG <- kruskal.test(GSSG ~ Inflammation, data = data)
# KW_TotalGSH <- kruskal.test(TotalGSH ~ Inflammation, data = data)
# KW_GSHGSSGRatio <- kruskal.test(GSHGSSGRatio ~ Inflammation, data = data)
# KW_Eh <- kruskal.test(Eh ~ Inflammation, data = data)
# KW_NADH <- kruskal.test(NADH ~ Inflammation, data = data)
# KW_NADP <- kruskal.test(NADP ~ Inflammation, data = data)
# KW_NADPH <- kruskal.test(NADPH ~ Inflammation, data = data)
# KW_NADPNADPHRatio <- kruskal.test(NADPNADPHRatio ~ Inflammation, data = data)
# KW_AST <- kruskal.test(AST ~ Inflammation, data = data)
# KW_ALT <- kruskal.test(ALT ~ Inflammation, data = data)
# KW_BUN <- kruskal.test(BUN ~ Inflammation, data = data)
# KW_Glucose <- kruskal.test(Glucose ~ Inflammation, data = data)
# KW_LiverWeight <- kruskal.test(LiverWeight ~ Inflammation, data = data)
# KW_HydropicDegeneration <- kruskal.test(HydropicDegeneration ~ Inflammation, data = data)
# 
# #to get them all in one print out
# KWtests = list(KW_GSH, KW_GSSG, KW_TotalGSH, KW_GSHGSSGRatio, KW_Eh, KW_NADH, KW_NADP, KW_NADPH, KW_NADPNADPHRatio, KW_ALT, KW_AST, KW_BUN, KW_Glucose, KW_LiverWeight, KW_HydropicDegeneration)
# KWtests
# KWtests_corrected <- data.frame(unlist(KWtests))
# #setwd
# write.table(KWtests_corrected, "~/Rqtl2-Glutathione-Genetics/correlations/KWtests_Inflammation.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# 
# 
# 
# 
