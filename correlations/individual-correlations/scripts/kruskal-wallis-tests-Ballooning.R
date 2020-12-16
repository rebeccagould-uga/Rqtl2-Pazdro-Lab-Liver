#Ballooning Statistical Analysis
#Kruskal Wallis tests + summary statistics
#created by Becca Gould
#updated December 2020

#helpful link: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r


data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/correlations/data.csv")
head(data, 6)

data$group <- ordered(data$Ballooning,
                         levels = c("0", "1", "2", "3", "4", "5"))


#install.packages("dplyr")
library(dplyr)

#summary stats - GSH by Ballooning grade
#can do for any phenotype of interest
group_by(data, group) %>%
  summarise(
    count = n(),
    mean = mean(GSH, na.rm = TRUE),
    sd = sd(GSH, na.rm = TRUE),
    median = median(GSH, na.rm = TRUE),
    IQR = IQR(GSH, na.rm = TRUE)
  )


library("ggpubr")
#install.packages("ggpubr")
## Install
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

# Box plots
# ++++++++++++++++++++
# Plot GSH by Ballooning grade and color by grade
p1 <- ggboxplot(data, x = "Ballooning", y = "GSH", 
          color = "Ballooning", palette = "npg",
          order = c("0", "1", "2", "3", "4", "5"),
          ylab = "GSH (nmol/mg)", xlab = "Ballooning Grade")


# Mean plots
# ++++++++++++++++++++
# Plot GSH by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p2 <- ggline(data, x = "Ballooning", y = "GSH", 
       color = "Black",
       #color = "Ballooning", palette = "npg",
       add = c("mean_se", "jitter"), 
       order = c("0", "1", "2", "3", "4", "5"),
       ylab = "GSH (nmol/mg)", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot GSSG by Ballooning grade and color by grade
p3 <- ggboxplot(data, x = "Ballooning", y = "GSSG", 
                color = "Ballooning", palette = "npg",
                order = c("0", "1", "2", "3", "4", "5"),
                ylab = "GSSG (nmol/mg)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot GSSG by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p4 <- ggline(data, x = "Ballooning", y = "GSSG", 
             color = "Black",
             #color = "Ballooning", palette = "npg",
             add = c("mean_se", "jitter"), 
             order = c("0", "1", "2", "3", "4", "5"),
             ylab = "GSSG (nmol/mg)", xlab = "Ballooning Grade")


##################################################################################################



# Box plots
# ++++++++++++++++++++
# Plot Total GSH by Ballooning grade and color by grade
p5 <- ggboxplot(data, x = "Ballooning", y = "TotalGSH", 
                color = "Ballooning", palette = "npg",
                order = c("0", "1", "2", "3", "4", "5"),
                ylab = "Total GSH (nmol/mg)", xlab = "Ballooning Grade")



# Mean plots
# ++++++++++++++++++++
# Plot Total GSH  by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p6 <- ggline(data, x = "Ballooning", y = "TotalGSH", 
             color = "Black",
             #color = "Ballooning", palette = "npg",
             add = c("mean_se", "jitter"), 
             order = c("0", "1", "2", "3", "4", "5"),
             ylab = "Total GSH (nmol/mg)", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot GSH/GSSG by Ballooning grade and color by grade
p5 <- ggboxplot(data, x = "Ballooning", y = "GSHGSSGRatio", 
                color = "Ballooning", palette = "npg",
                order = c("0", "1", "2", "3", "4", "5"),
                ylab = "GSH/GSSG", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot GSH/GSSG  by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p6 <- ggline(data, x = "Ballooning", y = "GSHGSSGRatio", 
             color = "Black",
             #color = "Ballooning", palette = "npg",
             add = c("mean_se", "jitter"), 
             order = c("0", "1", "2", "3", "4", "5"),
             ylab = "GSH/GSSG", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot Eh GSSG/2GSH by Ballooning grade and color by grade
p7 <- ggboxplot(data, x = "Ballooning", y = "EhGSSG2GSH", 
                color = "Ballooning", palette = "npg",
                order = c("0", "1", "2", "3", "4", "5"),
                ylab = "Eh GSSG/2GSH", xlab = "Ballooning Grade")


# Mean plots
# ++++++++++++++++++++
# Plot Eh GSSG/2GSH  by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p8 <- ggline(data, x = "Ballooning", y = "EhGSSG2GSH", 
             color = "Black",
             #color = "Ballooning", palette = "npg",
             add = c("mean_se", "jitter"), 
             order = c("0", "1", "2", "3", "4", "5"),
             ylab = "Eh GSSG/2GSH", xlab = "Ballooning Grade")

##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot NADH by Ballooning grade and color by grade
p9 <- ggboxplot(data, x = "Ballooning", y = "NADH", 
                color = "Ballooning", palette = "npg",
                order = c("0", "1", "2", "3", "4", "5"),
                ylab = "NADH (pmol/ug)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot NADH by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p10 <- ggline(data, x = "Ballooning", y = "NADH", 
             color = "Black",
             #color = "Ballooning", palette = "npg",
             add = c("mean_se", "jitter"), 
             order = c("0", "1", "2", "3", "4", "5"),
             ylab = "NADH (pmol/ug)", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot NADP by Ballooning grade and color by grade
p11 <- ggboxplot(data, x = "Ballooning", y = "NADP", 
                color = "Ballooning", palette = "npg",
                order = c("0", "1", "2", "3", "4", "5"),
                ylab = "NADP (pmol/ug)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot NADP by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p12 <- ggline(data, x = "Ballooning", y = "NADP", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "NADP (pmol/ug)", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot NADPH by Ballooning grade and color by grade
p13 <- ggboxplot(data, x = "Ballooning", y = "NADPH", 
                color = "Ballooning", palette = "npg",
                order = c("0", "1", "2", "3", "4", "5"),
                ylab = "NADPH (pmol/ug)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot NADPH by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p14 <- ggline(data, x = "Ballooning", y = "NADPH", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "NADPH (pmol/ug)", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot NADP/NADPH by Ballooning grade and color by grade
p15 <- ggboxplot(data, x = "Ballooning", y = "NADPNADPHRatio", 
                 color = "Ballooning", palette = "npg",
                 order = c("0", "1", "2", "3", "4", "5"),
                 ylab = "NADP/NADPH", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot NADP/NADPH by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p16 <- ggline(data, x = "Ballooning", y = "NADPNADPHRatio", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "NADP/NADPH", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot AST by Ballooning grade and color by grade
p17 <- ggboxplot(data, x = "Ballooning", y = "AST", 
                 color = "Ballooning", palette = "npg",
                 order = c("0", "1", "2", "3", "4", "5"),
                 ylab = "AST (U/L)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot AST by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p18 <- ggline(data, x = "Ballooning", y = "AST", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "AST (U/L)", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot ALT by Ballooning grade and color by grade
p19 <- ggboxplot(data, x = "Ballooning", y = "ALT", 
                 color = "Ballooning", palette = "npg",
                 order = c("0", "1", "2", "3", "4", "5"),
                 ylab = "ALT (U/L)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot ALT by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p20 <- ggline(data, x = "Ballooning", y = "ALT", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "ALT (U/L)", xlab = "Ballooning Grade")


##################################################################################################


# Box plots
# ++++++++++++++++++++
# Plot BUN by Ballooning grade and color by grade
p21 <- ggboxplot(data, x = "Ballooning", y = "BUN", 
                 color = "Ballooning", palette = "npg",
                 order = c("0", "1", "2", "3", "4", "5"),
                 ylab = "BUN (mg/dL)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot BUN by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p22 <- ggline(data, x = "Ballooning", y = "BUN", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "BUN (mg/dL)", xlab = "Ballooning Grade")


##################################################################################################

# Box plots
# ++++++++++++++++++++
# Plot Glucose by Ballooning grade and color by grade
p23 <- ggboxplot(data, x = "Ballooning", y = "Glucose", 
                 color = "Ballooning", palette = "npg",
                 order = c("0", "1", "2", "3", "4", "5"),
                 ylab = "Glucose (mg/dL)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot Glucose by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p24 <- ggline(data, x = "Ballooning", y = "Glucose", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "Glucose (mg/dL)", xlab = "Ballooning Grade")


##################################################################################################

# Box plots
# ++++++++++++++++++++
# Plot Liver Weight by Ballooning grade and color by grade
p25 <- ggboxplot(data, x = "Ballooning", y = "LiverWeight", 
                 color = "Ballooning", palette = "npg",
                 order = c("0", "1", "2", "3", "4", "5"),
                 ylab = "Liver Weight (g)", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot Liver Weight by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p26 <- ggline(data, x = "Ballooning", y = "LiverWeight", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "Liver Weight (g)", xlab = "Ballooning Grade")

##################################################################################################

# Box plots
# ++++++++++++++++++++
# Plot Steatosis by Ballooning grade and color by grade
p27 <- ggboxplot(data, x = "Ballooning", y = "Steatosis", 
                 color = "Ballooning", palette = "npg",
                 order = c("0", "1", "2", "3", "4", "5"),
                 ylab = "Steatosis", xlab = "Ballooning Grade")

# Mean plots
# ++++++++++++++++++++
# Plot Steatosis by Ballooning grade
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
p28 <- ggline(data, x = "Ballooning", y = "Steatosis", 
              color = "Black",
              #color = "Ballooning", palette = "npg",
              add = c("mean_se", "jitter"), 
              order = c("0", "1", "2", "3", "4", "5"),
              ylab = "Steatosis", xlab = "Ballooning Grade")


library(ggplot2)
pdf(file = "Ballooning-plots.pdf")
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, ncol = 1, nrow = 1)
dev.off()


#individual KW tests for reference
KW_GSH <- kruskal.test(GSH ~ Ballooning, data = data)
KW_GSSG <- kruskal.test(GSSG ~ Ballooning, data = data)
KW_TotalGSH <- kruskal.test(TotalGSH ~ Ballooning, data = data)
KW_GSHGSSGRatio <- kruskal.test(GSHGSSGRatio ~ Ballooning, data = data)
KW_Eh <- kruskal.test(EhGSSG2GSH ~ Ballooning, data = data)
KW_NADH <- kruskal.test(NADH ~ Ballooning, data = data)
KW_NADP <- kruskal.test(NADP ~ Ballooning, data = data)
KW_NADPH <- kruskal.test(NADPH ~ Ballooning, data = data)
KW_NADPNADPHRatio <- kruskal.test(NADPNADPHRatio ~ Ballooning, data = data)
KW_AST <- kruskal.test(AST ~ Ballooning, data = data)
KW_ALT <- kruskal.test(ALT ~ Ballooning, data = data)
KW_BUN <- kruskal.test(BUN ~ Ballooning, data = data)
KW_Glucose <- kruskal.test(Glucose ~ Ballooning, data = data)
KW_LiverWeight <- kruskal.test(LiverWeight ~ Ballooning, data = data)
KW_Steatosis <- kruskal.test(Steatosis ~ Ballooning, data = data)


#to get them all in one print out
KWtests = list(KW_GSH, KW_GSSG, KW_TotalGSH, KW_GSHGSSGRatio, KW_Eh, KW_NADH, KW_NADP, KW_NADPH, KW_NADPNADPHRatio, KW_ALT, KW_AST, KW_BUN, KW_Glucose, KW_LiverWeight, KW_Steatosis)
KWtests
KWtests_corrected <- data.frame(unlist(KWtests))
#setwd
write.table(KWtests_corrected, "~/Rqtl2-Glutathione-Genetics/correlations/KWtests_Ballooning.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)




