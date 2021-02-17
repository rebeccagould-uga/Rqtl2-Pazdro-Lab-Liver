#Plotting Spearmans Rank (rho) correlations in Correlation Matrix Heat Map
#R01 GSH Liver qRT-PCR data - Socs1 GSH/GSSG Experiment

#Hprt REFERENCE GENE*****************



#Removed all outliers
#Saved data in CSV format
#GSH/GSSG data for the Liver using adjusted GSH/GSSG values 

#plot info reference: https://www.displayr.com/how-to-create-a-correlation-matrix-in-r/

#This section provides a simple function for formatting a correlation matrix into a table with 4 columns containing :
#Column 1 : row names (variable 1 for the correlation test)
#Column 2 : column names (variable 2 for the correlation test)
#Column 3 : the correlation coefficients
#Column 4 : the p-values of the correlations
#The custom function below can be used :

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

###################################################

#Low GSH/GSSG group
  rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/PCR-Socs1/Hprt-ReferenceGene/Low-GSHGSSGRatio-Targ-Ref.csv", check.names = FALSE)
  head(rawdata, 6)
  rownames(rawdata) <- rawdata$Mouse
  #remove column 1-3 (Mouse ID and sex)
  data = subset(rawdata, select = -c(1,2))
  library("Hmisc")
  
  low_data_correlations = rcorr(as.matrix(data), type = c("spearman"))

  low_data_coefficients = low_data_correlations$r
  low_data_p = low_data_correlations$P
  
  flattenCorrMatrix(low_data_coefficients, low_data_p)
  

#High GSH/GSSG group
  rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/PCR-Socs1/Hprt-ReferenceGene/High-GSHGSSGRatio-Targ-Ref.csv", check.names = FALSE)
  head(rawdata, 6)
  rownames(rawdata) <- rawdata$Mouse
  #remove column 1-3 (Mouse ID and sex)
  data = subset(rawdata, select = -c(1,2))
  library("Hmisc")
  
  high_data_correlations = rcorr(as.matrix(data), type = c("spearman"))

  high_data_coefficients = high_data_correlations$r
  high_data_p = high_data_correlations$P
  
  flattenCorrMatrix(high_data_coefficients, high_data_p)
  

#All mice
  rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/PCR-Socs1/Hprt-ReferenceGene/GSHGSSGRatio-Targ-Ref.csv", check.names = FALSE)
  head(rawdata, 6)
  rownames(rawdata) <- rawdata$Mouse
  #remove column 1-3 (Mouse ID and sex)
  data = subset(rawdata, select = -c(1,2))
  library("Hmisc")
  
  all_data_correlations = rcorr(as.matrix(data), type = c("spearman"))
  
  all_data_coefficients = all_data_correlations$r
  all_data_p = all_data_correlations$P
  
  flattenCorrMatrix(all_data_coefficients, all_data_p)




####################################################
## exporting data
####################################################
#set working directory

low_pvalues <- as.data.frame(low_data_p)
low_correlations <- as.data.frame(low_data_coefficients)

high_pvalues <- as.data.frame(high_data_p)
high_correlations <- as.data.frame(high_data_coefficients)

all_pvalues <- as.data.frame(all_data_p)
all_correlations <- as.data.frame(all_data_coefficients)

library(writexl)

write_xlsx(list("All Correlations" = all_correlations, 
                "All P Values" = all_pvalues,
                "Low Correlations" = low_correlations, 
                "Low P Values" = low_pvalues,
                "High Correlations" = high_correlations, 
                "High P Values" = high_pvalues),
                "Hprt-correlations and p values.xlsx")

####################################################
## assessing significance between groups
## Wilcoxan Rank Sum test (Mann Whitney test)
####################################################

#High vs Low GSH/GSSG
  my_data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/PCR-Socs1/Hprt-ReferenceGene/HighvsLow.csv", check.names = FALSE)

  LowGroup <- my_data[1:51,2]
  HighGroup <- my_data[52:100,2]
  
  res <- wilcox.test(LowGroup,HighGroup,exact = FALSE)
  res
  res$p.value
  