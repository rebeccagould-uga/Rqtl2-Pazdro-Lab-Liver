#Gathering descriptive statisitcs on data
#R01 GSH Liver Data


#Removed all outliers
#Saved data in CSV format
#Glutathione data for the Liver and using Adjusted GSSG values 


#load my data into the session using read.csv function
rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data.csv")
head(rawdata, 6)

rownames(rawdata) <- rawdata$Mouse
#checking pheno file
rawdata[1:10,]

#remove column 1-3 (Mouse ID, sex, and generation)
data = subset(rawdata, select = -c(1,2))
data$Sex[data$Sex == "M"] <- 1
data$Sex[data$Sex == "F"] <- 2


##########

#to get mean by gender
stats <- aggregate(data,by=list(sex=data$Sex), mean, median, na.rm=TRUE)

#overall means
sapply(data, mean, na.rm=TRUE)

#descriptive statistics
library(Hmisc)
describe(data)
# n, nmiss, unique, mean, 5,10,25,50,75,90,95th percentiles
# 5 lowest and 5 highest scores


############

#descriptive statistics
library(pastecs)
statistics <- stat.desc(data)

# nbr.val, nbr.null, nbr.na, min max, range, sum,
# median, mean, SE.mean, CI.mean, var, std.dev, coef.var

library(psych)
statistics_sex <- describeBy(data, data$Sex)



#save excel sheet
library(writexl)
write_xlsx(list("statistics" = statistics,
                "by sex" = statistics_sex),
           "descriptive-statistics.xlsx")





