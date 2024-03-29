library(survival)
library(survminer)
library(dplyr)
library(tidyverse)

#This is a log rank test of inhibitors on Jlat and Jurkat cells with an HIV-Promotor. 
#The cells were activated with aCD3/aCD28 and PMA/Ionomycin and commited to a "TIME-TO-REBOUND" experiment.
#Cells were then challenged with inhibotor molecules to inhibit the HIV-1 Promotor after activation. The batch flow .fcs files
#where analyzed and gated using the FlowGating R code within this repository then GFP MFI is recorded. 

#The STI dataset contains over 3,000 inhibitor molecules that were challenged by a automated robot cytometer.
#This code analyzes the MFI and builds a survival plot and calculates PVAL and FDR of the successful inhibitors.


#import and shape data
read = read.csv('G2.csv')
df1 = read[1:24,]

#for loop performs 3 actions:
#1. Calculates median of each parameter then runs if statement to label timepoint as "above" or "below median.
#2. fits divided data into survival model .
#3. Builds name and pastes into plot then saves jpg from "dev.off".

for (i in 2:length(df1)){
  df = df1[,c(1,i)]
  df <- df %>%
    mutate(Median = ifelse(df[,2] > median(df[,2]), "Above", "Below"))
  fit = survfit(Surv(df$Time.to.viral.Rebound) ~ Median, data = df)
  mypath <- paste(colnames(df[2]), ".jpg")
  jpeg(mypath)
  print(ggsurvplot(fit, data = df, pval = T, risk.table = F,conf.int = F,
             surv.median.line = "hv") +
    ggtitle(colnames(df[2])) +
    ylab(""))
  dev.off()
}

#for loop performs 3 actions:
#1. Calculates median of each parameter then runs if statement to label timepoint as "above" or "below median.
#2. fits divided data into survival model then retrieves p-value from model matrix into df "tab".
#3. adds calculated FDR from dataset from adjusted p-value (q-value)command.
tab = data.frame("Species" = character(), "Pval" = numeric())
colnames(tab) = c("Species", "Pval")

for (i in 2:length(df1)){
  df = df1[,c(1,i)]
  df <- df %>%
    mutate(Median = ifelse(df[,2] > median(df[,2]), "Above", "Below"))
  fit = survfit(Surv(df$Time.to.viral.Rebound) ~ Median, data = df)
  record = data.frame(Species = colnames(df[2]),
                      Pval = surv_pvalue(fit)$pval,
                      stringsAsFactors = FALSE)
  tab = rbind(tab, record)

}
#calculate FDR (adjusted p-value (q-value))
tab  = tab %>% mutate(FDR = p.adjust(tab$Pval, method = "BH"))
write.csv(tab, "G2Analysis.csv")











