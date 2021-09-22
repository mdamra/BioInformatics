#importing lectin data and rendering a heatmap and calculating signficance.


library(tidyverse)
library(dply)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(ggplot2)

mydata = read.csv("data.csv")
lectins = mydata$Name
rownames(mydata) = mydata$Name
mydata = mydata[,-c(1)]
mydata$Name = lectins
mydata = t(mydata)

model = lm(aov(X1~X2,data = mydata))

ggplot(mydata, aes(X2, X1, color = "blue")) + geom_point() + geom_smooth(method = "lm", se = F, col = "red")



mydata = as.numeric(mydata)

mydata = mydata[,2:5]
mydata1 = t(mydata)
colnames(mydata) = mydata[1,]
mydata = mydata[2:5,]

mydata = read.csv("/Users/zeitgeist/Desktop/lectin/TimTissue/TIMCOLON-LECTIN.csv")
mydata = mydata[,-c(1)]
rownames(mydata) = mydata$Name
mydata = mydata[,-c(1)]
mydata = aggregate(mydata[,2:47], list(mydata$GROUP),mean)

mydata = as.numeric(mydata)
rownames(mydata) = mydata[,1]
colnames(mydata) = mydata[1,]

mydata = mydata[2:5,]

mydata1 = as.data.frame(mydata1)
mydata = t(mydata)

aggSIG = aggregate(sig[, 3:7], list(sig$Group), mean)
rownames(aggSIG) = aggSIG[,1]
aggSIG = aggSIG[,-c(1)]
aggSIG = t(aggSIG)

y = mydata[,1:6]
cols.cor <- cor(mydata, use = "pairwise.complete.obs", method = "spearman")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(mydata), use = "pairwise.complete.obs", method = "spearman")

pheatmap(
  mydata, scale = "row", cluster_rows = T, cluster_cols = F, angle_col=90,
  #clustering_distance_cols = as.dist(1 - cols.cor),
  #clustering_distance_rows = as.dist(1 - rows.cor)
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~P-VAL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#t.test(b,a, data = df)$p.value
lectin = factor(mydata$Name)
sink("Output.txt")
mydata$Group <- factor(mydata$Group, levels=c("Y","S","PP"))
for (i in 3:ncol(mydata)) {
  formula <- as.formula(paste(colnames(mydata)[i], " ~ lectins", sep=""))
  model <- aov(formula, data = mydata)
  print(model)
  cat("\n-----\n\n")
  cat(colnames(mydata)[i])
  cat("\n")
  print(summary(model))
}
sink()

sink("Output.txt")
#My Code
sink()

#TrustTest
model = aov(lm(SNA~Name, mydata))
model
summary(model)[[1]][["Pr(>F)"]]

library(ggplot2)
model = lm(UEA_I~Group, data = mydata)
ggplot(model, aes(Group, UEA_I, col=Group)) + geom_point()+ geom_smooth(stat = 'smooth', color = 'Red')


#select significant

sig = mydata[,c("Sample","Group","AAL","TJA.I","ACG", "TxLC_I", "STL" )]


data = read.csv("/Users/zeitgeist/Desktop/LAB/lectin/TimTissue/TIMCOLON-LECTIN.csv")
regDF = data[,-c(1)]
regDF = aggregate(regDF[,1:46], list(regDF$GROUP),mean)
regDF = regDF[,-c(2)]
regDF
colnames(regDF)
