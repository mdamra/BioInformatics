library(dbplyr)
library(dplyr)
library(corrplot)
library(tidyverse)
library(tidyselect)
library(tidyr)
#import dataset and set to dataframe variable
dataR =read.csv("R.csv")
R = data.matrix(dataR)
rownames(R) = dataR$X
R = R[,2:9]
colnames(R) = c("L-Glutamic acid",	"AAL",	"LCA",	"UDA",
                "STL",	"ABA",	"MPA",	"Pyruvic Acid")
rownames(R) = c("CA-DNA", "CA-RNA ",
                "Intact HIV DNA" ,"Defective HIV DNA", "Hypermut HIV DNA")


dataP =read.csv("P.csv")
P = data.matrix(dataP)
rownames(P) = dataP$X
P = P[,2:9]
colnames(P) = c("L-Glutamic acid",	"AAL",	"LCA",	"UDA",
                "STL",	"ABA",	"MPA",	"Pyruvic Acid")
rownames(P) = c("CA-DNA", "CA-RNA ",
                "Intact HIV DNA" ,"Defective HIV DNA", "Hypermut HIV DNA")



col3 <- colorRampPalette(c("blue", "white", "red"))
corrplot(t(R), method="circle", is.corr=FALSE,
         sig.level = c(.001, .01, .05),
         pch.col = "white",p.mat = t(P),insig = "label_sig" ,
         tl.col = "Black",
         tl.srt = 45, col = col3(100))

