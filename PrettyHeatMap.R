library(tidyverse)
library(pheatmap)

#Create pretty heatmap with annotated key
#Read CSV file to variable
data = read.csv("JO.csv")
#Assign rownames and adjust df dimensions
rownames(data) = data[[1]]
colnames(data) = c("Name","RMG1",	"RMG1-KO",	"RMG1 VIRAL(CONTROL)",	"RMG1 shARIDIA-KD",	"OVCA429",	"OVCA429-KO",	"RMG1 GALNT7",	"RMG1 KO GALNT7",	"RMG1 GCNT3",	"RMG1 KO GCNT3",	"RMG1 KO CNTRL",	"RMG1 CNTRL")
data = data[1:45, 2:13]
roworder = c("Calsepa","PHA(E)","AAL","AOL","LCA","PSA","TJA-II","UEA_I","EEL","TxLC_",
            "ECA","RCA120","DBA","HPA","VVA","WFA","PTL_I","GSL_I_A4","GSL_I_B4","SBA",
            "GSL-II","STL","PWM","DSA","LEL","PHA(L)","WGA","ConA","LTL","GNA","HHL",
           "NPA","UDA","MAH","MAL_I","SNA","SSA","TJA-I","ACG","ACA","MPA","PNA","ABA",
           "Jacalin","BPL")
data = data[roworder,]

#Calculate Z-Score for each x in row
r = apply(data[,1:6], 1, function(x){(x - mean(x))/sd(x)})
r = t(r)
r1 = apply(data[,7:12], 1, function(x){(x - mean(x))/sd(x)})
r1 = t(r1)


r2 = cbind(r,r1)

#Label Sample Types
colannot = data.frame(Sample = c(rep("RMA GALNT7",2),
                                 rep("RMA GCNT3",2),
                                 rep("CONTROL",2)
                                 ))

#Label Specificity
rowannot = data.frame(Specificity = c(rep("Bisecting GlcNAC",2),
                                      rep("Fucose",6),
                                      rep("Gal, Fucose", 1),
                                      rep("Gal, GalNac",1),
                                      rep("Gal, LacNac",2),
                                      rep("GalNac", 4),
                                      rep("GalNac, Fucose, Gal", 1),
                                      rep("GalNac, Gal",1),
                                      rep("Gal",1),
                                      rep("GalNac, Gal",1),
                                      rep("GlcNac",2),
                                      rep("GlcNac, Gal, Mannose",1),
                                      rep("GlcNac, LacNac",3),
                                      rep("GlcNac, Sialic Acid",1),
                                      rep("Mannose",1),
                                      rep("Mannose, Fucose",1),
                                      rep("Mannose, Gal",1),
                                      rep("Mannose, Gal, GlcNac",2),
                                      rep("Mannose, GlcNac",1),
                                      rep("Sialic Acid",5),
                                      rep("Sialic Acid, Gal",1),
                                      rep("T-Antigen",3),
                                      rep("T-Antigen, GlcNAc",1),
                                      rep("T-Antigen, GlcNac, Mannose",1),
                                      rep("T-Antigen, LacNAc, Gal",1)))



#Assign annotations to columns and rows
rownames(rowannot) = rownames(r)
row.names(colannot) = colnames(r)

#Plot heatmap with annotations and breaks
pheatmap(r2, cluster_rows = F,cluster_cols = F,gaps_row = c(2,8,12,27,33,39), cellwidth = 20,
         annotation_row = rowannot, annotation_col = colannot)



