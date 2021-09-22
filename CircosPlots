Create Circos Plots

library(dplyr)
library(dbplyr)
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

#import data and convert to matrix
#Shape Corr data
r = read.csv("IGG1.csv")
r1 = r[,1:3]
r1 = r1 %>% spread(key = To, value = R, fill = 0)
rownames(r1) = r1$From
r1 = r1[,2:11]
mat = data.matrix(r1)

#Shape significance data (secondardy table for second track)
r2 = r[,c("From","To","P")]
r2 = r2 %>% spread(key = To, value = P, fill = 0)
rownames(r2) = r2$From
r2 = r2[,2:11]
matP = data.matrix(r2)

# Order elements in inner circle
col.order = c("sCD163","LBP","MCP-2","MIP-1a","SDF-1a", "TNF-a","CP_SCORE","SISCORE" ,
              "TOTSC_AS", "TOTSC_VS")
row.order = c("Agalactosylated","Bisected","G0","G0F","G0FB","G1F","G2F","LBP","sCD163",                  "Term.Galactosylated","Tot.Galactosylated","Tot.Sialyated"  )
mat = mat[row.order, col.order]

#convert matrix into vectored data.frame and set directional relations for arrows
df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)

#filter NA..if zeros are present after compile, replace with filter(value != 0)
df1 <- df %>% filter(!is.na(value))

#set elements color codes
grid.col = c("Agalactosylated"="grey","Bisected"="grey","G0"="grey","G0F"="grey",               "G0FB"="grey","G1F"="grey","G2F"="grey","sCD163"="grey",
             "Term.Galactosylated"="grey",  "Tot.Galactosylated"="grey",
             "Tot.Sialyated"="grey","CP_SCORE"="grey","CP_SCORE"="grey", "MCP-2"="grey",
             "MIP-1a"="grey",   "sCD163"="grey",   "SDF-1a"="grey", "SISCORE"="grey",
             "TNF-a"="grey","TOTSC_AS"="grey", "TOTSC_VS"="grey", "LBP"="grey",
             "MCP-2"="grey",    "MIP-1a"="grey",     "SDF-1a"="grey",
             "SISCORE"="grey",  "TNF-a"="grey",    "TOTSC_AS"="grey", "TOTSC_VS" ="grey"  )

#draw plot and customize
col_text = c("orange", "orange", "orange","orange","orange","orange","light blue",
"light blue","light blue","light blue","green","red","green","green", "red")

#Initialize circos plot and set as data frame for annotation
circos.par(start.degree = 197,canvas.ylim = c(-1.0,1.0))
cdm_res = chordDiagram(df1,order = c(colnames(mat), rownames(mat)),
             col = ifelse(df1$value > 0, "red", "blue"),
             directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.length = 0.1,symmetric = T,
             link.arr.type = "big.arrow", diffHeight = -uh(1, "mm"),
             scale = F,grid.col = grid.col ,annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = uh(3, "mm"),
                                      track.margin = c(uh(3, "mm"), 0)))
head(cdm_res)

#plot circos
circos.clear()
circos.par(start.degree = 200)
chordDiagram(df1,order = c(colnames(mat), rownames(mat)),
             col = ifelse(df1$value > 0, "red", "blue"),
             directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.length = 0.1,symmetric = T,
             link.arr.type = "big.arrow", diffHeight = -uh(1, "mm"),
             scale = F,grid.col = grid.col ,annotationTrack = c("grid"),
             preAllocateTracks = list(track.height = uh(1, "mm"),
                                      track.margin = c(uh(1, "mm"), 0)))

#For loop that expands label range to fit
for (i in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.3, sector.index = i, track.index = 2) }

#add labels for element track

circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.5, niceFacing = TRUE)
}, bg.border = NA, track.height = uh(5,"mm"))

#set y-coordinates
ylim = get.cell.meta.data("ylim", sector.index = rownames(mat)[1], track.index = 1)
y1 = ylim[1] + (ylim[2] - ylim[1])*0.4
y2 = ylim[2]
abs_max = quantile(abs(c(matP) - 0.5), 0.95, na.rm = TRUE)

#creates categorical highlights(outter most track)
highlight.sector(c("Agalactosylated","Bisected","G0","G0F","G0FB"),
                 track.index = 1, col = "purple",facing = "bending.inside",
                 text = "Pro-Inflammatory Glycans", cex = 0.8,
                 padding = c(0, 0, 1.5, 0),
                 text.col = "black", niceFacing = TRUE)
highlight.sector(c("sCD163","LBP","MCP-2","MIP-1a","SDF-1a", "TNF-a"),
                  track.index = 1, col = "orange",
                 text = "Plasma Inflammatory Markers", cex = 0.8,
                 padding = c(0, 0, 1.5, 0),facing = "bending.inside",
                 text.col = "black", niceFacing = TRUE)
highlight.sector(c("Term.Galactosylated","Tot.Galactosylated","Tot.Sialyated","G1F","G2F"), track.index = 1, col = "green",
                 padding = c(0, 0, 1.5, 0),facing = "bending.inside",
                 text = "Anti-Inflammatory Glycans", cex = 0.8,
                 text.col = "black", niceFacing = TRUE)
highlight.sector(c("CP_SCORE","SISCORE" ,
                   "TOTSC_AS", "TOTSC_VS" ),
                 track.index = 1, col = "light blue",
                 padding = c(0, 0, 1.5, 0),facing = "bending.inside",
                 text = "CVD Scores", cex = 0.8,
                 text.col = "black", niceFacing = TRUE)

# set color pallete/gamut for significance meter
col_fun = colorRamp2(c(0.5 - abs_max, 0.5, 0.5 + abs_max),
                     c("red", "white", "blue"))

# annotates second dataframe data to categories and elements
for (i in seq_len(nrow(cdm_res))) {
  if (cdm_res$value1[i] != 0) {

    #annotates row data from cdm_res
    circos.rect(cdm_res[i, "x1"], y1, cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]),
                y1 + (y2 - y1)*0.45,
                col = col_fun(matP[cdm_res$rn[i], cdm_res$cn[i]]),
                border = col_fun(matP[cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$rn[i], track.index = 1)

    #annotates column from cdm_res
    circos.rect(cdm_res[i, "x2"], y1, cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]),
                y1 + (y2 - y1)*0.45,
                col = col_fun(matP[cdm_res$rn[i], cdm_res$cn[i]]),
                border = col_fun(matP[cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$cn[i], track.index = 1)
  }
}

#clear plot
circos.clear()

