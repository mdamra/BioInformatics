library(flowCore)
library(dplyr)
library(tidyverse)

#This code import batch flow files and sets gates (in code as blob function)

ff = read.FCS("PBMC + A301_ACH2__1_019.fcs")
ff = autocomp(ff)
ff = doTransform(ff, cols = 8:23)
niceNames = function(ff) {
  des = parameters(ff)$desc
  idx = which(!is.na(des))

  colnames(ff)[idx] = des[idx]

  ff
}
ff2 = niceNames(ff)
bb_4 = blob.boundary(ff2, parameters = c("CD4", "CD8"),
                     location = bx(c(4000, 0)), height = .1)
bb_8 = blob.boundary(ff2, parameters = c("CD4", "CD8"),
                     location = bx(c(0, 20000)), height = .1)
bb_dn = blob.boundary(ff2, parameters = c("CD4", "CD8"),
                      location = bx(c(0, 0)), height = .1)
lines(bb_4, lwd = 3, col = 'dodgerblue2')
lines(bb_8, lwd = 3, col = 'indianred2')
lines(bb_dn, lwd = 3, col = 'black')
bb_4 = inflate.contour(blob = bb_4, dist = .25)
bb_8 = inflate.contour(blob = bb_8, dist = .25)
bb_dn = inflate.contour(blob = bb_dn, dist = .25)

# make gates from the polygons
gate_4 = polygonGate(.gate = bb_4)
gate_8 = polygonGate(.gate = bb_8)
gate_dn = polygonGate(.gate = bb_dn)
lines(bb_4, lwd = 3, col = 'red')
lines(bb_8, lwd = 3, col = 'indianred2')
lines(bb_dn, lwd = 3, col = 'black')
n_cd4 = nrow(Subset(ff2, gate_4))
n_cd8 = nrow(Subset(ff2, gate_8))
n_dn = nrow(Subset(ff2, gate_dn))
pplot(ff2, c("CD4","CD8"), xlim = c(-1, 5), ylim = c(-1, 5))
lines(bb_4, lwd = 3)
lines(bb_8, lwd = 3)
lines(bb_dn, lwd = 3)
n_tot = nrow(ff2)
text(x = mean(bb_4[, 1]) + 1, y = mean(bb_4[, 2]), labels = sprintf("%.1f%%", 100 * n_cd4 / n_tot), pos = 4, cex = 2)
text(x = mean(bb_8[, 1]), y = mean(bb_8[, 2]) + 1, labels = sprintf("%.1f%%", 100 * n_cd8 / n_tot), pos = 4, cex = 2)
text(x = mean(bb_dn[, 1]), y = mean(bb_dn[, 2]) + 1, labels = sprintf("%.1f%%", 100 * n_dn / n_tot), pos = 4, cex = 2)
n_cd4 = nrow(Subset(ff2, gate_4))
n_cd8 = nrow(Subset(ff2, gate_8))
n_dn = nrow(Subset(ff2, gate_dn))
lines(bb_4, lwd = 3)
lines(bb_8, lwd = 3)
lines(bb_dn, lwd = 3)
n_tot = nrow(ff2)
text(x = mean(bb_4[, 1]) + 1, y = mean(bb_4[, 2]),
     labels = sprintf("%.1f%%", 100 * n_cd4 / n_tot), pos = 4, cex = 2)
text(x = mean(bb_8[, 1]), y = mean(bb_8[, 2]) + 1,
     labels = sprintf("%.1f%%", 100 * n_cd8 / n_tot), pos = 4, cex = 2)
text(x = mean(bb_dn[, 1]), y = mean(bb_dn[, 2]) + 1,
     labels = sprintf("%.1f%%", 100 * n_dn / n_tot), pos = 4, cex = 2)
thresh_cd4 = ibx(deGate(ff2,channel = "CD4"))
thresh_cd8 = ibx(deGate(ff2,channel = "CD8"))
xline(bx(thresh_cd4), lwd = 3)
yline(bx(thresh_cd8), lwd = 3)
qg = quadGate(filterId="quad", "CD4" = bx(thresh_cd4), "CD8" = bx(thresh_cd8))
sapply(split(ff2, qg), FUN = function(x){nrow(x)}) / n_tot
thresh_cd3 = ibx(deGate(f = ff2, channel = "CD3"))
thresh_dump = ibx(deGate(f = ff2, channel = "DUMP"))
xline(bx(thresh_dump), lwd = 3)
yline(bx(thresh_cd3), lwd = 3)
qg_tcell = quadGate(filterId="tcell", "DUMP" = bx(thresh_dump), "CD3" = bx(thresh_cd3))
tcell = split(ff2, qg_tcell)$`DUMP-CD3+`
n_tcell = nrow(tcell)
pplot(ff2, c("CD4", "CD8"), xlim = c(-1, 5), ylim = c(-1, 5), main = "Ungated")
xline(bx(thresh_cd4), lwd = 3)
yline(bx(thresh_cd8), lwd = 3)

thresh_cd4 = ibx(deGate(tcell,channel = "CD4"))
thresh_cd8 = ibx(deGate(tcell,channel = "CD8"))
xline(bx(thresh_cd4), lwd = 3)
yline(bx(thresh_cd8), lwd = 3)
qg2 = quadGate(filterId="quad", "CD4" = bx(thresh_cd4), "CD8" = bx(thresh_cd8))
sapply(split(ff2, qg), FUN = function(x){nrow(x)}) / n_tot
sapply(split(tcell, qg2), FUN = function(x){nrow(x)}) / n_tcell













