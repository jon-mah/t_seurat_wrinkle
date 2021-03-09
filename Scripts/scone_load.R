# scone_load.R

library(SeuratData)
library(scone)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
InstallData('pbmcsca')
data('pbmcsca')

gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
memory.limit()
memory.limit(size=156000)
my_scone = readRDS(file='scone_5000.rds')

pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)


bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)

points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
       pch = 1, col = "blue", cex = 1)
points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
       pch = 1, col = "blue", cex = 1.5)

arrows(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][1],
       bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][2],
       bp_obj[1,][1],
       bp_obj[1,][2],
       lty = 2, lwd = 2)
