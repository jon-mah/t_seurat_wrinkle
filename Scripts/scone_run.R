
library(SeuratData)
library(scone)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
InstallData('pbmcsca')
data('pbmcsca')

# cc <- c(brewer.pal(9, "Set1"))
# batch = factor(pbmcsca$Method)
# qc = pbmcsca$nUMI
# qc <- sapply(qc, as.numeric)
# o = order(qc)[order(batch[order(qc)])]
# barplot(qc[o], col=cc[batch][o], 
#         border=cc[batch][o], main="UMI counts")
# legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)
# 
# 
# expr <- GetAssayData(pbmcsca)
# # experiment_name <- '10x Chromium (v2)'
# # col_select <- pbmcsca$Method == experiment_name
# # col_select_1 <- pbmcsca$Experiment == 'pbmc2'
# # expr = expr[, col_select]
# # expr <- expr[which(apply(expr > 0, 1, any)), ]
# # pbmcsca_small <- pbmcsca[, (pbmcsca$Method == experiment_name)]

# Define Housekeeping object
data(housekeeping)
hk = intersect(housekeeping$V1,rownames(pbmcsca))

# Mean log10(x+1) expression
# Assumed False Negatives




# num_reads = quantile(expr[expr > 0])[4]
# num_cells = 0.25*ncol(expr)
# is_common = rowSums(expr >= num_reads ) >= num_cells

# Metric-based Filtering
# mfilt = metric_sample_filter(as.matrix(expr),
#                             gene_filter = is_common,
#                             pos_controls = rownames(expr) %in% hk,
#                             zcut = 3, mixture = FALSE,
#                             plot = TRUE)



#mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
#goodDat = expr[, mfilt]
# data <- subset(pbmcsca, idents='pbmc1')

# Redefine memory limit size for Rstudio on Windows
memory.limit(size=300000)

# Grab expression data
expr = GetAssayData(pbmcsca)
n=100 # cell number, default to 1000 for bugtesting
set.seed(1)
sampled = sample(colnames(expr), size=n, replace=F)
expr = expr[, sampled]
expr = expr[which(apply(expr > 0, 1, any)), ]
data = pbmcsca[, sampled]
print(dim(data))

num_reads = quantile(expr[expr > 0])[4]
num_cells = 0.25*ncol(expr)
is_common = rowSums(expr >= num_reads ) >= num_cells

# Metric-based Filtering
mfilt = metric_sample_filter(as.matrix(expr),
                             gene_filter = is_common,
                             hard_nreads=2,
                             # pos_controls = rownames(expr) %in% hk,
                             zcut = 3, mixture = FALSE,
                             plot = TRUE)
mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
goodDat = expr[, mfilt]
# data <- subset(pbmcsca, idents='pbmc1')

print(dim(data))
goodDat = as.matrix(expr)
qc_metrics = c('nCount_RNA', 'nFeature_RNA', 'nGene', 'nUMI', 'percent.mito')
qc = as.matrix(data[[qc_metrics]])
qc = apply(qc, 2, as.numeric)
ppq = scale(qc[,apply(qc,2,sd) > 0],center = TRUE,scale = TRUE)
negcon = intersect(rownames(goodDat),hk)




my_scone <- SconeExperiment(goodDat,
                            qc=ppq,
                            bio=as.factor(data$CellType),
                            batch=as.factor(data$Method),
                            negcon_ruv=rownames(goodDat) %in% negcon,
)



EFF_FN = function (ei)
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

## ----- Scaling Argument -----

scaling=list(none=identity, # Identity - do nothing
             
             eff = EFF_FN, # User-defined function
             
             sum = SUM_FN, # SCONE library wrappers...
             tmm = TMM_FN, 
             uq = UQ_FN,
             fq = FQT_FN,
             deseq = DESEQ_FN)



# my_scone <- scone(my_scone,
#                 scaling=scaling,
#                 k_qc=3, k_ruv = 0,
#                 adjust_bio="no",
#                 run=FALSE)



BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution
rm(data)
rm(expr)
rm(goodDat)
rm(pbmcsca)
gc()

my_scone <- scone(my_scone,
                  scaling=scaling,
                  run=TRUE,
                  # k_ruv=0,
                  eval_kclust = 2:6,
                  stratified_pam = FALSE,
                  return_norm = "in_memory",
                  zero = "postadjust")


pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)




saveRDS(my_scone, file="scone_100.rds")





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


