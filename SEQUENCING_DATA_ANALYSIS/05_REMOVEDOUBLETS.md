# REMOVE DOUBLETS USING DOUBLETFINDER

## SEURAT CLUSTERING BEFORE RUNNING DOUBLETFINDER
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl

```{R}
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

##------------------------------------
## Load CellBender retained data
##------------------------------------
load("CLOCKSEQ_P07_P56_WT_HU_KO_CELLBENDER_RETAINED.RData")
# dt.cb

##------------------------------------
## Drop the bad samples based on CellBender
# YL_P07_WT_04
# YL_P07_HU_17
# YL_P07_KO_06
dt.cb.removed <- subset(dt.cb, subset = SampleIndex %in% c("YL_P07_WT_04", "YL_P07_HU_17", "YL_P07_KO_06"), invert = FALSE) # 33
dt.cb.cleaned <- subset(dt.cb, subset = SampleIndex %in% c("YL_P07_WT_04", "YL_P07_HU_17", "YL_P07_KO_06"), invert = TRUE) # 140,605
dt.cb.cleaned.meta <- as.data.frame(dt.cb.cleaned@meta.data)

##------------------------------------
## fetch count table
DefaultAssay(dt.cb.cleaned) <- "RNA"
cbCleanedCounts <- GetAssayData(dt.cb.cleaned, slot = "counts")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = cbCleanedCounts, project = "CLOCK_SEQ")

##------------------------------------
## update meta data
seuObj@meta.data <- merge(seuObj@meta.data, dt.cb.cleaned.meta, by = "row.names")
row.names(seuObj@meta.data) <- seuObj@meta.data$Row.names
seuObj@meta.data$Row.names <- NULL

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)

ggsave(filename = "SEURAT_CLOCK_SEQ_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_CLOCK_SEQ_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)


##------------------------------------
## Data Scaling
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_CLOCK_SEQ_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_CLOCK_SEQ_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)


##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "CLOCKSEQ_P07_P56_WT_HU_KO_CELLBENDER_RETAINED_CLUST.RData")
```


## DOUBLETFINDER
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl

```{R}
##------------------------------------
## DOUBLET FINDER
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("CLOCKSEQ_P07_P56_WT_HU_KO_CELLBENDER_RETAINED_CLUST.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

# ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
# sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
# bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seuObj@meta.data$CellType2) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "DF.classifications_0.25_0.005_14060", sct = FALSE)

table(seuObjDF$SampleIndex, seuObjDF$DF.classifications_0.25_0.005_14060)        
#                Doublet Singlet
#   YL_P07_HU_03     262    9796
#   YL_P07_HU_05     365    7842
#   YL_P07_HU_12     400    8352
#   YL_P07_KO_01     869    6201
#   YL_P07_KO_09     283    7985
#   YL_P07_KO_18     287    3994
#   YL_P07_WT_02     514    4875
#   YL_P07_WT_11     297    8633
#   YL_P07_WT_16     344    3550
#   YL_P56_HU_02     720    5580
#   YL_P56_HU_05     958    4183
#   YL_P56_HU_09    1070    8285
#   YL_P56_HU_11     957    4610
#   YL_P56_KO_03     540    5593
#   YL_P56_KO_04     932    5250
#   YL_P56_KO_07     694    2303
#   YL_P56_KO_13     725    4871
#   YL_P56_WT_01    1036    4933
#   YL_P56_WT_08     783    7477
#   YL_P56_WT_10    1189    7247
#   YL_P56_WT_12     835    4985

table(seuObjDF$SampleIndex, seuObjDF$DF.classifications_0.25_0.005_11883)             
#                Doublet Singlet
#   YL_P07_HU_03    9796     262
#   YL_P07_HU_05    2087    6120
#   YL_P07_HU_12       0    8752
#   YL_P07_KO_01       0    7070
#   YL_P07_KO_09       0    8268
#   YL_P07_KO_18       0    4281
#   YL_P07_WT_02       0    5389
#   YL_P07_WT_11       0    8930
#   YL_P07_WT_16       0    3894
#   YL_P56_HU_02       0    6300
#   YL_P56_HU_05       0    5141
#   YL_P56_HU_09       0    9355
#   YL_P56_HU_11       0    5567
#   YL_P56_KO_03       0    6133
#   YL_P56_KO_04       0    6182
#   YL_P56_KO_07       0    2997
#   YL_P56_KO_13       0    5596
#   YL_P56_WT_01       0    5969
#   YL_P56_WT_08       0    8260
#   YL_P56_WT_10       0    8436
#   YL_P56_WT_12       0    5820


##------------------------------------
## Save RData
save(seuObjDF, file = "CLOCKSEQ_P07_P56_WT_HU_KO_CELLBENDER_RETAINED_CLUST_DOUBLETFINDER.RData")


seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_14060 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_11883 == "Singlet",]

## Before Doublets Removal
table(seuMeta$SampleIndex)
# YL_P07_HU_03 YL_P07_HU_05 YL_P07_HU_12 YL_P07_KO_01 YL_P07_KO_09 YL_P07_KO_18 
#        10058         8207         8752         7070         8268         4281 
# YL_P07_WT_02 YL_P07_WT_11 YL_P07_WT_16 YL_P56_HU_02 YL_P56_HU_05 YL_P56_HU_09 
#         5389         8930         3894         6300         5141         9355 
# YL_P56_HU_11 YL_P56_KO_03 YL_P56_KO_04 YL_P56_KO_07 YL_P56_KO_13 YL_P56_WT_01 
#         5567         6133         6182         2997         5596         5969 
# YL_P56_WT_08 YL_P56_WT_10 YL_P56_WT_12 
#         8260         8436         5820

## After Doublets Removal
table(seuMetaFilt$SampleIndex)
# YL_P07_HU_05 YL_P07_HU_12 YL_P07_KO_01 YL_P07_KO_09 YL_P07_KO_18 YL_P07_WT_02 
#         5755         8352         6201         7985         3994         4875 
# YL_P07_WT_11 YL_P07_WT_16 YL_P56_HU_02 YL_P56_HU_05 YL_P56_HU_09 YL_P56_HU_11 
#         8633         3550         5580         4183         8285         4610 
# YL_P56_KO_03 YL_P56_KO_04 YL_P56_KO_07 YL_P56_KO_13 YL_P56_WT_01 YL_P56_WT_08 
#         5593         5250         2303         4871         4933         7477 
# YL_P56_WT_10 YL_P56_WT_12 
#         7247         4985

save(seuMeta, seuMetaFilt, file = "CLOCKSEQ_P07_P56_WT_HU_KO_CELLBENDER_RETAINED_CLUST_DOUBLETFINDER_META.RData")


seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))
df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

table(seuObjDF$DoubletFinder)
# Discarded  Retained 
#     25943    114662

```

-----
