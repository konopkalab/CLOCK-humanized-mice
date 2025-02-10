
## Create Seurat Objects
```{R}
## START
## R Libraries
rm(list = ls())
.libPaths("/project/RESOURCES/R/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(rhdf5)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(dittoSeq)
library(scCustomize)
set.seed(10)

## Function to Run Seurat pipeline
process_sample <- function(sample_name, cbpath) {
print(paste("Processing ", sample_name, sep = ""))

## Load dataset from CellBender run
cbData <- Read_CellBender_h5_Mat(cbpath, use.names = TRUE, unique.features = TRUE, h5_group_name = NULL, feature_slot_name = "features")
print(dim(cbData))
colnames(cbData) <- paste(sample_name, colnames(cbData), sep = "_")
cbDataTemp <- as.data.frame(cbData)
## cbData2 <- cbData[rowSums(cbData[]) > 0,] removes genes with 0 row sums
cbData <- cbDataTemp[, colSums(cbDataTemp[]) > 0] ## removes cells with 0 col sums

## Create seurat object
seuObj <- CreateSeuratObject(counts = cbData, project = sample_name)
seuObj[["Sample_Name"]] <- sample_name
print(dim(seuObj))

## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^MT-")

## Set identities
Idents(seuObj) <- "Sample_Name"

## Visualize Data QC
p1 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 4, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = paste(sample_name, "_Seurat_QC_1.pdf", sep = ""), plot = p1, width = 12, height = 4, units = "in", dpi = 150)
ggsave(filename = paste(sample_name, "_Seurat_QC_2.pdf", sep = ""), plot = p2, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = paste(sample_name, "_Seurat_QC_3.pdf", sep = ""), plot = p3, width = 6, height = 4, units = "in", dpi = 150)

## Data Normalization
seuObj <- NormalizeData(seuObj)

## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## Plot variable features with and without labels
p5 <- VariableFeaturePlot(seuObj)
p6 <- LabelPoints(plot = p5, points = top90, repel = TRUE)
ggsave(filename = paste(sample_name, "_Seurat_Var_Genes_1.pdf", sep = ""), plot = p5, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = paste(sample_name, "_Seurat_Var_Genes_2.pdf", sep = ""), plot = p6, width = 8, height = 4, units = "in", dpi = 150)

## Data Scaling
## all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p7 <- ElbowPlot(seuObj, ndims = 100)
p8 <- JackStrawPlot(seuObj, dims = 1:100)
ggsave(filename = paste(sample_name, "_Seurat_PCA_1.pdf", sep = ""), plot = p7, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = paste(sample_name, "_Seurat_PCA_2.pdf", sep = ""), plot = p8, width = 12, height = 6, units = "in", dpi = 150)

## Identify significant PCs
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
if (length(pcScores) < 1) {
selpcs <- 30
} else {
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
}

## Data Clustering
seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.8))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

## Save Seurat object RData
save(seuObj, file = paste(sample_name, "_Seurat.RData", sep = ""))
print(paste("Processing ", sample_name, " completed", sep = ""))
}

## Call the function above for each sample from Chimpanzee, Human and Macaque samples across Jorstad et al., 2023 and Caglayan et al., 2023 datasets.
## END
```


## Run DoubletFinder
```{R}
## START
## R Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

## Function to Run DoubletFinder
process_sample <- function(sample_name) {
print(paste("Processing ", sample_name, sep = ""))
seuobj_file <- paste(sample_name, "_Seurat.RData", sep = "")

## Load cellbender seurat object
load(seuobj_file)
# seuObj

## Identify significant PCs
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)

## pK Identification (no ground-truth)
sweep.res.list_str <- paramSweep(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)

pdf(paste(sample_name, "pK_Sweep.pdf", sep = "_"), width=7, height=6)
plot(bcmvn_str$pK, bcmvn_str$MeanBC, type = "l", xlab = "pK", ylab = "BCMVN", main = "BCMVN vs pK")
dev.off()

selpk <- 0.005

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8)
nExp_poi <- round(0.15*nrow(seuObj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seuObjDF <- doubletFinder(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seuObjDF <- doubletFinder(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = colnames(seuObjDF@meta.data)[grepl("^pANN", colnames(seuObjDF@meta.data))], sct = FALSE)
save(seuObjDF, selpcs, selpk, nExp_poi.adj, file = paste(sample_name, "_Seurat_DF_Temp.RData", sep = ""))

## Remove doublets
seuMeta <- as.data.frame(seuObjDF@meta.data)
dfscores <- colnames(seuMeta)[grepl("^DF.classifications_", colnames(seuMeta))]
seuMetaFiltTemp <- seuMeta[seuMeta[,dfscores[[1]]] == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp[,dfscores[[2]]] == "Singlet",]

seuMetaFilt$DoubletFinder <- rep("Singlet", nrow(seuMetaFilt))
df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)
seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Doublets"

## Save seurat object as RData before removing Doublets
save(seuObjDF, file = paste(sample_name, "_Seurat_DF_Classified.RData", sep = ""))
print(table(seuObjDF$DoubletFinder))

## seurat object after dropping doublets
seuObj_doublets_removed <- subset(seuObjDF, subset = DoubletFinder == "Singlet")

## Save doublets removed seurat object as RData
save(seuObj_doublets_removed, file = paste(sample_name, "_Seurat_DF.RData", sep = ""))
}

## Call the function above for each sample from Chimpanzee, Human and Macaque samples across Jorstad et al., 2023 and Caglayan et al., 2023 datasets.
## END
```
