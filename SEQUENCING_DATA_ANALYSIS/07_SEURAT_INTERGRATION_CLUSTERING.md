# SEURAT INTEGRATION AND CLUSTERING

```{R}
###################################
## LOAD BIOHPC MODULES
## BEFORE STARTING R SESSION
## ON A TERMINAL
###################################
module add python
module load hdf5_18/1.8.17
module load R/4.0.2-gccmkl


###################################
## AFTER STARTING R SESSION
## LOAD LIBRARIES
###################################
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(WGCNA)

##-------------------------------------------------------
## LOAD DATASETS TO COMBINE
##-------------------------------------------------------
load("pathway/name.RData")
wt <- seuObjAll
rm(seuObjAll)

load("pathway/name.RData")
ko <- seuObjAll
rm(seuObjAll)

load("pathway/name.RData")
hu <- seuObjAll
rm(seuObjAll)

##-------------------------------------------------------
## DATA INTEGRATION
##-------------------------------------------------------
## list of seurat objects and pre-processing
clock.list <- c(wt, ko, hu)

clock.list <- lapply(X = clock.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

## perform integration
clock.anchors <- FindIntegrationAnchors(object.list = clock.list, dims = 1:20)
clock.combined <- IntegrateData(anchorset = clock.anchors, dims = 1:20)

## perform clustering on integrated object
DefaultAssay(clock.combined) <- "integrated"

## scale data, run pca, run jackstraw
clock.combined <- ScaleData(clock.combined, verbose = FALSE)
clock.combined <- RunPCA(clock.combined, npcs = 100, verbose = FALSE)
clock.combined <- JackStraw(clock.combined, dims = 100)
clock.combined <- ScoreJackStraw(clock.combined, dims = 1:100)

## elbow plot to determine number of PC for clustering
pEl <- ElbowPlot(object = clock.combined, ndims = 100)
ggsave(filename = "ElbowPlot.pdf" , plot = pEl, width = 8, height = 6, units = "in", dpi = 150)

## run umap and identify clusters, selpcs is the number of PC decided by using elbow plot
clock.combined <- RunUMAP(clock.combined, reduction = "pca", dims = 1:selpcs)
clock.combined <- FindNeighbors(clock.combined, reduction = "pca", dims = 1:selpcs)
clock.combined <- FindClusters(clock.combined, resolution = c(1.2))

##-------------------------------------------------------
## GENERATE BASIC PLOTS
##-------------------------------------------------------
## plot umap clusters
umap <- DimPlot(clock.combined, reduction = "umap", label = TRUE)
ggsave(filename = "umap.pdf", plot = umap, width = 8, height = 6, units = "in", dpi = 150)

## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE
genoPerCluster <- as.data.frame.matrix(table(clock.combined@meta.data$integrated_snn_res.1.2, clock.combined@meta.data$Genotype))
genoPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(genoPerCluster))), sep = "_")
write.table(ClockPerCluster, "composition.txt", row.names = F, col.names = T, quote = F, sep = "\t")
genoPerCluster2 <- melt(genoPerCluster)
colnames(genoPerCluster2) <- c("CLUSTER", "GENOTYPE", "CELLS")
composition <- ggplot(genoPerCluster2) +
  geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOTYPE), stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NULL
ggsave(filename = "name.pdf", plot = composition, width = 8, height = 4, units = "in", dpi = 150)

##-------------------------------------------------------
## SAVE INTEGRATED and CLUSTERED DATASET
##-------------------------------------------------------
save(clock.combined, file = "CLUSTERED.RData")
```

-----
