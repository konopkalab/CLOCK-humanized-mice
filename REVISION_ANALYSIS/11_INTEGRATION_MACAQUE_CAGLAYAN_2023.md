# SEURAT INTEGRATION, MACAQUE, CAGLAYAN 2023

## Filter CellBender+DoubletFinder Datasets
```{R}
## START
## Load R Libraries
rm(list = ls())
.libPaths("/project/RESOURCES/FROM_HOME/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(dittoSeq)
library(scCustomize)
library(ShinyCell)
set.seed(10)


```


## Seurat Preprocessing
```{R}
## FUNCTION TO FILTER CELLBARCODES (CB-RETAINED) AND GENES (REF-PROTEIN-CODING)
filterCounts <- function(cbPath, dfPath, samplePrefix) {
print(samplePrefix)

## load cb data
cbData <- Read_CellBender_h5_Mat(cbPath, use.names = TRUE, unique.features = TRUE, h5_group_name = NULL, feature_slot_name = "features")
print(dim(cbData))
colnames(cbData) <- paste(samplePrefix, colnames(cbData), sep = "_")
cbData <- as.data.frame(cbData)

## load df data
load(dfPath)
dfRetained <- row.names(seuObj_doublets_removed@meta.data)

## filter cb data
cbdfData <- cbData[, colnames(cbData) %in% dfRetained]
print(dim(cbdfData))

## return
return(cbdfData)

## LOAD CELLBENDER H5 FILE
cbdf_m1 <- filterCounts("/project/C_CELLBENDER_KONOPKA/CELLBENDER_M1_Out_filtered.h5",
"/project/D_DOUBLETFINDER_KONOPKA/M1_Seurat_DF.RData",
"M1")

cbdf_m2 <- filterCounts("/project/C_CELLBENDER_KONOPKA/CELLBENDER_M2_Out_filtered.h5",
"/project/D_DOUBLETFINDER_KONOPKA/M2_Seurat_DF.RData",
"M2")

cbdf_m3 <- filterCounts("/project/C_CELLBENDER_KONOPKA/CELLBENDER_M3_Out_filtered.h5",
"/project/D_DOUBLETFINDER_KONOPKA/M3_Seurat_DF.RData",
"M3")

cbdf_m4 <- filterCounts("/project/C_CELLBENDER_KONOPKA/CELLBENDER_M4_Out_filtered.h5",
"/project/D_DOUBLETFINDER_KONOPKA/M4_Seurat_DF.RData",
"M4")

save(cbdf_m1, cbdf_m2, cbdf_m3, cbdf_m4, file = "KONOPKA_MACAQUE_CBDF_COUNTS.RData")

rm(list = ls())
load("KONOPKA_MACAQUE_CBDF_COUNTS.RData")

## FUNCTION TO RUN SEURAT PREPROCESSING PER SAMPLE
preprocess_seurat <- function(countsTable, samplePrefix, sampleSpecies) {
print(samplePrefix)
seuobj <- CreateSeuratObject(counts = countsTable)
seuobj <- PercentageFeatureSet(seuobj, pattern = "^MT-", col.name = "pMito_RNA")
seuobj$Sample <- samplePrefix
seuobj$Species <- sampleSpecies
return(seuobj)
}

seu_m1 <- preprocess_seurat(cbdf_m1, "M1", "MACAQUE")
seu_m2 <- preprocess_seurat(cbdf_m2, "M2", "MACAQUE")
seu_m3 <- preprocess_seurat(cbdf_m3, "M3", "MACAQUE")
seu_m4 <- preprocess_seurat(cbdf_m4, "M4", "MACAQUE")

seu_list <- list(M1 = seu_m1, M2 = seu_m2, M3 = seu_m3, M4 = seu_m4)
save(seu_list, file = "KONOPKA_MACAQUE_CBDF_SEURAT.RData")
}
```


## Merge & Integration
```{R}
## load list of seurat objects
rm(list = ls())
load("KONOPKA_MACAQUE_CBDF_SEURAT.RData")

## merge list of seurat objects into a single seurat object
seuObjMerged <- merge(x = seu_list$M1,
y = c(seu_list$M2,
seu_list$M3,
seu_list$M4
),
project = "KONOPKA_MACAQUE")

## run standard anlaysis workflow
seuObjMerged <- NormalizeData(seuObjMerged)
seuObjMerged <- FindVariableFeatures(seuObjMerged)
seuObjMerged <- ScaleData(seuObjMerged)
seuObjMerged <- RunPCA(seuObjMerged)

seuObjMerged <- FindNeighbors(seuObjMerged, dims = 1:30, reduction = "pca")
seuObjMerged <- FindClusters(seuObjMerged, resolution = 0.8, cluster.name = "merged_clusters")

## save merged object
save(seuObjMerged, file = "KONOPKA_MACAQUE_CBDF_SEURAT_MERGED.RData")

## perform integration
seuObjInt <- IntegrateLayers(object = seuObjMerged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
seuObjInt[["RNA"]] <- JoinLayers(seuObjInt[["RNA"]])
seuObjInt <- FindNeighbors(seuObjInt, reduction = "integrated.cca", dims = 1:30)
seuObjInt <- FindClusters(seuObjInt, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObjInt <- RunUMAP(seuObjInt, dims = 1:30, reduction = "integrated.cca", min.dist = 0.1, spread = 1, n.epochs = 1000)

## Save integrated object
save(seuObjInt, file = "KONOPKA_MACAQUE_CBDF_SEURAT_INTEGRATED.RData")
```


## Generate Shiny App
```{R}
rm(list = ls())
load("KONOPKA_MACAQUE_CBDF_SEURAT_INTEGRATED.RData")
# seuObjInt

seuObjInt[["RNA3"]] <- as(object = seuObjInt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjInt)
makeShinyApp(seuObjInt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "KONOPKA_MACAQUE/", shiny.title = "KONOPKA_MACAQUE")

##END
```


## Filter Data
```{R}
## Load Integrated Seurat Object
rm(list = ls())
load("KONOPKA_MACAQUE_CBDF_SEURAT_INTEGRATED.RData")
# seuObjInt

DefaultAssay(seuObjInt) <- "RNA"
Idents(seuObjInt) <- "RNA_snn_res.0.8"

## filtering criterion
nuhi <- 4000
pmhi <- 10

print(dim(seuObjInt))
# 19605 56135

filt_temp <- subset(seuObjInt, subset = nCount_RNA < nuhi & pMito_RNA < pmhi)
print(dim(filt_temp))
# 19605 55726

seuObjIntFilt <- subset(filt_temp, idents = c(8, 18, 21), invert = TRUE)
print(dim(seuObjIntFilt))
# 19605 55522

## Re-cluster
seuObjIntFilt$RNA_snn_res.0.2 <- NULL
seuObjIntFilt$RNA_snn_res.0.4 <- NULL
seuObjIntFilt$RNA_snn_res.0.6 <- NULL
seuObjIntFilt$RNA_snn_res.0.8 <- NULL
seuObjIntFilt$RNA_snn_res.1 <- NULL
seuObjIntFilt$RNA_snn_res.1.2 <- NULL 
seuObjIntFilt$RNA_snn_res.1.4 <- NULL
seuObjIntFilt$RNA_snn_res.1.6 <- NULL
seuObjIntFilt$RNA_snn_res.1.8 <- NULL
seuObjIntFilt$RNA_snn_res.2 <- NULL
seuObjIntFilt$seurat_clusters <- NULL

seuObjIntFilt <- FindNeighbors(seuObjIntFilt, reduction = "integrated.cca", dims = 1:30)
seuObjIntFilt <- FindClusters(seuObjIntFilt, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObjIntFilt <- RunUMAP(seuObjIntFilt, reduction = "integrated.cca", dims = 1:30, min.dist = 0.1, spread = 1, n.epochs = 1000)

seutree <- clustree(seuObjIntFilt, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "KONOPKA_MACAQUE_CBDF_SEURAT_INTEGRATED_FILT_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## Save integrated object
save(seuObjIntFilt, file = "KONOPKA_MACAQUE_CBDF_SEURAT_INTEGRATED_FILT.RData")
```


## Generate Shiny App
```{R}
rm(list = ls())
load("KONOPKA_MACAQUE_CBDF_SEURAT_INTEGRATED_FILT.RData")
# seuObjIntFilt

seuObjIntFilt[["RNA3"]] <- as(object = seuObjIntFilt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjIntFilt)
makeShinyApp(seuObjIntFilt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "KONOPKA_MACAQUE_FILT/", shiny.title = "KONOPKA_MACAQUE_FILT")

## END
```
