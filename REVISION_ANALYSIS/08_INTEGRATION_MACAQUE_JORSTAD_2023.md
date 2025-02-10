
## Filter CellBender+DoubletFinder Datasets
```{R}
## START
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
}

## LOAD CELLBENDER H5 FILE
cbdf_nwtx_0023_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_7_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_7_Seurat_DF.RData",
                              "NW_TX0023_7")

cbdf_nwtx_0023_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_8_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_8_Seurat_DF.RData",
                              "NW_TX0023_8")

cbdf_nwtx_0023_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_9_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_9_Seurat_DF.RData",
                              "NW_TX0023_9")

cbdf_nwtx_0023_10 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_10_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_10_Seurat_DF.RData",
                              "NW_TX0023_10")

cbdf_nwtx_0023_11 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_11_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_11_Seurat_DF.RData",
                              "NW_TX0023_11")

cbdf_nwtx_0023_12 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_12_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_12_Seurat_DF.RData",
                              "NW_TX0023_12")

cbdf_nwtx_0023_13 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_13_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_13_Seurat_DF.RData",
                              "NW_TX0023_13")

cbdf_nwtx_0023_14 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0023_14_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0023_14_Seurat_DF.RData",
                              "NW_TX0023_14")

cbdf_nwtx_0024_2 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_2_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_2_Seurat_DF.RData",
                              "NW_TX0024_2")

cbdf_nwtx_0024_3 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_3_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_3_Seurat_DF.RData",
                              "NW_TX0024_3")

cbdf_nwtx_0024_4 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_4_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_4_Seurat_DF.RData",
                              "NW_TX0024_4")

cbdf_nwtx_0024_5 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_5_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_5_Seurat_DF.RData",
                              "NW_TX0024_5")

cbdf_nwtx_0032_4 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0032_4_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0032_4_Seurat_DF.RData",
                              "NW_TX0032_4")

cbdf_nwtx_0032_5 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0032_5_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0032_5_Seurat_DF.RData",
                              "NW_TX0032_5")

cbdf_nwtx_0032_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0032_6_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0032_6_Seurat_DF.RData",
                              "NW_TX0032_6")

cbdf_nwtx_0032_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0032_7_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0032_7_Seurat_DF.RData",
                              "NW_TX0032_7")

cbdf_nwtx_0032_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0032_8_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0032_8_Seurat_DF.RData",
                              "NW_TX0032_8")

cbdf_nwtx_0032_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0032_9_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0032_9_Seurat_DF.RData",
                              "NW_TX0032_9")

cbdf_nwtx_0033_16 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0033_16_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0033_16_Seurat_DF.RData",
                              "NW_TX0033_16")

cbdf_nwtx_0053_5 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0053_5_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0053_5_Seurat_DF.RData",
                              "NW_TX0053_5")

cbdf_nwtx_0053_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0053_6_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0053_6_Seurat_DF.RData",
                              "NW_TX0053_6")

cbdf_nwtx_0053_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0053_7_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0053_7_Seurat_DF.RData",
                              "NW_TX0053_7")

cbdf_nwtx_0073_10 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0073_10_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0073_10_Seurat_DF.RData",
                              "NW_TX0073_10")

cbdf_nwtx_0074_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0074_6_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0074_6_Seurat_DF.RData",
                              "NW_TX0074_6")

cbdf_nwtx_0074_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0074_7_Out_filtered.h5", 
                              "/project/D_DOUBLETFINDER_ALLEN/NW_TX0074_7_Seurat_DF.RData",
                              "NW_TX0074_7")

save(cbdf_nwtx_0023_7, cbdf_nwtx_0023_8, cbdf_nwtx_0023_9, cbdf_nwtx_0023_10, cbdf_nwtx_0023_11, cbdf_nwtx_0023_12, cbdf_nwtx_0023_13, cbdf_nwtx_0023_14, 
     cbdf_nwtx_0024_2, cbdf_nwtx_0024_3, cbdf_nwtx_0024_4, cbdf_nwtx_0024_5, 
     cbdf_nwtx_0032_4, cbdf_nwtx_0032_5, cbdf_nwtx_0032_6, cbdf_nwtx_0032_7, cbdf_nwtx_0032_8, cbdf_nwtx_0032_9, 
     cbdf_nwtx_0033_16, 
     cbdf_nwtx_0053_5, cbdf_nwtx_0053_6, cbdf_nwtx_0053_7, 
     cbdf_nwtx_0073_10, 
     cbdf_nwtx_0074_6, cbdf_nwtx_0074_7, 
file = "ALLEN_MACAQUE_CBDF_COUNTS.RData")

```

## Seurat Preprocessing
```{R}
rm(list = ls())
load("ALLEN_MACAQUE_CBDF_COUNTS.RData")

## FUNCTION TO RUN SEURAT PREPROCESSING PER SAMPLE
preprocess_seurat <- function(countsTable, samplePrefix, sampleSpecies) {

    # countsTable <- cbdf_h1
    # samplePrefix <- "C1"
    # sampleSpecies <- "CHIMPANZEE"

    print(samplePrefix)
    seuobj <- CreateSeuratObject(counts = countsTable)
    seuobj <- PercentageFeatureSet(seuobj, pattern = "^MT-", col.name = "pMito_RNA")
    seuobj$Sample <- samplePrefix
    seuobj$Species <- sampleSpecies

    return(seuobj)
}

seu_nwtx_0023_7 <- preprocess_seurat(cbdf_nwtx_0023_7, "NW_TX0023_7", "MACAQUE")
seu_nwtx_0023_8 <- preprocess_seurat(cbdf_nwtx_0023_8, "NW_TX0023_8", "MACAQUE")
seu_nwtx_0023_9 <- preprocess_seurat(cbdf_nwtx_0023_9, "NW_TX0023_9", "MACAQUE")
seu_nwtx_0023_10 <- preprocess_seurat(cbdf_nwtx_0023_10, "NW_TX0023_10", "MACAQUE")
seu_nwtx_0023_11 <- preprocess_seurat(cbdf_nwtx_0023_11, "NW_TX0023_11", "MACAQUE")
seu_nwtx_0023_12 <- preprocess_seurat(cbdf_nwtx_0023_12, "NW_TX0023_12", "MACAQUE")
seu_nwtx_0023_13 <- preprocess_seurat(cbdf_nwtx_0023_13, "NW_TX0023_13", "MACAQUE")
seu_nwtx_0023_14 <- preprocess_seurat(cbdf_nwtx_0023_14, "NW_TX0023_14", "MACAQUE")
seu_nwtx_0024_2 <- preprocess_seurat(cbdf_nwtx_0024_2, "NW_TX0024_2", "MACAQUE")
seu_nwtx_0024_3 <- preprocess_seurat(cbdf_nwtx_0024_3, "NW_TX0024_3", "MACAQUE")
seu_nwtx_0024_4 <- preprocess_seurat(cbdf_nwtx_0024_4, "NW_TX0024_4", "MACAQUE")
seu_nwtx_0024_5 <- preprocess_seurat(cbdf_nwtx_0024_5, "NW_TX0024_5", "MACAQUE")
seu_nwtx_0032_4 <- preprocess_seurat(cbdf_nwtx_0032_4, "NW_TX0032_4", "MACAQUE")
seu_nwtx_0032_5 <- preprocess_seurat(cbdf_nwtx_0032_5, "NW_TX0032_5", "MACAQUE")
seu_nwtx_0032_6 <- preprocess_seurat(cbdf_nwtx_0032_6, "NW_TX0032_6", "MACAQUE")
seu_nwtx_0032_7 <- preprocess_seurat(cbdf_nwtx_0032_7, "NW_TX0032_7", "MACAQUE")
seu_nwtx_0032_8 <- preprocess_seurat(cbdf_nwtx_0032_8, "NW_TX0032_8", "MACAQUE")
seu_nwtx_0032_9 <- preprocess_seurat(cbdf_nwtx_0032_9, "NW_TX0032_9", "MACAQUE")
seu_nwtx_0033_16 <- preprocess_seurat(cbdf_nwtx_0033_16, "NW_TX0033_16", "MACAQUE")
seu_nwtx_0053_5 <- preprocess_seurat(cbdf_nwtx_0053_5, "NW_TX0053_5", "MACAQUE")
seu_nwtx_0053_6 <- preprocess_seurat(cbdf_nwtx_0053_6, "NW_TX0053_6", "MACAQUE")
seu_nwtx_0053_7 <- preprocess_seurat(cbdf_nwtx_0053_7, "NW_TX0053_7", "MACAQUE")
seu_nwtx_0073_10 <- preprocess_seurat(cbdf_nwtx_0073_10, "NW_TX0073_10", "MACAQUE")
seu_nwtx_0074_6 <- preprocess_seurat(cbdf_nwtx_0074_6, "NW_TX0074_6", "MACAQUE")
seu_nwtx_0074_7 <- preprocess_seurat(cbdf_nwtx_0074_7, "NW_TX0074_7", "MACAQUE")

seu_list <- list("NW_TX0023_7" = seu_nwtx_0023_7, 
                 "NW_TX0023_8" = seu_nwtx_0023_8, 
                 "NW_TX0023_9" = seu_nwtx_0023_9, 
                 "NW_TX0023_10" = seu_nwtx_0023_10, 
                 "NW_TX0023_11" = seu_nwtx_0023_11, 
                 "NW_TX0023_12" = seu_nwtx_0023_12, 
                 "NW_TX0023_13" = seu_nwtx_0023_13, 
                 "NW_TX0023_14" = seu_nwtx_0023_14, 
                 "NW_TX0024_2" = seu_nwtx_0024_2, 
                 "NW_TX0024_3" = seu_nwtx_0024_3, 
                 "NW_TX0024_4" = seu_nwtx_0024_4, 
                 "NW_TX0024_5" = seu_nwtx_0024_5, 
                 "NW_TX0032_4" = seu_nwtx_0032_4, 
                 "NW_TX0032_5" = seu_nwtx_0032_5, 
                 "NW_TX0032_6" = seu_nwtx_0032_6, 
                 "NW_TX0032_7" = seu_nwtx_0032_7, 
                 "NW_TX0032_8" = seu_nwtx_0032_8, 
                 "NW_TX0032_9" = seu_nwtx_0032_9, 
                 "NW_TX0033_16" = seu_nwtx_0033_16, 
                 "NW_TX0053_5" = seu_nwtx_0053_5, 
                 "NW_TX0053_6" = seu_nwtx_0053_6, 
                 "NW_TX0053_7" = seu_nwtx_0053_7, 
                 "NW_TX0073_10" = seu_nwtx_0073_10, 
                 "NW_TX0074_6" = seu_nwtx_0074_6,
                 "NW_TX0074_7" = seu_nwtx_0074_7)
save(seu_list, file = "ALLEN_MACAQUE_CBDF_SEURAT.RData")

```

## Merge & Integration
```{R}
rm(list = ls())
load("ALLEN_MACAQUE_CBDF_SEURAT.RData")

## merge list of seurat objects into a single seurat object
seuObjMerged <- merge(x = seu_list$NW_TX0023_7, 
                      y = c(seu_list$NW_TX0023_8, 
                            seu_list$NW_TX0023_9, 
                            seu_list$NW_TX0023_10, 
                            seu_list$NW_TX0023_11, 
                            seu_list$NW_TX0023_12, 
                            seu_list$NW_TX0023_13, 
                            seu_list$NW_TX0023_14,
                            seu_list$NW_TX0024_2,
                            seu_list$NW_TX0024_3,
                            seu_list$NW_TX0024_4,
                            seu_list$NW_TX0024_5,
                            seu_list$NW_TX0032_4,
                            seu_list$NW_TX0032_5,
                            seu_list$NW_TX0032_6,
                            seu_list$NW_TX0032_7,
                            seu_list$NW_TX0032_8,
                            seu_list$NW_TX0032_9,
                            seu_list$NW_TX0033_16,
                            seu_list$NW_TX0053_5,
                            seu_list$NW_TX0053_6,
                            seu_list$NW_TX0053_7,
                            seu_list$NW_TX0073_10,
                            seu_list$NW_TX0074_6,
                            seu_list$NW_TX0074_7
                            ), 
                      project = "ALLEN_MACAQUE")

## run standard anlaysis workflow
seuObjMerged <- NormalizeData(seuObjMerged)
seuObjMerged <- FindVariableFeatures(seuObjMerged)
seuObjMerged <- ScaleData(seuObjMerged)
seuObjMerged <- RunPCA(seuObjMerged)

seuObjMerged <- FindNeighbors(seuObjMerged, dims = 1:30, reduction = "pca")
seuObjMerged <- FindClusters(seuObjMerged, resolution = 0.8, cluster.name = "merged_clusters")

## save merged object
save(seuObjMerged, file = "ALLEN_MACAQUE_CBDF_SEURAT_MERGED.RData")

## perform integration
seuObjInt <- IntegrateLayers(object = seuObjMerged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
seuObjInt[["RNA"]] <- JoinLayers(seuObjInt[["RNA"]])
seuObjInt <- FindNeighbors(seuObjInt, reduction = "integrated.cca", dims = 1:30)
seuObjInt <- FindClusters(seuObjInt, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObjInt <- RunUMAP(seuObjInt, dims = 1:30, reduction = "integrated.cca", min.dist = 0.1, spread = 1, n.epochs = 1000)

## Save integrated object
save(seuObjInt, file = "ALLEN_MACAQUE_CBDF_SEURAT_INTEGRATED.RData")

```

## Generate Shiny App
```{R}
rm(list = ls())
load("ALLEN_MACAQUE_CBDF_SEURAT_INTEGRATED.RData")
# seuObjInt

seuObjInt[["RNA3"]] <- as(object = seuObjInt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjInt)
makeShinyApp(seuObjInt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "ALLEN_MACAQUE/", shiny.title = "ALLEN_MACAQUE") 

## END
```

## Filter Data
```{R}
## Load Integrated Seurat Object
rm(list = ls())
load("ALLEN_MACAQUE_CBDF_SEURAT_INTEGRATED.RData")
# seuObjInt

DefaultAssay(seuObjInt) <- "RNA"
Idents(seuObjInt) <- "RNA_snn_res.0.8"

## filtering criterion
nuhi <- 10000
pmhi <- 5

print(dim(seuObjInt))
# 19623 199352

filt_temp <- subset(seuObjInt, subset = nCount_RNA < nuhi & pMito_RNA < pmhi)
print(dim(filt_temp))
# 19605 198976

seuObjIntFilt <- subset(filt_temp, idents = c(23, 30, 32, 33), invert = TRUE)
print(dim(seuObjIntFilt))
# 19605 196499

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
ggsave(filename = "ALLEN_MACAQUE_CBDF_SEURAT_INTEGRATED_FILT_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## Save integrated object
save(seuObjIntFilt, file = "ALLEN_MACAQUE_CBDF_SEURAT_INTEGRATED_FILT.RData")
```


## Generate Shiny App
```{R}
rm(list = ls())
load("ALLEN_MACAQUE_CBDF_SEURAT_INTEGRATED_FILT.RData")
# seuObjIntFilt

seuObjIntFilt[["RNA3"]] <- as(object = seuObjIntFilt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjIntFilt)
makeShinyApp(seuObjIntFilt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "ALLEN_MACAQUE_FILT/", shiny.title = "ALLEN_MACAQUE_FILT")

## END
```

