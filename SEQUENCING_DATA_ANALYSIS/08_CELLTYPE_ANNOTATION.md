# CELLTYPE ANNOTATION

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

#load saved dataset
load("CLUSTERED..RData")

##-------------------------------------------------------
## IDENTIFY CLUSTER MARKERS
##-------------------------------------------------------
seuMarkers <- FindAllMarkers(clock.combined, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)

#Generate a new dataframe that contains marker genes of each cluster
table <- read.table("composition.txt", sep = "\t", header = T)
tab <- seuMarkers
gene_exp <- length(unique(tab$gene))
tab <- tab[tab$p_val_adj <= 0.05 & tab$pct.1 >= 0.3,]
tab <- tab[c(6,7)]
colnames(tab) <- c("Cluster", "Gene")
all <- merge(table, tab, by = "Cluster")
tab <- all[c(9,10)]
colnames(tab) <- c("DEFINITION","Gene")
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes <- as.data.frame(table(tab$DEFINITION))

#load marker genes of each cell type in Yao et al., 2021, which could be downloaded from the supplementary results of the paper
biccn_master <- read.table("Allen_mouseCortex_snRNASeq_annotation.csv", sep = ",", header = T)
biccn_master_up <- biccn_master[biccn_master$direction == "up",]
biccn_master_up2 <- biccn_master_up[,c("genes", "cl1_cluster_label")]
colnames(biccn_master_up2) <- c("gene", "celltype")
biccn_list <- split(biccn_master_up2, biccn_master_up2$celltype)
biccn_list2 <- lapply(biccn_list, function(x) { x[1:50,] } )
GeneSets <- biccn_list2

# or load marker genes of each cell type in Hrvatin et al., 2018, which could be downloaded from the supplementary results of the paper
hrvatin.markers <- read.table("Greenberg_VisualCortex_Data_QC_Filt_Norm_PCA_Clust_DEG.txt", sep = "\t", header = T)
hrvatinAllSig <- hrvatin.markers[hrvatin.markers$p_val_adj <= 0.05,]
hrvatinAllSig2 <- hrvatinAllSig[hrvatinAllSig$pct.1 >= 0.5,]
hrvatinAllSig3 <- hrvatinAllSig2[c(7,6)]
vcGenes <- list(Astro = hrvatinAllSig3[hrvatinAllSig3$cluster == "Astro",],
                Endo_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Endo_1",],
                Endo_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Endo_2",],
                Macrophage = hrvatinAllSig3[hrvatinAllSig3$cluster == "Macrophage",],
                Micro_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Micro_1",],
                Micro_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Micro_2",],
                Neu_ExcL23 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL23",],
                Neu_ExcL4 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL4",],
                Neu_ExcL5_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL5_1",],
                Neu_ExcL5_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL5_2",],
                Neu_ExcL5_3 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL5_3",],
                Neu_ExcL6 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL6",],
                Neu_Int_Cck = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Cck",],
                Neu_Int_Npy = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Npy",],
                Neu_Int_Pv = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Pv",],
                Neu_Int_Sst_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Sst_1",],
                Neu_Int_Sst_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Sst_2",],
                Neu_Int_Vip = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Vip",],
                Olig_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_1",],
                Olig_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_2",],
                Olig_3 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_3",],
                Olig_4 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_4",],
                Olig_5 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_5",],
                Olig_6 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_6",],
                Olig_7 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_7",],
                OPC_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "OPC_1",],
                OPC_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "OPC_2",],
                Pericyte = hrvatinAllSig3[hrvatinAllSig3$cluster == "Pericyte",],
                SM_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "SM_1",],
                SM_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "SM_2",])
GeneSets <- vcGenes

#Re-arrange data
for(i in 1:length(GeneSets))
{
  colnames(GeneSets[[i]])[1] <- "Gene"
}

ln <- length(GeneSets)
cl <- length(Genes$Var1)
TEMP <- list()
INT <- list()
for (i in 1:ln)
{
  TEMP[[i]] <- tab[tab$Gene %in% GeneSets[[i]]$Gene,]
  INT[[i]] <- as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT) <- names(GeneSets)
names(TEMP) <- names(GeneSets)

######################### Create the matrix for each GeneSet
NROWS <- sapply(GeneSets, nrow)

#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
  INT[[i]]$b <- NROWS[[i]] - INT[[i]]$Freq
  INT[[i]]$c <- Genes$Freq - INT[[i]]$Freq
  INT[[i]]$d <- 10000 - Genes$Freq - nrow(GeneSets[[i]])
}

######################### Function for Fisher's exact test
RunFisher <- function(row, alt = 'greater', cnf = 0.85) 
{
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value),
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

######################### Run Fisher's exact test
FisherMat=list()
for (i in 1:length(INT))
{
  FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
  rownames(FisherMat[[i]]) <- INT[[i]]$Var1
  FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat) <- names(INT)
#save(FisherMat, TEMP, file = "SCN_SEURAT_FILT_NORM_PCA_CLUST_FISHER_OUTPUT_HRVATIN.RData")


######################### Arrange a matrix for P-val
tmp <- list()
FisherP <- matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)

for (i in 1:length(INT))
{
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
  FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

######################### Arrange a matrix for OR
tmp <- list()
FisherOR <- matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
  FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

######################### Compute adjusted P-val
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj > 0.05] <- 1
FisherOR[FisherOR < 1] <- 0

df <- -log10(FisherAdj)
df <- t(df)



######################### Plot
pdf("annotation.pdf",width=21, height=16, pointsize=23)
par(mar = c(7, 9, 2, 2))
LabelMat <- paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df,
               xLabels = colnames(df),
               yLabels = rownames(df),
               colorLabels =FALSE,
               colors=colorRampPalette(c("white", "red"))(50),
               #textMatrix=LabelMat,
               setStdMargins = FALSE,
               cex.text = 0.5,
               xLabelsAngle = 90)
dev.off()
```

-----