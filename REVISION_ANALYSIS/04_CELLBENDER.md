
## Prepare raw CellRanger data as input for CellBender
```{R}
## START

## R Libraries
rm(list = ls())
.libPaths("/project/RESOURCES/R/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(DropletQC)
library(Matrix)
set.seed(10)

prepareCB <- function(raw_bc, sample) {
DATA_CB <- Read10X_h5(raw_bc)
DATA_CB_RNA <- DATA_CB
DATA_CB_GENES <- as.data.frame(rownames(DATA_CB_RNA))
colnames(DATA_CB_GENES) <- "GeneID"
DATA_CB_GENES$GeneName <- DATA_CB_GENES$GeneID
DATA_CB_GENES$Category <- rep("Gene Expression", nrow(DATA_CB_GENES))
DATA_CB_BARCODES <- as.data.frame(colnames(DATA_CB_RNA))
colnames(DATA_CB_BARCODES) <- "Barcode"
dir.create(sample)
setwd(sample)

write.table(DATA_CB_GENES, "genes.tsv", col.names = F, sep = '\t', quote = F, row.names = F)
write.table(DATA_CB_BARCODES, "barcodes.tsv", col.names = F, sep = '\t', quote = F, row.names = F)
writeMM(DATA_CB_RNA, "matrix.mtx")

setwd("./../")
}

prepareCB("NW_TX0017-12_raw_feature_bc_matrix.h5", "NW_TX0017-12")

## Repeat for all Chimpanzee, Human and Macaque Samples.
## Repeat for Jorstad et al., 2023 and Caglayan et al., 2023 datasets.
## END
```

## Run CellBender
```{bash}
## START

echo "===> NW_TX0017-12"
/project/RESOURCES/TOOLS/CONDACB/cellbendergpuak/bin/cellbender remove-background \
--cuda \
--input NW_TX0017-12/ \
--output CELLBENDER_NW_TX0017_12_Out.h5 \
--expected-cells 8000 \
--total-droplets-included 50000

## Repeat for all Chimpanzee, Human and Macaque Samples.
## Repeat for Jorstad et al., 2023 and Caglayan et al., 2023 datasets.
## END
```

