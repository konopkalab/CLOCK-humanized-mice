# UPDATE RAW COUNT TABLES

## WT SAMPLES
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
####################
## LOAD LIBRARIES ##
####################
rm(list = ls())
library(Seurat)
library(Matrix)
library(SeuratDisk)


####################
## REF GENES      ##
####################
refGenes <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", sep = "\t", header = T)
refGenesSymbol <- as.data.frame(unique(sort(refGenes$GeneSymbol)))
colnames(refGenesSymbol) <- "genes"
row.names(refGenesSymbol) <- refGenesSymbol$genes
print(dim(refGenesSymbol))
# 21967     1


####################
## LOOP OVER ALL  ##
####################
myFiles <- list.files(path = "./COUNTS_DIR_MOUSE/", pattern = "*.gz")
myLib <- unlist(lapply(myFiles, function(x) { gsub("_Counts.tsv.gz", "", x) } ))
myLibWT <- myLib[grepl("_WT_", myLib)]

for (i in 1:length(myLibWT))
	{
	mylib <- myLibWT[[i]]
	print(paste("######################### ", mylib, sep = ""))

	print(paste("######################### ", "Reading count table for Mouse ", mylib, sep = ""))
	exp.data.m <- read.table(paste("./COUNTS_DIR_MOUSE/", mylib, "_Counts.tsv.gz", sep = ""), header = T, sep = "\t", row.names = 1)
	dim(exp.data.m)
	exp.data.m$genes <- row.names(exp.data.m)

	libname <- mylib 

	print(paste("######################### ", "Adding missing rows", sep = ""))
	mergedTemp <- list(REFGENES = refGenesSymbol, EXPDATA = exp.data.m)
	mergedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "genes") } , mergedTemp)
	row.names(mergedData) <- mergedData$genes
	mergedData$genes <- NULL
	print(dim(mergedData))

	print(paste("######################### ", "Fixing NAs", sep = ""))
	mergedData[is.na(mergedData)] <- 0

	counttable <- as(as.matrix(mergedData), "sparseMatrix") 

	mygenes <- as.data.frame(row.names(counttable))
	colnames(mygenes) <- "GeneID"
	mygenes$GeneName <- mygenes$GeneID
	mygenes$Category <- rep("Gene Expression", nrow(mygenes))
	mycellbarcodes <- colnames(counttable)

	print(paste("######################### ", "Saving output", sep = ""))	
	mydir <- paste(mylib, "CB", sep = "_")
	dir.create(mydir)
	write.table(mygenes, paste(mydir,"/features.tsv", sep = ""), col.names = F, sep = '\t', quote = F, row.names = F)
	write.table(mycellbarcodes, paste(mydir,"/barcodes.tsv", sep = ""), col.names = F, sep = '\t', quote = F, row.names = F)
	writeMM(counttable, paste(mydir, "/matrix.mtx", sep = ""))

	# write.table(mergedData, paste(libname, "_Counts_NovaSeq.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)
	save(counttable, file = paste(libname, "_Counts_NovaSeq.RData", sep = ""))
	}

#########
## END ##
#########
```
<br></br>



## KO SAMPLES
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
####################
## LOAD LIBRARIES ##
####################
rm(list = ls())
library(Seurat)
library(Matrix)
library(SeuratDisk)


####################
## REF GENES      ##
####################
refGenes <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", sep = "\t", header = T)
refGenesSymbol <- as.data.frame(unique(sort(refGenes$GeneSymbol)))
colnames(refGenesSymbol) <- "genes"
row.names(refGenesSymbol) <- refGenesSymbol$genes
print(dim(refGenesSymbol))
# 21967     1


####################
## LOOP OVER ALL  ##
####################
myFiles <- list.files(path = "./COUNTS_DIR_MOUSE/", pattern = "*.gz")
myLib <- unlist(lapply(myFiles, function(x) { gsub("_Counts.tsv.gz", "", x) } ))
myLibKO <- myLib[grepl("_KO_", myLib)]

for (i in 1:length(myLibKO))
	{
	mylib <- myLibKO[[i]]
	print(paste("######################### ", mylib, sep = ""))

	print(paste("######################### ", "Reading count table for Mouse ", mylib, sep = ""))
	exp.data.m <- read.table(paste("./COUNTS_DIR_MOUSE/", mylib, "_Counts.tsv.gz", sep = ""), header = T, sep = "\t", row.names = 1)
	dim(exp.data.m)
	exp.data.m$genes <- row.names(exp.data.m)

	libname <- mylib 

	print(paste("######################### ", "Adding missing rows", sep = ""))
	mergedTemp <- list(REFGENES = refGenesSymbol, EXPDATA = exp.data.m)
	mergedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "genes") } , mergedTemp)
	row.names(mergedData) <- mergedData$genes
	mergedData$genes <- NULL
	print(dim(mergedData))

	print(paste("######################### ", "Fixing NAs", sep = ""))
	mergedData[is.na(mergedData)] <- 0

	counttable <- as(as.matrix(mergedData), "sparseMatrix") 

	mygenes <- as.data.frame(row.names(counttable))
	colnames(mygenes) <- "GeneID"
	mygenes$GeneName <- mygenes$GeneID
	mygenes$Category <- rep("Gene Expression", nrow(mygenes))
	mycellbarcodes <- colnames(counttable)

	print(paste("######################### ", "Saving output", sep = ""))	
	mydir <- paste(mylib, "CB", sep = "_")
	dir.create(mydir)
	write.table(mygenes, paste(mydir,"/features.tsv", sep = ""), col.names = F, sep = '\t', quote = F, row.names = F)
	write.table(mycellbarcodes, paste(mydir,"/barcodes.tsv", sep = ""), col.names = F, sep = '\t', quote = F, row.names = F)
	writeMM(counttable, paste(mydir, "/matrix.mtx", sep = ""))

	# write.table(mergedData, paste(libname, "_Counts_NovaSeq.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)
	save(counttable, file = paste(libname, "_Counts_NovaSeq.RData", sep = ""))
	}

#########
## END ##
#########
```
<br></br>



## HU SAMPLES
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
####################
## LOAD LIBRARIES ##
####################
rm(list = ls())
library(Seurat)
library(Matrix)
library(SeuratDisk)


####################
## REF GENES      ##
####################
refGenes <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", sep = "\t", header = T)
refGenesSymbol <- as.data.frame(unique(sort(refGenes$GeneSymbol)))
colnames(refGenesSymbol) <- "genes"
row.names(refGenesSymbol) <- refGenesSymbol$genes
print(dim(refGenesSymbol))
# 21967     1


####################
## LOOP OVER ALL  ##
####################
myFiles <- list.files(path = "./COUNTS_DIR_MOUSE/", pattern = "*.gz")
myLib <- unlist(lapply(myFiles, function(x) { gsub("_Counts.tsv.gz", "", x) } ))
myLibHU <- myLib[grepl("_HU_", myLib)]

for (i in 1:length(myLibHU))
	{
	mylib <- myLibHU[[i]]
	print(paste("######################### ", mylib, sep = ""))

	print(paste("######################### ", "Reading count table for Mouse ", mylib, sep = ""))
	exp.data.m <- read.table(paste("./COUNTS_DIR_MOUSE/", mylib, "_Counts.tsv.gz", sep = ""), header = T, sep = "\t", row.names = 1)
	dim(exp.data.m)
	exp.data.m$genes <- row.names(exp.data.m)

	libname <- mylib 

	print(paste("######################### ", "Adding missing rows", sep = ""))
	mergedTemp <- list(REFGENES = refGenesSymbol, EXPDATA = exp.data.m)
	mergedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "genes") } , mergedTemp)
	row.names(mergedData) <- mergedData$genes
	mergedData$genes <- NULL
	print(dim(mergedData))

	print(paste("######################### ", "Fixing NAs", sep = ""))
	mergedData[is.na(mergedData)] <- 0

	counttable <- as(as.matrix(mergedData), "sparseMatrix") 

	mygenes <- as.data.frame(row.names(counttable))
	colnames(mygenes) <- "GeneID"
	mygenes$GeneName <- mygenes$GeneID
	mygenes$Category <- rep("Gene Expression", nrow(mygenes))
	mycellbarcodes <- colnames(counttable)

	print(paste("######################### ", "Saving output", sep = ""))	
	mydir <- paste(mylib, "CB", sep = "_")
	dir.create(mydir)
	write.table(mygenes, paste(mydir,"/features.tsv", sep = ""), col.names = F, sep = '\t', quote = F, row.names = F)
	write.table(mycellbarcodes, paste(mydir,"/barcodes.tsv", sep = ""), col.names = F, sep = '\t', quote = F, row.names = F)
	writeMM(counttable, paste(mydir, "/matrix.mtx", sep = ""))

	# write.table(mergedData, paste(libname, "_Counts_NovaSeq.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)
	save(counttable, file = paste(libname, "_Counts_NovaSeq.RData", sep = ""))
	}

#########
## END ##
#########
```

-----
