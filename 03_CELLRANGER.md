
## Chimpanzee | Jorstad et al., 2023
```{bash}
## START
echo "Loading Modules"
module purge && module load shared slurm
module load cellranger/8.0.1

echo "Setting Paths"
txome="/project/RESOURCES/GENOME_INDEX/CELLRANGER_panTro6/CR_PANTRO6_LIFTOFF"
fqdir="/project/A_ALLEN_LAB/CHIMPANZEE"

## Samples
# NW_TX0017-12
# NW_TX0017-13
# NW_TX0020-6
# NW_TX0020-7
# NW_TX0020-10
# NW_TX0020-13
# NW_TX0025-10
# NW_TX0025-11
# NW_TX0025-12
# NW_TX0027-6
# NW_TX0027-7
# NW_TX0027-8
# NW_TX0027-9
# NW_TX0028-1
# NW_TX0053-8
# NW_TX0053-9
# NW_TX0055-7
# NW_TX0055-8
# NW_TX0055-9
# NW_TX0067-8
# NW_TX0067-9
# NW_TX0067-10
# NW_TX0068-8
# NW_TX0068-9
# NW_TX0070-1
# NW_TX0070-6

echo "Running CellRanger Count | NW_TX0017-12"
cellranger count --id NW_TX0017-12 --transcriptome ${txome} --fastqs ${fqdir} --sample NW_TX0017-12 --include-introns false --create-bam true

## Repeat for all samples listed above.
## END
```

## Human | Jorstad et al., 2023
```{bash}
## START
echo "Loading Modules"
module purge && module load shared slurm
module load cellranger/8.0.1

echo "Setting Paths"
txome="/project/RESOURCES/GENOME_INDEX/CELLRANGER_GRCh38p12_GENCODEvH44/CR_GRCH38_GENCODEH44"
fqdir="/project/A_ALLEN_LAB/HUMAN"

## Samples
# NW_TX0024-8
# NW_TX0024-9
# NW_TX0024-12
# NW_TX0024-13
# 10X359-4
# 10X359-5
# 10X216-7
# 10X216-8
# 10X240-1
# 10X240-2
# 10X240-3
# 10X240-4
# 10X240-5
# 10X240-6
# 10X240-7
# 10X240-8
# 10X241-1
# 10X241-2
# 10X241-3
# 10X241-4
# 10X241-5
# 10X241-6
# 10X241-7
# 10X241-8

echo "Running CellRanger Count | NW_TX0024-8"
cellranger count --id NW_TX0024-8 --transcriptome ${txome} --fastqs ${fqdir} --sample NW_TX0024-8 --include-introns false --create-bam true

## Repeat for all samples listed above.
## END
```

## Macaque | Jorstad et al., 2023
```{bash}
## START
echo "Loading Modules"
module purge && module load shared slurm
module load cellranger/8.0.1

echo "Setting Paths"
txome="/project/RESOURCES/GENOME_INDEX/CELLRANGER_rheMac10/CR_RHEMAC10_LIFTOFF"
fqdir="/project/A_ALLEN_LAB/MACAQUE"

## Samples
# NW_TX0023-7
# NW_TX0023-8
# NW_TX0023-9
# NW_TX0023-10
# NW_TX0023-11
# NW_TX0023-12
# NW_TX0023-13
# NW_TX0023-14
# NW_TX0024-2
# NW_TX0024-3
# NW_TX0024-4
# NW_TX0024-5
# NW_TX0032-4
# NW_TX0032-5
# NW_TX0032-6
# NW_TX0032-7
# NW_TX0032-8
# NW_TX0032-9
# NW_TX0033-16
# NW_TX0053-5
# NW_TX0053-6
# NW_TX0053-7
# NW_TX0073-10
# NW_TX0074-6
# NW_TX0074-7

echo "Running CellRanger Count | NW_TX0073-10"
cellranger count --id NW_TX0073-10 --transcriptome ${txome} --fastqs ${fqdir} --sample NW_TX0073-10 --include-introns false --create-bam true

## Repeat for all samples listed above.
## END
```

## Chimpanzee | Caglayan et al., 2023
```{bash}
## START
echo "Loading Modules"
module purge && module load shared slurm
module load cellranger/8.0.1

echo "Setting Paths"
txome="/project/RESOURCES/GENOME_INDEX/CELLRANGER_panTro6/CR_PANTRO6_LIFTOFF"
fqdir="/project/A_KONOPKA_LAB/01_FASTQ_RNA"

## Samples
# C1
# C2
# C3
# C4

echo "Running CellRanger Count | C1"
cellranger count --id C1 --transcriptome ${txome} --fastqs ${fqdir} --sample C1 --include-introns false --create-bam true

echo "Running CellRanger Count | C2"
cellranger count --id C2 --transcriptome ${txome} --fastqs ${fqdir} --sample C2 --include-introns false --create-bam true

echo "Running CellRanger Count | C3"
cellranger count --id C3 --transcriptome ${txome} --fastqs ${fqdir} --sample C3 --include-introns false --create-bam true

echo "Running CellRanger Count | C4"
cellranger count --id C4 --transcriptome ${txome} --fastqs ${fqdir} --sample C4 --include-introns false --create-bam true

echo "Done"
## END
```

## Human | Caglayan et al., 2023
```{bash}
## START
echo "Loading Modules"
module purge && module load shared slurm
module load cellranger/8.0.1

echo "Setting Paths"
txome="/project/RESOURCES/GENOME_INDEX/CELLRANGER_GRCh38p12_GENCODEvH44/CR_GRCH38_GENCODEH44"
fqdir="/project/A_KONOPKA_LAB/01_FASTQ_RNA"

## Samples
# H1
# H2
# H3
# H4

echo "Running CellRanger Count | H1"
cellranger count --id H1 --transcriptome ${txome} --fastqs ${fqdir} --sample H1 --include-introns false --create-bam true

echo "Running CellRanger Count | H2"
cellranger count --id H2 --transcriptome ${txome} --fastqs ${fqdir} --sample H2 --include-introns false --create-bam true

echo "Running CellRanger Count | H3"
cellranger count --id H3 --transcriptome ${txome} --fastqs ${fqdir} --sample H3 --include-introns false --create-bam true

echo "Running CellRanger Count | H4"
cellranger count --id H4 --transcriptome ${txome} --fastqs ${fqdir} --sample H4 --include-introns false --create-bam true

echo "Done"
## END
```

## Macaque | Caglayan et al., 2023
```{bash}
## START
echo "Loading Modules"
module purge && module load shared slurm
module load cellranger/8.0.1

echo "Setting Paths"
txome="/project/RESOURCES/GENOME_INDEX/CELLRANGER_rheMac10/CR_RHEMAC10_LIFTOFF"
fqdir="/project/A_KONOPKA_LAB/01_FASTQ_RNA"

## Samples
# M1
# M2
# M3
# M4

echo "Running CellRanger Count | M1"
cellranger count --id M1 --transcriptome ${txome} --fastqs ${fqdir} --sample M1 --include-introns false --create-bam true

echo "Running CellRanger Count | M2"
cellranger count --id M2 --transcriptome ${txome} --fastqs ${fqdir} --sample M2 --include-introns false --create-bam true

echo "Running CellRanger Count | M3"
cellranger count --id M3 --transcriptome ${txome} --fastqs ${fqdir} --sample M3 --include-introns false --create-bam true

echo "Running CellRanger Count | M4"
cellranger count --id M4 --transcriptome ${txome} --fastqs ${fqdir} --sample M4 --include-introns false --create-bam true

echo "Done"
## END
```

