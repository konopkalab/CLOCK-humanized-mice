# BUILD CELLRANGER REFERENCE

## Chimpanzee
```{bash}
## START
echo "Loading Modules"
module load cellranger/8.0.1

echo "Setting Paths"
genomefa="/project/RESOURCES/GENOME_INDEX/panTro6/panTro6_selected.fa"
annotgtf="/project/RESOURCES/GENOME_INDEX/panTro6/hg38_liftoff_to_pantro6.gtf"

echo "Generating Genome Index for CellRanger"
cellranger mkref --genome CR_PANTRO6_LIFTOFF --fasta ${genomefa} --genes ${annotgtf}

echo "Done"
## END
```

## Human
```{bash}
## START
echo "Loading Modules"
module load cellranger/8.0.1

echo "Setting Paths"
genomefa="/project\RESOURCES/GENOME_INDEX/GRCh38p12/GRCh38.p12.genome.selected.fa"
annotgtf="/project/RESOURCES/GENOME_INDEX/Gencode_vH44/gencode.v44.basic.annotation.protein_coding.gtf"

echo "Generating Genome Index for CellRanger"
cellranger mkref --genome CR_GRCH38_GENCODEH44 --fasta ${genomefa} --genes ${annotgtf}

echo "Done"
## END
```

## Macaque
```{bash}
## START
echo "Loading Modules"
module load cellranger/8.0.1

echo "Setting Paths"
genomefa="/project/RESOURCES/GENOME_INDEX/rheMac10/rheMac10_selected.fa"
annotgtf="/project/RESOURCES/GENOME_INDEX/rheMac10/hg38_liftoff_to_rhemac10.gtf"

echo "Generating Genome Index for CellRanger"
cellranger mkref --genome CR_RHEMAC10_LIFTOFF --fasta ${genomefa} --genes ${annotgtf}

echo "Done"
## END
```

