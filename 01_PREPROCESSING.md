# 01 - PRE-PROCESSING


## Extracting FASTQ files from BCL
- bcl2fastq (2.19.1)
- cellranger (3.0.2)
````{sh}
<path_to_cellranger>/cellranger-3.0.2/cellranger mkfastq --csv=samplesheet.csv --ignore-dual-index --run=<path_to_uncompressed_bcl_folder>
````

'samplesheet.csv' (only an example for one sample is shown).

| Lane | Sample | Index |
| ----- | ----- | ----- |
| 1 | YLP0701 | SI-GA-A1 |

<br></br>


## FASTQC
- fastqc (0.11.5)
```{sh}
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
cd ${FQDIR}
ls *R*.gz | xargs -I % -n 1 -P 24 sh -c 'echo %; fastqc -t 48 -q %'
```
<br></br>


## CELLBARCODE WHITELIST
- python (3.6.1-2-anaconda)
- UMI-tools (0.5.4)
```{sh}
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
cd ${FQDIR}
ls *R1*.gz | sed "s/_R1_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; <path_to_umi_tools_executable>/umi_tools whitelist --stdin=%_R1_001.fastq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --set-cell-number=50000 --plot-prefix=%_Expect_Whitelist --log=%_Whitelist_log.out --error=%_Whitelist_log.err --stdout=%_Whitelist.txt'
```
<br></br>


## EXTRACT READS FOR CELLBARCODE WHITELIST
- python (3.6.1-2-anaconda)
- UMI-tools (0.5.4)
```{sh}
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
cd ${FQDIR}
ls *R1*.gz | sed "s/_R1_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; <path_to_umi_tools_executable>/umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --stdin=%_R1_001.fastq.gz --stdout=%_R1_extracted.fastq.gz --read2-in=%_R2_001.fastq.gz --read2-out=%_R2_extracted.fastq.gz --filter-cell-barcode --whitelist=%_Whitelist.txt --log=%_Extract_log.out --error=%_Extract_log.err'
```
<br></br>


## GENERATE INDEX FOR REFERENCE MOUSE GENOME & MOUSE ANNOTATION
- star (2.5.2b)
- MM10 GRCm38p6 reference genome FA
- GENCODE vM17 reference annotation GTF
```{sh}
STAR --runMode genomeGenerate \
     --runThreadN 48 \
     --genomeDir <path_to_ref_mouse_index_output> \
     --genomeFastaFiles <path_to_reference_mouse_genome_fasta> \
     --sjdbGTFfile <path_to_ref_annotation_mouse_gtf> \
     --sjdbOverhang 100
```
<br></br>


## GENERATE INDEX FOR REFERENCE HUMAN GENOME & HUMAN CLOCK ANNOTATION
- star (2.5.2b)
- HG38 GRCh38p12 reference genome FA
- GENCODE vH28 reference annotation GTF for CLOCK
```{sh}
STAR --runMode genomeGenerate \
     --runThreadN 48 \
     --genomeDir <path_to_ref_human_index_output> \
     --genomeFastaFiles <path_to_reference_human_genome_fasta> \
     --sjdbGTFfile <path_to_ref_annotation_human_gtf_for_clock_gene> \
     --sjdbOverhang 100
```
<br></br>


## ALIGN READS TO REFERENCE MOUSE INDEX
- samtools (1.6)
- star (2.5.2b)
- Reference mouse index for STAR
```{sh}
EXTDIR=<path_to_directory_containing_extracted_fastq_files_corresponding_to_whitelisted_cellbarcodes>
cd ${EXTDIR}

for FQFILE in `ls *R2_extracted*.gz`
 do
  prefx=`echo ${FQFILE} | sed "s/_R2_extracted.fastq.gz//g"`
  echo "Processing" ${prefx}

  STAR --runThreadN 48 \
       --genomeDir <path_to_ref_mouse_index_output> \
       --readFilesIn ${FQFILE} \
       --readFilesCommand zcat \
       --sjdbGTFfile <path_to_ref_annotation_mouse_gtf> \
       --outFilterType BySJout  \
       --outFilterMismatchNoverReadLmax 0.04 \
       --outFilterMultimapNmax 10 \
       --alignSJoverhangMin 10 \
       --alignSJDBoverhangMin 1 \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterMismatchNmax 5 \
       --twopassMode Basic \
       --outFileNamePrefix ${prefx}_STAR_
 done
```
<br></br>


## ALIGN READS TO REFERENCE HUMAN INDEX
- samtools (1.6)
- star (2.5.2b)
- Reference human index for STAR
```{sh}
EXTDIR=<path_to_directory_containing_extracted_fastq_files_corresponding_to_whitelisted_cellbarcodes>
cd ${EXTDIR}

for FQFILE in `ls *R2_extracted*.gz`
 do
  prefx=`echo ${FQFILE} | sed "s/_R2_extracted.fastq.gz//g"`
  echo "Processing" ${prefx}

  STAR --runThreadN 48 \
       --genomeDir <path_to_ref_human_index_output> \
       --readFilesIn ${FQFILE} \
       --readFilesCommand zcat \
       --sjdbGTFfile <path_to_ref_annotation_human_gtf_for_clock_gene> \
       --outFilterType BySJout  \
       --outFilterMismatchNoverReadLmax 0.04 \
       --outFilterMultimapNmax 10 \
       --alignSJoverhangMin 10 \
       --alignSJDBoverhangMin 1 \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterMismatchNmax 5 \
       --twopassMode Basic \
       --outFileNamePrefix ${prefx}_STAR_
 done
```
<br></br>



## ASSIGN UNIQUELY MAPPED READS TO GENES FOR MOUSE
- featureCounts (Subread-1.6.2)
- GENCODE vM17 reference annotation GTF
```{sh}
STARDIR=<path_to_star_output>
cd ${STARDIR}
## STRANDED PRIMARY
for BAMFILE in `ls *_STAR_Aligned.sortedByCoord.out.bam`
 do
  echo ${BAMFILE}
  prefx=`echo ${BAMFILE} | sed "s/_STAR_Aligned.sortedByCoord.out.bam//g"`
  echo ${prefx}
  
  <path_to_feartureCounts_executable>/featureCounts \
    --primary \
    -R BAM \
    -T 48 \
    -s 1 \
    -t gene \
    -g gene_name \
    -a <path_to_ref_annotation_mouse_gtf> \
    -o ${prefx}_Primary_Gene_Assigned \
    ${BAMFILE}
 done
```
<br></br>



## ASSIGN UNIQUELY MAPPED READS TO GENES FOR HUMAN
- featureCounts (Subread-1.6.2)
- GENCODE vM17 reference annotation GTF
```{sh}
STARDIR=<path_to_star_output>
cd ${STARDIR}
## STRANDED PRIMARY
for BAMFILE in `ls *_STAR_Aligned.sortedByCoord.out.bam`
 do
  echo ${BAMFILE}
  prefx=`echo ${BAMFILE} | sed "s/_STAR_Aligned.sortedByCoord.out.bam//g"`
  echo ${prefx}
  
  <path_to_feartureCounts_executable>/featureCounts \
    --primary \
    -R BAM \
    -T 48 \
    -s 1 \
    -t gene \
    -g gene_name \
    -a <path_to_ref_annotation_human_gtf_for_clock_gene> \
    -o ${prefx}_Primary_Gene_Assigned \
    ${BAMFILE}
 done
```
<br></br>



## SORTING AND INDEXING BAM FOR MOUSE AND HUMAN
- samtools (1.6)
```{sh}
FCDIR=<path_to_featureCount_output_directory>
cd ${FCDIR}

## sorting bam
ls *featureCounts.bam | sed "s/.sortedByCoord.out.bam.featureCounts.bam//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; samtools sort %.sortedByCoord.out.bam.featureCounts.bam -o %_Assigned_Sorted.bam'

## indexing bam
ls *Assigned_Sorted.bam | xargs -I % -n 1 -P 48 sh -c 'echo %; samtools index %'
```
<br></br>



## GENERATE RAW COUNTS TABLE FOR MOUSE AND HUMAN
- python (3.6.1-2-anaconda)
- UMI-tools (0.5.4)
```{sh}
SBDIR=<path_to_directory_containing_sorted_indexed_bam>
cd ${SBDIR}

ls *_STAR_Aligned_Assigned_Sorted.bam | sed "s/_STAR_Aligned_Assigned_Sorted.bam//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; <path_to_umi_tools_executable>/umi_tools count --per-gene --gene-tag=XT --per-cell --stdin=%_STAR_Aligned_Assigned_Sorted.bam --stdout=%_Counts.tsv.gz --log=%_Counts.log --error=%_Counts.err --wide-format-cell-counts'
```
<br></br>
<br></br>

