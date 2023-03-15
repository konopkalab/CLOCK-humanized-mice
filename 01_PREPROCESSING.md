# H1 01 PRE-PROCESSING

## H2 Extracting FASTQ files from BCL
## bcl2fastq(2.19.1)
## cellranger(3.0.2)
````
<path_to_cellranger>/cellranger-3.0.2/cellranger mkfastq --csv=samplesheet.csv --ignore-dual-index --run=<path_to_uncompressed_bcl_folder>
````

## samplesheet.csv
## only an example for one sample is shown
| Lane | Sample | Index |
| 1 | YLP0701 | SI-GA-A1 |

## H2 FASTQC
## fastqc(0.11.5)
```
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
for wdir in YLP0701
    do
        echo "=============> Processing" ${FQDIR}"/"${wdir}
        cd ${FQDIR}"/"${wdir}
        ls *R*.gz | xargs -I % -n 1 -P 24 sh -c 'echo %; fastqc -t 24 -q %'
    done
```

## H2 CELLBARCODE WHITELIST
## python(3.6.1-2-anaconda)
## UMI-tools(0.5.4)
```
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
cd ${FQDIR}
ls *R1*.gz | sed "s/_L002_R1_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; <path_to_umi_tools_executable>/umi_tools whitelist --stdin=%_L002_R1_001.fastq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --expect-cells=10000 --plot-prefix=%_Expect_Whitelist --log=%_Whitelist_log.out --error=%_Whitelist_log.err --stdout=%_Whitelist.txt'
```

## H2 EXTRACT READS FOR CELLBARCODE WHITELIST
## python(3.6.1-2-anaconda)
## UMI-tools(0.5.4)
```
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
cd ${FQDIR}
ls *R1*.gz | sed "s/_L002_R1_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; <path_to_umi_tools_executable>/umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --stdin=%_L002_R1_001.fastq.gz --stdout=%_R1_extracted.fastq.gz --read2-in=%_L002_R2_001.fastq.gz --read2-out=%_R2_extracted.fastq.gz --filter-cell-barcode --whitelist=%_Whitelist.txt --log=%_Extract_log.out --error=%_Extract_log.err'
```


## H2 ALIGN READS TO REFERENCE GENOME & ANNOTATION
## samtools(1.6)
## star(2.5.2b)
```
EXTDIR=<path_to_directory_containing_extracted_fastq_files_corresponding_to_whitelisted_cellbarcodes>
cd ${EXTDIR}

## alignment to get unfiltered reads
for FQFILE in `ls *R2_extracted*.gz`
 do
  # echo ${FQFILE}
  prefx=`echo ${FQFILE} | sed "s/_R2_extracted.fastq.gz//g"`
  echo "Processing" ${prefx}

  STAR --runThreadN 48 \
       --genomeDir /work/Neuroinformatics_Core/akulk1/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_STAR \
       --readFilesIn ${FQFILE} \
       --readFilesCommand zcat \
       --sjdbGTFfile /work/Neuroinformatics_Core/akulk1/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_STAR/gencode.vM17.protein_coding.gtf \
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





## python(3.6.1-2-anaconda)
## samtools(1.6)
## star(2.5.2b)
## UMI-tools(0.5.4)
