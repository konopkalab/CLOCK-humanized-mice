# 01 | PRE-PROCESSING

## Extracting FASTQ files from BCL
- bcl2fastq(2.19.1)
- cellranger(3.0.2)
````
<path_to_cellranger>/cellranger-3.0.2/cellranger mkfastq --csv=samplesheet.csv --ignore-dual-index --run=<path_to_uncompressed_bcl_folder>
````

'samplesheet.csv' (only an example for one sample is shown).

| Lane | Sample | Index |
| ----- | ----- | ----- |
| 1 | YLP0701 | SI-GA-A1 |



## FASTQC
- fastqc(0.11.5)
```
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
for wdir in YLP0701
    do
        echo "=============> Processing" ${FQDIR}"/"${wdir}
        cd ${FQDIR}"/"${wdir}
        ls *R*.gz | xargs -I % -n 1 -P 24 sh -c 'echo %; fastqc -t 24 -q %'
    done
```


## CELLBARCODE WHITELIST
- python(3.6.1-2-anaconda)
- UMI-tools(0.5.4)
```
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
cd ${FQDIR}
ls *R1*.gz | sed "s/_L002_R1_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; <path_to_umi_tools_executable>/umi_tools whitelist --stdin=%_L002_R1_001.fastq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --expect-cells=10000 --plot-prefix=%_Expect_Whitelist --log=%_Whitelist_log.out --error=%_Whitelist_log.err --stdout=%_Whitelist.txt'
```


## EXTRACT READS FOR CELLBARCODE WHITELIST
- python(3.6.1-2-anaconda)
- UMI-tools(0.5.4)
```
FQDIR=<path_to_cellranger_mkfastq_output_directory_containing_fastq_files>
cd ${FQDIR}
ls *R1*.gz | sed "s/_L002_R1_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; <path_to_umi_tools_executable>/umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --stdin=%_L002_R1_001.fastq.gz --stdout=%_R1_extracted.fastq.gz --read2-in=%_L002_R2_001.fastq.gz --read2-out=%_R2_extracted.fastq.gz --filter-cell-barcode --whitelist=%_Whitelist.txt --log=%_Extract_log.out --error=%_Extract_log.err'
```

