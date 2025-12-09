# Welcome!
This is a pipeline that preprocesses eDNA and clusters it into OTUs. This is designed specifically for fungi samples.

## Quick User Manual

### Pipeline V1
* Created for processing **one .fastq file** into OTUs.
* Inputs neccesary are the full paths of the **.fastq file** and of **a preferably empty output directory**.

Usage: 
```
sbatch pipeline_v1.sh <full/path/file.fastq> <full/path/empty_output_dir>
```

### Pipeline V2
* Created for processing **all .fastq files in a specified directory** into OTUs (one OTU_output per file).
* Inputs neccesary are full paths of the **directory with .fastq samples** and of **a preferably empty output directory**.

>NOTE: All .fastq files in the specified directory will be run by the pipeline. Although the code only selects and runs on fastq formats, ideally the directory should not have any other formats or subdirectories.*


Usage: 

```
sbatch pipeline_v2.sh <full/path/directory_with_samples> <full/path/empty_output_dir>
```

## Output format
### Pipeline V1
```
/output_dir
├── intermediate
│   └── all_intermediate_files
├── multiqc_report.html
├── OTU_centroids.fasta
├── OTU_consensus.fasta
└── QC
    └── all_fastqc.html
```
+ **intermediate** directory stores intermediate sample files generated after every step;
+ **multiqc_report.html** documents the changes after every processing step
+ **OTU_centroids** stores the most represented sequence for each OTU
+ **OTU_consensus** stores the average sequence for each OTU (like a "majority vote" of all sequences). If a nucleotide is ambiguous, it's represented with N;
+ **QC** directory stores FastQC analyses after each step. However all this information is found in the multiqc report.

### Pipeline V2
```
/output_dir
├── multiqc
│   └── {sample_name}_multiqc_reports.html
├── OTUs
│   ├── {sample_name}_OTU_centroids.fasta
│   └── {sample_name}_OTU_consensus.fasta
└── retained_reads
    └── {sample_name}_retained_ITS_region_reads.fasta
```
+ **multiqc** directory one multiqc report per sample, showing the changes after every processing step;
+ **OTUs** directory has one centroids and one consensus file per sample;
+ **retained_reads** directory contains one fasta with ITS regions per sample. This is the file used to create OTUs, and can be used for example to compare how many reads were retained.

>NOTE: Because of the possible massive amount of intermediate files, pipeline_v2 stores all intermediate files in $TMPDIR, which is local scratch memory. This is a form of node-local storage made for temporary files that are deleted after the job finishes running. **It's possible that this will not work if running the job outside of RACKHAM**.
