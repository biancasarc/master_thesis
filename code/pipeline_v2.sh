#!/bin/bash
#SBATCH -A uppmax2025-2-64
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH -J OTU_pipeline
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

set -euo pipefail

######## HANDLE ARGUMENTS ######## 

WORKDIR=$1
INPUTDIR=$2

# Check if arguments are missing
if [ -z "$WORKDIR" ] || [ -z "$INPUTDIR" ]; then
    echo "Usage: sbatch pipeline.sh <workdir> <inputdir>"
    exit 1
fi

echo "Workdir:  $WORKDIR"
echo "Input:    $INPUTDIR"

############# CONDA ############# 

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate b_thesis_env


# ============ File management ============

cd "$WORKDIR" || { echo "Cannot enter $WORKDIR"; exit 1; }
mkdir multiqc       #here will be one multiqc/sample, to show the progress of the processing

shopt -s nullglob   #stops pipeline from running if no files are found

for sample in "$INPUTDIR"/*.fastq.gz; do

    echo "Processing $sample"

    name="$(basename "${sample%.fastq.gz}")"      #extracts the name of the sample without fastq.gz

    # ============ Internal directory management ============

    mkdir "intermediate_${name}"
    mkdir "QC_${name}"

    ############# CODE ############## 

    echo "================== FastQC 1 =======================" 

    fastqc "$sample" -o "QC_${name}" -t 4

    echo "===================================================" 
    echo "============= 1. Trimming adapters ================" 
    echo "================ with cutadapt ===================="
    echo "===================================================" 

    cutadapt --cores 4 \
    -a TCCGTAGGTGAACCTGC...CGAAGTTTCCCTCAGGA \
    --revcomp \
    --error-rate 0.1 \
    --overlap 10 \
    --no-indels \
    --action=trim \
    --discard-untrimmed \
    -o "intermediate_${name}/01_trimmed_reads.fastq.gz" \
    "$sample"

    echo "================== FastQC 2 =======================" 

    fastqc "intermediate_${name}/01_trimmed_reads.fastq.gz" -o "QC_${name}" -t 4

    echo "===================================================" 
    echo "============= 2.  Quality filtering ===============" 
    echo "================== with vsearch ==================="
    echo "===================================================" 

    vsearch \
    --fastq_filter "intermediate_${name}/01_trimmed_reads.fastq.gz" \
    --fastq_maxee 12 \
    --fastq_minlen 1000 \
    --fastq_maxlen 2400 \
    --fastq_maxns 0 \
    --fastqout "intermediate_${name}/02_q_filtered_reads.fastq"

    echo "================== FastQC 3 =======================" 

    fastqc "intermediate_${name}/02_q_filtered_reads.fastq" -o "QC_${name}" -t 4

    echo "===================================================" 
    echo "====== 3. Deduplication. Abundance counting. ======" 
    echo "================== with vsearch ==================="
    echo "===================================================" 

    vsearch \
    --fastx_uniques "intermediate_${name}/02_q_filtered_reads.fastq" \
    --sizeout \
    --fastqout "intermediate_${name}/03_dereplicated_reads.fastq"

    vsearch --sortbysize "intermediate_${name}/03_dereplicated_reads.fastq" --output "intermediate_${name}/04_sorted_dereplicated_reads.fasta"


    echo "================== FastQC 4 =======================" 


    fastqc "intermediate_${name}/03_dereplicated_reads.fastq" -o "QC_${name}" -t 4

    echo "================== MultiQC =======================" 

    multiqc "QC_${name}" -o .
    mv multiqc_report.html "multiqc/${name}_multiqc_report.html"   #moves and renames multiqc report to specified directory
    rm -r multiqc_data                  #removes intermediate file that multiqc creates automatically

    echo "===================================================" 
    echo "============= 4. Chimera filtering ================" 
    echo "================== with vsearch ==================="
    echo "===================================================" 

    vsearch \
    --uchime_denovo "intermediate_${name}/04_sorted_dereplicated_reads.fasta" \
    --nonchimeras "intermediate_${name}/nonchimeras.fasta" \
    --abskew 2 \
    --mindiv 0.28 \
    --minh 0.28 


    echo "===================================================" 
    echo "=========== 5. ITS region extraction ==============" 
    echo "=================== with ITSx ====================="
    echo "===================================================" 

    ITSx \
    -i "intermediate_${name}/nonchimeras.fasta" \
    -o "intermediate_${name}/ITSx_out" \
    -t F \
    --cpu 4 \
    -E 1e-2 \
    --partial 50 \
    --complement F \
    --only_full T \
    --truncate T \
    --preserve T \
    --complement T \
    --save_regions ITS2



    echo "===================================================" 
    echo "============ 6. Clustering into OTUs ==============" 
    echo "================== with vsearch ==================="
    echo "===================================================" 

    vsearch \
    --cluster_fast "intermediate_${name}/ITSx_out.ITS2.full_and_partial.fasta" \
    --id 0.985 \
    --iddef 1 \
    --strand plus \
    --centroids "${name}_OTU_centroids.fasta" \
    --consout "${name}_OTU_consensus.fasta" \
    --sizeout \
    --threads 4 \
    --qmask dust 


    echo "================== Removing intermediate directories =======================" 

    rm -r "intermediate_${name}"
    rm -r "QC_${name}"
done
