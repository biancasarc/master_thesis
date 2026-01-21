#!/bin/bash
#SBATCH -A uppmax2026-1-17
#SBATCH -p pelle
#SBATCH --ntasks=1
#SBATCH -c 4q
#SBATCH --cpus-per-task=4
#SBATCH -t 4:00:00
#SBATCH -J OTU_pipeline
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

set -euo pipefail

######## HANDLE ARGUMENTS ######## 

INPUTFILE=$1
WORKDIR=$2

# Check if arguments are missing
if [ -z "$INPUTFILE" ] || [ -z "$WORKDIR" ]; then
    echo "Usage: sbatch pipeline.sh <inputfile> <workdir>"
    exit 1
fi

echo "Input:    $INPUTFILE"
echo "Workdir:  $WORKDIR"

############# CONDA ############# 

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate /proj/toband/conda_envs/biancas_pipeline_env/env/b_thesis_env

# ============ File management ============

cd "$WORKDIR" || { echo "Cannot enter $WORKDIR"; exit 1; }
mkdir intermediate
mkdir QC

############# CODE ############## 

echo "================== FastQC 1 =======================" 

fastqc "$INPUTFILE" -o QC/ -t 4

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
  -o intermediate/01_trimmed_reads.fastq.gz \
  "$INPUTFILE"

echo "================== FastQC 2 =======================" 

fastqc intermediate/01_trimmed_reads.fastq.gz -o QC/ -t 4

echo "===================================================" 
echo "============= 2.  Quality filtering ===============" 
echo "================== with vsearch ==================="
echo "===================================================" 

vsearch \
 --fastq_filter intermediate/01_trimmed_reads.fastq.gz \
 --fastq_qmax 93 \
 --fastq_maxee 12 \
 --fastq_minlen 1000 \
 --fastq_maxlen 2400 \
 --fastq_maxns 0 \
 --fastqout intermediate/02_q_filtered_reads.fastq

echo "================== FastQC 3 =======================" 

fastqc intermediate/02_q_filtered_reads.fastq -o QC/ -t 4

echo "===================================================" 
echo "====== 3. Deduplication. Abundance counting. ======" 
echo "================== with vsearch ==================="
echo "===================================================" 

vsearch \
 --fastx_uniques intermediate/02_q_filtered_reads.fastq \
 --sizeout \
 --fastqout intermediate/03_dereplicated_reads.fastq

vsearch --sortbysize intermediate/03_dereplicated_reads.fastq --output intermediate/04_sorted_dereplicated_reads.fasta


echo "================== FastQC 4 =======================" 


fastqc intermediate/03_dereplicated_reads.fastq -o QC/ -t 4

echo "===================================================" 
echo "============= 4. Chimera filtering ================" 
echo "================== with vsearch ==================="
echo "===================================================" 

vsearch \
 --uchime_denovo intermediate/04_sorted_dereplicated_reads.fasta \
 --nonchimeras intermediate/nonchimeras.fasta \
 --abskew 2 \
 --mindiv 0.28 \
 --minh 0.28 

echo "================== MultiQC =======================" 

multiqc QC/ -o .
rm -r multiqc_data   #remove directory created automatically by MultiQC
echo "===================================================" 
echo "=========== 5. ITS region extraction ==============" 
echo "=================== with ITSx ====================="
echo "===================================================" 

ITSx \
 -i intermediate/nonchimeras.fasta \
 -o intermediate/ITSx_out \
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
 --cluster_fast intermediate/ITSx_out.ITS2.full_and_partial.fasta \
 --id 0.985 \
 --iddef 1 \
 --strand plus \
 --centroids OTU_centroids.fasta \
 --consout OTU_consensus.fasta \
 --sizeout \
 --threads 4 \
 --qmask dust 

