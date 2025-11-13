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
INPUTFILE=$2

# Check if arguments are missing
if [ -z "$WORKDIR" ] || [ -z "$INPUTFILE" ]; then
    echo "Usage: sbatch run_job.sh <workdir> <inputfile>"
    exit 1
fi

echo "Workdir:  $WORKDIR"
echo "Input:    $INPUTFILE"

############# CONDA ############# 

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate b_thesis_env

# ============ File management ============

cd "$WORKDIR" || { echo "Cannot enter $WORKDIR"; exit 1; }
mkdir -p intermediate
mkdir -p QC

############# CODE ############## 

echo "================== FastQC 1 =======================" 

fastqc "$INPUTFILE" -o "$WORKDIR/QC" -t 4

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
  -o "$WORKDIR/intermediate/cutadapt_trimmed_reads.fastq.gz" \
  "$INPUTFILE"

echo "================== FastQC 2 =======================" 

fastqc "$WORKDIR/intermediate/cutadapt_trimmed_reads.fastq.gz" -o "$WORKDIR/QC" -t 4

echo "===================================================" 
echo "============= 2.  Quality filtering ===============" 
echo "================== with vsearch ==================="
echo "===================================================" 

vsearch \
 --fastq_filter "$WORKDIR/intermediate/cutadapt_trimmed_reads.fastq.gz" \
 --fastq_maxee 12 \
 --fastq_minlen 1000 \
 --fastq_maxlen 2400 \
 --fastq_maxns 0 \
 --fastqout "$WORKDIR/intermediate/q_filtered_reads.fastq"

echo "================== FastQC 3 =======================" 

fastqc "$WORKDIR/intermediate/q_filtered_reads.fastq" -o "$WORKDIR/QC" -t 4

echo "===================================================" 
echo "====== 3. Deduplication. Abundance counting. ======" 
echo "================== with vsearch ==================="
echo "===================================================" 

vsearch \
 --fastx_uniques "$WORKDIR/intermediate/q_filtered_reads.fastq" \
 --sizeout \
 --fastqout "$WORKDIR/intermediate/dereplicated_reads.fastq"

vsearch --sortbysize "$WORKDIR/intermediate/dereplicated_reads.fastq" --output "$WORKDIR/intermediate/sorted_dereplicated_reads.fasta"


echo "================== FastQC 4 =======================" 


fastqc "$WORKDIR/intermediate/dereplicated_reads.fastq" -o "$WORKDIR/QC" -t 4

echo "===================================================" 
echo "============= 4. Chimera filtering ================" 
echo "================== with vsearch ==================="
echo "===================================================" 

vsearch \
 --uchime_denovo "$WORKDIR/intermediate/sorted_dereplicated_reads.fasta" \
 --nonchimeras "$WORKDIR/intermediate/nonchimeras.fasta" \
 --abskew 2 \
 --mindiv 0.28 \
 --minh 0.28 

echo "================== MultiQC =======================" 

multiqc "$WORKDIR/QC" -o "$WORKDIR"

echo "===================================================" 
echo "=========== 5. ITS region extraction ==============" 
echo "=================== with ITSx ====================="
echo "===================================================" 

ITSx \
 -i "$WORKDIR/intermediate/nonchimeras.fasta" \
 -o "$WORKDIR/intermediate/ITSx_out" \
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
 --cluster_fast "$WORKDIR/intermediate/ITSx_out.ITS2.full_and_partial.fasta" \
 --biomout OTU_table.biom \
 --id 0.985 \
 --iddef 1 \
 --strand plus \
 --centroids "$WORKDIR/OTU_centroids.fasta" \
 --consout "$WORKDIR/OTU_consensus.fasta" \
 --sizeout \
 --threads 4 \
 --qmask dust 


