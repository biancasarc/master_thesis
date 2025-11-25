#!/bin/bash
#SBATCH -A uppmax2025-2-64
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J annot
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

set -euo pipefail

############# CONDA ############# 

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate b_thesis_env

vsearch --sintax /proj/toband/bianca/workdir/OTU_centroids.fasta \
 --db /proj/toband/bianca/UNITE_db/prepared_UNITE_db.fasta \
 --sintax_cutoff 0.8 \
 --tabbedout /proj/toband/bianca/workdir/taxon_table_2.tsv \
 --threads 4
