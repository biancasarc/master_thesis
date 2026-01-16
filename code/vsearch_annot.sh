#!/bin/bash
#SBATCH -A uppmax2026-1-17
#SBATCH -p pelle
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 1:00:00
#SBATCH -J annotation
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

set -euo pipefail

############# CONDA #############

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate b_thesis_env

vsearch --sintax /proj/toband/bianca/workdir/OTU_centroids.fasta \
 --db /proj/toband/bianca/uncorrected_UNITE_db_samples/prepared_UNITE_db.fasta \
 --sintax_cutoff 0.8 \
 --tabbedout /proj/toband/bianca/workdir/taxon_table.tsv \
 --threads 4

# --sintax_cutoff option is used to set a minimum level of bootstrap support for the taxonomic ranks to be reported
