#!/bin/bash
#SBATCH -A uppmax2026-1-17
#SBATCH -p pelle
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --cpus-per-task=1
#SBATCH -t 72:00:00
#SBATCH -J UNITE_renaming
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate /proj/toband/conda_envs/biancas_pipeline_env/env/b_thesis_env

python renaming.py --unite_db UNITE_public_19.02.2025.fasta --curated_db sh_general_release_dynamic_19.02.2025.fasta


