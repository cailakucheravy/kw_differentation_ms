#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --mem=50G
#SBATCH --job-name=merge_fastq
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# Go to directory with the trimmed reads
cd /scratch/cailak/killer_whales3/trimmed

# Make directory for merged file outputs
mkdir merged

# Repeat following code for each sample:

#1
cat KW_2021_PG_03_S17_*_R1_paired.fastq.gz > merged/KW_2021_PG_03_S17_R1_paired.fastq.gz
cat KW_2021_PG_03_S17_*_R2_paired.fastq.gz > merged/KW_2021_PG_03_S17_R2_paired.fastq.gz


