#!/bin/bash

#SBATCH --time=14-00:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=3
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --job-name=map_fastq
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# Load modules
module load StdEnv/2020 bwa/0.7.17 samtools/1.12

# Go to directory with final fastqs that have been trimmed and merged
cd /scratch/cailak/killer_whales2/trimmed/merged

# Map paired.fastq's to reference genome 
ls *_R1_paired.fastq.gz | sed 's/_R1_paired.fastq.gz$//' | parallel --jobs $SLURM_NTASKS 'bwa mem /scratch/cailak/ref_genome/KW_GCA_937001465.1_mOrcOrc1.1.scafname.fasta {}_R1_paired.fastq.gz {}_R2_paired.fastq.gz | samtools sort -T {}_tmp -o {}.sorted.bam'
