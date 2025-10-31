#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=10GB
#SBATCH --job-name=index_bams
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# Load modules
module load StdEnv/2020 bwa/0.7.17 samtools/1.12

# directory
cd /scratch/cailak/killer_whales3/3_bams_sorted

# Index all bam files in folder
ls *.sorted.bam | parallel --jobs $SLURM_NTASKS 'samtools index {}'
