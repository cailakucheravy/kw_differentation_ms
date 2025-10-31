#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --mem=20G
#SBATCH --cpus-per-task=6
#SBATCH --job-name=trim_fastq
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# Load modules
module load nixpkgs/16.09 trimmomatic/0.36

# Use TruSeq3-PE-2.fa to account for read-through, obtained at https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa

# Added to TruSeq3-PE-2.fa to account for Nextera transposase sequence:
	# >transposase1
	# GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG
	# >transposase1_rc
	# CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC
	# >transposase2
	# GCCTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG
	# >transposase2_rc
	# CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC	
	
	
# Go to directory with all the raw reads
cd /scratch/cailak/killer_whales3/0_fastqs_raw

# Make folder for output
mkdir trimmed

# Run trimmomatic (make sure path to TruSeq3-PE-2.fa is accurate)
ls ./*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' | parallel --jobs $SLURM_CPUS_PER_TASK 'java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE {}_R1_001.fastq.gz {}_R2_001.fastq.gz trimmed/{}_R1_paired.fastq.gz trimmed/{}_R1_unpaired.fastq.gz trimmed/{}_R2_paired.fastq.gz trimmed/{}_R2_unpaired.fastq.gz ILLUMINACLIP:/scratch/cailak/killer_whales2/fastqs_raw/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36'