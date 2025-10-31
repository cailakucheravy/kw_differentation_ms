#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --job-name=addRG_KW
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# Load modules
module load nixpkgs/16.09 gcc/8.3.0 picard/2.20.6 samtools/1.9

# Create reference genome dictionary (.dict file) - only needs to be done once.
#java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=/scratch/cailak/killer_whales/ref_genome/KW_GCA_937001465.1_mOrcOrc1.1.scafname.fasta O=/scratch/cailak/killer_whales/ref_genome/KW_reference.dict 

# Go to directory with sorted bam files
cd /scratch/cailak/killer_whales3/4_deDup


# Add read group information to each deDup bam
ls /scratch/cailak/killer_whales3/4_deDup/*deDup.bam | sed 's/.deDup.bam$//' | parallel --jobs $SLURM_NTASKS 'java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={}.deDup.bam O={}.deDupRG.bam RGID={} RGPL=illumina RGSM={} RGLB=lib1 RGPU=unit1'

# Index new bam files
ls /scratch/cailak/killer_whales3/4_deDup/*deDupRG.bam | sed 's/.deDupRG.bam$//' | parallel --jobs $SLURM_NTASKS 'samtools index {}.deDupRG.bam'
