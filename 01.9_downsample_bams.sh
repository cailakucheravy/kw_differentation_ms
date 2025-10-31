#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=downsample_batch1
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# Downsampling some of the bam files to a lower coverage to avoid SNP calling biases, aiming for 19x for the killer whale samples.
# 'keep' value proportion of reads to keep, and is based on target coverage divided by original coverage. For example if coverage=25x and want 19x, then keep=0.76
# Note that we are going by *modal* coverage, not mean.
# Second part of script is looking at new modal coverage, which will need bam_coverage.sh in directory

# List of sample name (including suffix from sequence lane)
# e.g., bam_ID1=OO_02_1_S2 for the file OO_02_1_S2.deDupRG.pp.bam. Here I put the current coverage as a comment next to sample ID
bam_ID1=KW_2021_PG_06_S18          #21x
bam_ID2=KW_2021_PG_09_S16          #22x
bam_ID3=KW_2021_PG_11_S15          #22x

# Keep1 corresponds to bam_ID1, keep2 corresponds to bam_ID2, etc...
# proportion of reads to keep (target coverage=19)
keep1=0.9047  #19/21
keep2=0.8636  #19/22
keep3=0.8636  #19/22

# Go to directory with .deDupRG.pp.bam files
cd /scratch/cailak/killer_whales3/6_samtools_filter

# Make a directory for downsampled outputs
mkdir downsampled

# Load gatk
module load nixpkgs/16.09 gatk/4.1.2.0

# Run DownsampleSam command
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID1.deDupRG.pp.bam -O downsampled/$bam_ID1.deDupRG.pp.downsampled.bam -P $keep1
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID2.deDupRG.pp.bam -O downsampled/$bam_ID2.deDupRG.pp.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID3.deDupRG.pp.bam -O downsampled/$bam_ID3.deDupRG.pp.downsampled.bam -P $keep3

### Check new modal coverage
# Load samtools
module load StdEnv/2020 samtools/1.11

./bam_coverage.sh downsampled/$bam_ID1.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID2.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID3.deDupRG.pp.downsampled.bam

# Index all new bams at the end
cd downsampled
for i in *.bam
do
echo "Indexing: "$i        
samtools index $i
done
