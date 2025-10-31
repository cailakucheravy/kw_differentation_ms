#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --mem=20GB
#SBATCH --job-name=bam_coverage
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# load modules
module load StdEnv/2020 bwa/0.7.17 samtools/1.12

# go to directory
cd /scratch/cailak/killer_whales3/6_samtools_filter

# run the loop for each file in the directory
for i in *.bam
do
echo "Running bam_coverage: "$i        
./bam_coverage.sh $i
samtools depth  $i  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > $i.samdepth
done