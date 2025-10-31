#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kucherac@myumanitoba.ca
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --job-name=deduplicate
#SBATCH --output=%x-%j.out

# Script from E. de Greef.

# Load modules
module load nixpkgs/16.09 picard/2.20.6

# Go to directory with sorted bam files
cd /scratch/cailak/killer_whales3/3_bams_sorted

# Remove duplicate reads
# if running out of memory may need to do 1 job at a time instead of parallel (try parallel first).
ls *.sorted.bam | sed 's/.sorted.bam$//' | parallel --jobs $SLURM_NTASKS 'java -Xmx50000m -jar $EBROOTPICARD/picard.jar MarkDuplicates I={}.sorted.bam O={}.deDup.bam M={}_deDupMetrics.txt REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
