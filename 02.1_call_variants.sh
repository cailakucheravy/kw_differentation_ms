#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --job-name=platypus
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 gcc/7.3.0 platypus/0.8.1

# Set up for killer whale

# Go to directory with final bam files (in hindsight can just make path to each bam file in the .txt list file instead of going to folder).
cd /scratch/edegreef/killerwhale/bams_for_SNPs

# Run platypus to call variants from bams & reference genome 
# Make sure all bam files are indexed beforehand. Also couldn't add Platypus.py to path for some reason so I inputed the path here in script
python /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/platypus/0.8.1/bin/Platypus.py callVariants \
--bamFiles=killerwhale_bams_list.txt \
--refFile=/scratch/edegreef/killerwhale/ref_genome/KW_GCA_937001465.1_mOrcOrc1.1.scafname.fasta \
--output=/scratch/edegreef/killerwhale/killerwhale_allvariantcalls.vcf \
--nCPU=32 --minReads=4