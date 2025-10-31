#!/bin/bash

# SNP filtering for population structure analysis. 
# Need to check for kinship, remove duplicates and close kin prior to filtering snps.

# Includes:
# Step 1: Indel and PASS filter
# Step 2: QUAL filter, MQ filter, QD filter (after this step, renamed to "filter1")
# Step 3: Missingness filter, biallelic filter
# Step 4: Removing small scaffolds (<100kb)
# Step 5: Removing sex-linked snps
# Step 6: Additional HWE, MAF filter - separate populations for HWE filter.
# Step 7: LD pruning

# edit_sampleID.sh should be done to vcf beforehand if needed.
# Need filter.hets.R file in directory for the HWE filter to work.
# Need list of scaffolds (for removing small scaffold step) and list of X and Y scaffolds. 

# Set up pathways and prefixes:
snps_prefix=killerwhale4_snps.ID
bcftools_path=/home/cailak/programs/bcftools-1.9
gatk_jar_path=/home/cailak/programs/gatk-4.1.9.0
plink_path=/home/cailak/programs
ref_genome_path=/home/cailak/ref_genome
ref_genome_prefix=KW_GCA_937001465.1_mOrcOrc1.1.scafname
min_scaf_list=mOrcOrc_scafmin_100kb #list of scaffolds to keep



######## Prior to SNP filtering: Remove duplicates and close kin 

# Need a .txt file with list of individuals IDs to remove.
# Remove duplicates from VCF file
vcftools --gzvcf $snps_prefix.vcf.gz --remove duplicates_to_remove.txt --recode --recode-INFO-all --out $snps_prefix.dups_removed

# Remove close kin from VCF file
vcftools --vcf $snps_prefix.dups_removed.recode.vcf --remove closekin_to_remove.txt --recode --recode-INFO-all --out $snps_prefix.dups_removed.kin_removed

# Rename
mv $snps_prefix.dups_removed.kin_removed.recode.vcf $snps_prefix.dups_kin_removed.vcf

# Zip & index
bgzip $snps_prefix.dups_kin_removed.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.vcf.gz



######## Step 1: Remove indels and SNPs that don't have the "PASS" filter. Start with vcf.gz file (with corresponding tbi)

# Count snps
vcftools --gzvcf $snps_prefix.dups_kin_removed.vcf.gz --out $snps_prefix.dups_kin_removed.snpcount

# Filter out non-PASS sites
$bcftools_path/bcftools view -f PASS $snps_prefix.dups_kin_removed.vcf.gz -o $snps_prefix.dups_kin_removed.PASS.vcf

# Zip & index
bgzip $snps_prefix.dups_kin_removed.PASS.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.PASS.vcf.gz

# Count snps
vcftools --gzvcf $snps_prefix.dups_kin_removed.PASS.vcf.gz --out $snps_prefix.dups_kin_removed.PASS.snpcount



######## Step 2: Quality filtering
# Need reference genome .dict file in same directory as reference fasta
# Create .dict file - only needs to be done once.
#java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar CreateSequenceDictionary -R $ref_genome_path/$ref_genome_prefix.fasta -O $ref_genome_path/$ref_genome_prefix.dict 

# Filter out QUAL < 50
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar VariantFiltration -R $ref_genome_path/$ref_genome_prefix.fasta -V $snps_prefix.dups_kin_removed.PASS.vcf.gz -O $snps_prefix.dups_kin_removed.PASS.QUAL.vcf.gz --filter-name "QUALlt50" --filter-expression "QUAL < 50" 

# Count snps
vcftools --gzvcf $snps_prefix.dups_kin_removed.PASS.QUAL.vcf.gz --out $snps_prefix.dups_kin_removed.PASS.QUAL.snpcount

# Filter out MQ < 40
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar SelectVariants -R $ref_genome_path/$ref_genome_prefix.fasta -V $snps_prefix.dups_kin_removed.PASS.QUAL.vcf.gz -O $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.vcf.gz -select "MQ >= 40.0"

# Count snps
vcftools --gzvcf $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.vcf.gz --out $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.snpcount

# Filter out QD < 4
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar SelectVariants -R $ref_genome_path/$ref_genome_prefix.fasta -V $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.vcf.gz -O $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.QD.vcf.gz -select "QD >= 4.0"

# Count snps
vcftools --gzvcf $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.QD.vcf.gz --out $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.QD.snpcount

# Rename with "filter1" and clean up files
mv $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.QD.vcf.gz $snps_prefix.dups_kin_removed.filter1.vcf.gz
mv $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.QD.vcf.gz.tbi $snps_prefix.dups_kin_removed.filter1.vcf.gz.tbi

rm $snps_prefix.dups_kin_removed.PASS.QUAL.vcf.gz $snps_prefix.dups_kin_removed.PASS.QUAL.vcf.gz.tbi $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.vcf.gz $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.vcf.gz.tbi $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.QD.vcf.gz $snps_prefix.dups_kin_removed.PASS.QUAL.MQ.QD.vcf.gz.tbi



######## Step 3: Missingness & bi-allelic filter. These two steps should automatically create the log files for snp counts

# Filter out snps with high missingness (in vcftools 1=no missing, 0=all missing)
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --out $snps_prefix.dups_kin_removed.filter1.miss

# Remove non-biallelic sites
vcftools --vcf $snps_prefix.dups_kin_removed.filter1.miss.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $snps_prefix.dups_kin_removed.filter1.miss.biallel

mv $snps_prefix.dups_kin_removed.filter1.miss.biallel.recode.vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.vcf

# Zip & index
bgzip $snps_prefix.dups_kin_removed.filter1.miss.biallel.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.vcf.gz



######## Step 4: Remove scaffolds less than 100kb in length. Needs input list of scaffold names to KEEP. I chose 100kb based on my dataset, but should be checked/adjusted for other data (example, if this filter loses a very large chunk of the genome, may need to use a smaller scaffold filter like 10kb).

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $min_scaf_list | paste -s -d, - > scaf_list_line

# Set up list
list=`cat scaf_list_line`

# Filter vcf for these scaffolds
$bcftools_path/bcftools filter --regions $list $snps_prefix.dups_kin_removed.filter1.miss.biallel.vcf.gz > $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.vcf

# Can make a list of scaffolds to double check
grep -v "^#" $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.vcf | cut -f1 | sort | uniq > filtered_contig_list_check.txt

# Zip & index
bgzip $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.vcf.gz

# Count snps
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.vcf.gz --out $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.snpcount



######## Step 5: Filter out sex-linked SNPS to create autosomal dataset

# Do the filter to remove snps on X and Y
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.vcf.gz --not-chr OW443365.1 --recode --recode-INFO-all --out $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes
mv $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.recode.vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.vcf

# Zip and index
bgzip $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.vcf.gz



######## Step 6: HWE, MAF filter

# 6.1: Make separate VCF files for each population

# Need a list of the samples in each population (ECAG1 and ECAG2).

# MAKE ECAG1 FILE: 
# Make directory for separate HWE filtering.
mkdir HWE_filtering_ECAG1

# Remove ECAG2 individuals: 
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.vcf.gz --remove ECAG2_sampleIDs.txt --recode --recode-INFO-all --out HWE_filtering_ECAG1/$snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG1_only

# Zip and index
bgzip HWE_filtering_ECAG1/$snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG1_only.recode.vcf
tabix -p vcf HWE_filtering_ECAG1/$snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG1_only.recode.vcf.gz


# MAKE ECAG2 FILE: 
# Make directory for separate HWE filtering.
mkdir HWE_filtering_ECAG2

# Remove ECAG1 individuals: 
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.vcf.gz --remove ECAG1_sampleIDs.txt --recode --recode-INFO-all --out HWE_filtering_ECAG2/$snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG2_only

# Zip and index
bgzip HWE_filtering_ECAG2/$snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG2_only.recode.vcf
tabix -p vcf HWE_filtering_ECAG2/$snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG2_only.recode.vcf.gz


# 6.2: Create out.hwe file for both populations and run filter.hets.R script to create list of snps with het > 0.6

# Copy filter.hets.R file to both folders
cp filter.hets.R /home/cailak/new_kw_snp_filtering/pop_structure/HWE_filtering_ECAG1
cp filter.hets.R /home/cailak/new_kw_snp_filtering/pop_structure/HWE_filtering_ECAG2

# ECAG1
cd /home/cailak/new_kw_snp_filtering/pop_structure/HWE_filtering_ECAG1

# Create out.hwe file 
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG1_only.recode.vcf.gz --hardy

# create list of snps with het > 0.6
R --vanilla < filter.hets.R

# ECAG2
cd /home/cailak/new_kw_snp_filtering/pop_structure/HWE_filtering_ECAG2

# Create out.hwe file 
vcftools --gzvcf k$snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.ECAG2_only.recode.vcf.gz --hardy

# create list of snps with het > 0.6
R --vanilla < filter.hets.R


##### Step 6.3: Combine filter files and filter out SNPS

# Return to main folder
cd /home/cailak/new_kw_snp_filtering/pop_structure

# Combine snps_het06_filter files from ECAG1 and ECAG2
cat HWE_filtering_ECAG1/snps_het06_filter HWE_filtering_ECAG2/snps_het06_filter > snps_het06_filter_combined

# Exclude positions
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.vcf.gz --exclude-positions snps_het06_filter_combined --recode --recode-INFO-all --out $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe

# Rename
mv $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.recode.vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.vcf

# Zip & Index
bgzip $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz

# Minor-allele frequency filter
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz --maf 0.05 --recode --recode-INFO-all --out $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf

mv $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.recode.vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf

bgzip $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz



######## Step 7: LD pruning

# Create list of snps to prune
# --double-id b/c multiple "_" in sample name
$plink_path/plink --vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --allow-extra-chr --indep-pairwise 50 5 0.8 --set-missing-var-ids @:#\$1,\$2 --out pruned.50window.r0.8 --double-id

# Preparing the prune-in file to filter vcf
sed 's/:/\t/g' pruned.50window.r0.8.prune.in > temp.in
sed 's/...$//' temp.in > temp2.in
mv temp2.in list.pruned.50window.r0.8.prune.in
rm temp.in

# Prune the snps
vcftools --gzvcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --positions list.pruned.50window.r0.8.prune.in --recode --recode-INFO-all --out $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

# Rename
mv $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.recode.vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf

# Zip & Index
bgzip $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf
tabix -p vcf $snps_prefix.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf.gz
 
