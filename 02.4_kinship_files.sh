#!/bin/bash

# Note that this script can be run directly in terminal; does not need to be submitted as a job.

# Convert vcf file to .bin .bed .fam files 
/home/cailak/programs/plink --allow-extra-chr --double-id --make-bed --vcf killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf.gz --set-missing-var-ids @:#\$1,\$2 --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

# Kinship files: 
/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08 --genome --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08 --missing --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08



# For ECAG1 only:
# Convert vcf file to .bin .bed .fam files 
/home/cailak/programs/plink --allow-extra-chr --double-id --make-bed --vcf killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.ECAG1_only.recode.vcf --set-missing-var-ids @:#\$1,\$2 --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.ECAG1_only.recode

# Kinship files: 
/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.ECAG1_only.recode --genome --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.ECAG1_only.recode

/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.ECAG1_only.recode --missing --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.ECAG1_only.recode


# For no kin removed: 
# Convert vcf file to .bin .bed .fam files 
/home/cailak/programs/plink --allow-extra-chr --double-id --make-bed --vcf killerwhale4_snps.ID.dups_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf.gz --set-missing-var-ids @:#\$1,\$2 --out killerwhale4_snps.ID.dups_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

# Kinship files: 
/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08 --genome --out killerwhale4_snps.ID.dups_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08 --missing --out killerwhale4_snps.ID.dups_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

# For lfmm files: 
# Convert vcf file to .bin .bed .fam files 
/home/cailak/programs/plink --allow-extra-chr --double-id --make-bed --vcf killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.maf.sm.imputed.thin1000.vcf.gz --set-missing-var-ids @:#\$1,\$2 --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.maf.sm.imputed.thin1000

# Kinship files: 
/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.maf.sm.imputed.thin1000 --genome --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.maf.sm.imputed.thin1000

/home/cailak/programs/plink --allow-extra-chr --double-id --bfile killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.maf.sm.imputed.thin1000 --missing --out killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.maf.sm.imputed.thin1000




