# Examining some stats on unfiltered snps/vcf to determine filtering thresholds

library(tidyverse)
library(readr)

setwd("~/Dropbox/killer_whale_genomics/snps4/snps_stats")

# Import unfiltered vcf information (.info file from bcftools)
vcfInfo <- read_table2("killerwhale4_snps.ID.info")
colnames(vcfInfo) <- c("CHROM", "POS", "REF", "ALT", "QUAL", "MQ", "QD", "TC", "FR")

# allele freq, missingness, hwe from vcftools
var_frq <- read_delim("killerwhale4_snps.ID.frq", delim="\t",col_names=c("CHR", "POS", "N_ALLELES", "N_CHR", "A1", "A2"), skip = 1)
var_miss <- read_delim("killerwhale4_snps.ID.lmiss", delim="\t")
#var_miss_updatedset <- read_delim("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.MISSINGNESS.lmiss", delim="\t")
ind_miss <- read_delim("killerwhale4_snps.ID.imiss", delim="\t")
#ind_miss_updatedset <- read_delim("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.MISSINGNESS.imiss", delim="\t")
hwe <- read_delim("killerwhale4_snps.ID.hwe", delim="\t", col_names=c("CHR", "POS", "OBS.HOM1.HET.HOM2", "E.HOM1.HET.HOM2", "ChiSq_HWE", "P_HWE", "P_HET_DEFICIT", "P_HET_EXCESS"), skip=1)

## View distributions. Red line setting filtering threshold
# QUAL
min(vcfInfo$QUAL)
max(vcfInfo$QUAL)

qual <- ggplot(vcfInfo, aes(x=QUAL))+
  geom_histogram(fill="gray60", bins=50)+
  theme_classic()+
  ggtitle("Quality Score")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
 # geom_vline(xintercept=20, col="red")+
  ylab("Count")

qual


# MQ
min(vcfInfo$MQ, na.rm=TRUE)
max(vcfInfo$MQ, na.rm=TRUE)

mq <- ggplot(vcfInfo, aes(x=MQ))+
  geom_histogram(fill="gray60", bins=40)+
  theme_classic()+
  ggtitle("Mapping Quality Score")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
#  geom_vline(xintercept=30, col="red")+
  ylab("Count")

mq


# QD
min(vcfInfo$QD, na.rm=TRUE)
max(vcfInfo$QD, na.rm=TRUE)
summary(vcfInfo$QD)

qd <- ggplot(vcfInfo, aes(x=QD))+
  geom_histogram(fill="gray60", bins=100)+
  theme_classic()+
  ggtitle("Quality By Depth")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
 # geom_vline(xintercept=2, col="red")+
  ylab("Count")

qd

# Minor allele frequency
var_frq$MAF <- var_frq %>% select(A1, A2) %>% apply(1, function(z) min(z))
summary(var_frq$MAF)

maf <- ggplot(var_frq, aes(x=MAF))+
  geom_histogram(fill="gray60", bins=80)+
  theme_classic()+
  ggtitle("Minor Allele Frequency")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
 # geom_vline(xintercept=0.041, col="red")+
  ylab("Count")

maf

# Variant missingness
variant_miss <- ggplot(var_miss, aes(x=F_MISS))+
  geom_histogram(fill="gray60", bins=50)+
  theme_classic()+
  ggtitle("Variant Missingness")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
 # geom_vline(xintercept=0.4, col="red")+
  xlab("Missingness")+
  ylab("Count")

variant_miss

# HWE for heterozygosity
summary(-log10(hwe$P_HWE))
hist(-log10(hwe$P_HWE), breaks=100, ylim=c(0,1000000), main="hwe p-value")

# Individual missingness
summary(ind_miss$N_MISS)

indiv_miss <- ggplot(data=ind_miss, aes(x=factor(INDV), y=F_MISS))+
  geom_bar(stat="identity", fill="gray60")+
  theme_classic()+
  xlab("Individual")+
  ylab("Missingness")+
  ggtitle("Individual Sample Missingness")+
  theme(plot.title=element_text(hjust=0.5, face="bold"),axis.text.x=element_blank())
 # geom_hline(yintercept=0.4, col="red")
#scale_fill_manual( values = c( "no"="gray", "yes"="orange" ), guide = FALSE )

indiv_miss


# If want to plot everything together in one figure:
library(patchwork)

(qual | mq ) / (qd | maf ) / (variant_miss | indiv_miss)

ggsave("KW4_SNPS_prefilter.png", width = 9, height = 9, dpi = 300)
