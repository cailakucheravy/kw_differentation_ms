# looking at kinship estimates

# Set working directory and load libraries
setwd("~/Dropbox/killer_whale_genomics/snps4/pop_structure_filtering")

library(plinkQC)
library(ggplot2)

# Set directory and file names for plinkQC:
directory="~/Dropbox/killer_whale_genomics/snps4/pop_structure_filtering/filtered_snps"
files="killerwhale4_snps.ID.dups_kin_removed.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08" ##for .genome, .imiss, .bin, .bed, .fam files
  
# Quick plot
#png("plinkQC_BOW.png", w=1000, h=500, res=100)
evaluate_check_relatedness(
  qcdir=directory,
  name=files,
  highIBDTh = 0.1875,
  imissTh = 0.03,
  interactive = FALSE,
  verbose = FALSE
)
#dev.off()
 
# Extract values for the related individuals from the .genome file (filtering by pi_hat values)
IBD <- read.table(paste(directory,"/", files,".genome", sep=""), header=T)
IBD_order_PH <- IBD[order(IBD[,10], decreasing=TRUE),]
write.csv(IBD_order_PH, "IBD_KW_PIHAT_round4_dups_ECAG1closekin_removed.csv")

# Make nicer histogram

# Order by pi_hat and then add nrow (for pair#s)
IBD_2 <- IBD[order(IBD[,10], decreasing=TRUE),]
IBD_2$pair <- 1:nrow(IBD_2)

ggplot(data=IBD_2, aes(x=PI_HAT))+
  geom_histogram(bins=125, color="steelblue", fill="steelblue")+
  theme_bw()+
  ylab("number of pairs")+
  xlab("pi hat")
#ggsave("IBD_KW_geomhistogram.png", width=6,height=4,dpi=300)


# PI_HAT value reference:
# 1st degree relative: 0.5
# 2nd degree relative: 0.25
# 3rd degree relative: 0.125

