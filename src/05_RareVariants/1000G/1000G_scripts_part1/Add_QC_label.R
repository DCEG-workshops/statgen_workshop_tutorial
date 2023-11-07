##########################################################
# Adds QC_label with all "PASS" to a post-QC GDS file
# Authors: Xihao Li, Zilin Li
##########################################################

library(gdsfmt)
library(SeqArray)

dir_geno <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/aGDS/"
gds_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
gds_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.gds"

for (chr in 1:22){
  print(paste("Chromosome:",chr))
  gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
  genofile<-seqOpen(gds.path, readonly = FALSE)
  #genofile
  
  position <- as.integer(seqGetData(genofile, "position"))
  length(position)
  Anno.folder <- index.gdsn(genofile, "annotation/info")
  add.gdsn(Anno.folder, "QC_label", val=factor(rep("PASS", length(position))), compress="LZMA_ra", closezip=TRUE)
  #genofile
  
  seqClose(genofile)
}

