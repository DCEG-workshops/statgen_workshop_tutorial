#############################################################
#     gds2bed
#     2023/09/05
#############################################################

rm(list=ls())
gc()

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(SNPRelate)
library(kinship2)

chr <- as.numeric(commandArgs(TRUE)[1])

## gds
gds_dir <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/GDS/"
gds_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
gds_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.gds"

gds.file <- paste0(gds_dir,gds_file_name_1,chr,gds_file_name_2)

## snpgds
snpgds_dir <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/BED/SNPGDS/"
snpgds_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
snpgds_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.gds"

snpgds.file <- paste0(snpgds_dir,snpgds_file_name_1,chr,snpgds_file_name_2)

## BED
bed_dir <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/BED/"
bed_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
bed_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel"

bed.file <- paste0(bed_dir,gds_file_name_1,chr,bed_file_name_2)

## input parameter
min.MAF <- 0.05
max.miss <- 0.05
removeSNPGDS <- TRUE


###############################################################
#          GDS2SNPGDS
###############################################################

gds <- seqOpen(gds.file)
nsamp <- length(seqGetData(gds, "sample.id"))

AC <- seqGetData(gds, "annotation/info/AC")
MAF <- AC/(2*nsamp)
MAF <- pmin(MAF,1-MAF)

AllFilter <- isSNV(gds,biallelic=TRUE) & MAF > 0.01

variant.id <- seqGetData(gds, "variant.id")
variant.sel<- variant.id[AllFilter]
seqSetFilter(gds, variant.sel=NULL, sample.sel=NULL, variant.id=variant.sel,verbose=TRUE)

seqGDS2SNP(gds, snpgds.file, verbose=TRUE)

seqClose(gds)


####################################################################
#            SNPGDS2BED
####################################################################


SNPgds <- snpgdsOpen(snpgds.file)
	
SNP.select <- snpgdsSelectSNP(SNPgds, maf=min.MAF, missing.rate=1-max.miss)
snpgdsGDS2BED(SNPgds, bed.file, snp.id=SNP.select, snpfirstdim=FALSE, verbose=TRUE)
snpgdsClose(SNPgds)
	
if(removeSNPGDS==TRUE)	
{
	unlink(snpgds.file)
}

#### Use SNP IDs in the format Chr:Pos_RefA/AltA
bim <- read.table(paste0(bed.file,".bim"),header=F)
bim[,2] <- paste0(bim[,1],":",bim[,4],"_",bim[,5],"/",bim[,6])

write.table(bim,paste0(bed.file,".bim"),row.names=F,col.names=F,quote=F)