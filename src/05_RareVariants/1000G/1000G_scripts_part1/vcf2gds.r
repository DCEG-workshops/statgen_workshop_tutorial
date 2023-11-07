#############################################################
#     convert 1000G WGS Data VCF files to GDS files 
#     2023/09/04
#############################################################

rm(list=ls())
gc()

library("SeqArray")

chr <- as.numeric(commandArgs(TRUE)[1])

## input
input_dir <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/VCF/"
input_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
input_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

## output_dir
output_dir <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/GDS/"
output_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
output_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.gds"


vcf.fn <- paste0(input_dir,input_file_name_1,chr,input_file_name_2)
out.fn <- paste0(output_dir,output_file_name_1,chr,output_file_name_2)

# modify the header 
h <- seqVCF_Header(vcf.fn)
# h$info
h$info$Number[h$info$ID=="SOURCE"] <- "."

seqVCF2GDS(vcf.fn, out.fn, header = h, genotype.var.name = "GT", info.import=NULL, fmt.import=NULL, 
ignore.chr.prefix="chr", raise.error=TRUE, verbose=TRUE)


