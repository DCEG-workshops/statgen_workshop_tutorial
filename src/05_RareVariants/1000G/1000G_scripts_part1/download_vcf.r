#############################################################
#     Download 1000G WGS Data in VCF
#     2023/09/04
#############################################################

rm(list=ls())
gc()

options(download.file.method="libcurl", url.method="libcurl",timeout = max(6000, getOption("timeout")))

######### Download VCF
## Input
input_url <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
input_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
input_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

## Output
output_dir <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/VCF/"
output_file_name_1 <- input_file_name_1
output_file_name_2 <- input_file_name_2


for(chr in 1:22)
{
	print(chr)
	download.file(url=paste0(input_url,input_file_name_1,chr,input_file_name_2), 
	destfile=paste0(output_dir,output_file_name_1,chr,output_file_name_2))
}

######### Download tbi
## Input
input_url <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
input_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
input_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"

## Output
output_dir <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/VCF/"
output_file_name_1 <- input_file_name_1
output_file_name_2 <- input_file_name_2


for(chr in 1:22)
{
	print(chr)
	download.file(url=paste0(input_url,input_file_name_1,chr,input_file_name_2), 
	destfile=paste0(output_dir,output_file_name_1,chr,output_file_name_2))
}



