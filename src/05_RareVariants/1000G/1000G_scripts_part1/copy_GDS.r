#############################################################
#     copy GDS files to aGDS file folder
#     2023/09/07
#############################################################

rm(list=ls())
gc()

### gds file
dir_gds <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/GDS/"
gds_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
gds_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.gds"

## agds file
dir_agds <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/aGDS/"

for(chr in 1:22)
{
	system(paste0("cp ",dir_gds,gds_file_name_1,chr,gds_file_name_2," ",dir_agds,gds_file_name_1,chr,gds_file_name_2))

}