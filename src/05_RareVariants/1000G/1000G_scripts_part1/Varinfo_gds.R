rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
file_DBsplit <- "/n/holylfs05/LABS/xlin/Lab/xihao_zilin/1000G/FAVORDB/FAVORdatabase_chrsplit.csv"
### Targeted GDS
dir_geno <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/GDS/"
gds_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
gds_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.gds"
### output
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/aGDS/Variant_Info/"

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function 
###########################################################################

### make directory
system(paste0("mkdir ",output_path,"chr",chr))

### R package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

### chromosome number
## read info
DB_info <- read.csv(file_DBsplit,header=TRUE)
DB_info <- DB_info[DB_info$Chr==chr,]

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
genofile <- seqOpen(gds.path)

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)

seqClose(genofile)

## Generate VarInfo
for(kk in 1:dim(DB_info)[1])
{
	print(kk)

	VarInfo <- VarInfo_genome[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
	VarInfo <- data.frame(VarInfo)
	
	write.csv(VarInfo,paste0(output_path,"chr",chr,"/VarInfo_chr",chr,"_",kk,".csv"),quote=FALSE,row.names = FALSE)
}

