#############################################################
#     Download FAVOR DB
#     2023/09/04
#############################################################

rm(list=ls())
gc()

options(download.file.method="libcurl", url.method="libcurl",timeout = max(60000, getOption("timeout")))

## Output
output_dir <- "/n/holylfs05/LABS/xlin/Lab/xihao_zilin/1000G/FAVORDB/"
output_file_name_1 <- "chr"
output_file_name_2 <- ".tar.gz"

## Input
input_url_seq <- c("https://dataverse.harvard.edu/api/access/datafile/6170506",
"https://dataverse.harvard.edu/api/access/datafile/6170501",
"https://dataverse.harvard.edu/api/access/datafile/6170502",
"https://dataverse.harvard.edu/api/access/datafile/6170521",
"https://dataverse.harvard.edu/api/access/datafile/6170511",
"https://dataverse.harvard.edu/api/access/datafile/6170516",
"https://dataverse.harvard.edu/api/access/datafile/6170505",
"https://dataverse.harvard.edu/api/access/datafile/6170513",
"https://dataverse.harvard.edu/api/access/datafile/6165867",
"https://dataverse.harvard.edu/api/access/datafile/6170507",
"https://dataverse.harvard.edu/api/access/datafile/6170517",
"https://dataverse.harvard.edu/api/access/datafile/6170520",
"https://dataverse.harvard.edu/api/access/datafile/6170503",
"https://dataverse.harvard.edu/api/access/datafile/6170509",
"https://dataverse.harvard.edu/api/access/datafile/6170515",
"https://dataverse.harvard.edu/api/access/datafile/6170518",
"https://dataverse.harvard.edu/api/access/datafile/6170510",
"https://dataverse.harvard.edu/api/access/datafile/6170508",
"https://dataverse.harvard.edu/api/access/datafile/6170514",
"https://dataverse.harvard.edu/api/access/datafile/6170512",
"https://dataverse.harvard.edu/api/access/datafile/6170519",
"https://dataverse.harvard.edu/api/access/datafile/6170504")

for(chr in 1:22)
{
	print(chr)
	
	input_url <- input_url_seq[chr]
	download.file(url=paste0(input_url), 
	destfile=paste0(output_dir,output_file_name_1,chr,output_file_name_2))
}


### untar
library(devtools)
for(chr in 1:22)
{
	print(chr)
	
	untar(paste0(output_dir,"chr",chr,".tar.gz"), exdir=output_dir)

}
