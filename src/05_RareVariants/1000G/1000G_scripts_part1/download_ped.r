#############################################################
#     Download 1000G pedigree file
#     2023/09/04
#############################################################

rm(list=ls())
gc()

download.file(url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt", 
destfile="/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/1kGP.3202_samples.pedigree_info.txt")
