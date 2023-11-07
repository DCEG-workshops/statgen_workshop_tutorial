rm(list=ls())
gc()

###############################################
#           Input
###############################################

## 1000G
dir.geno <- "/content/drive/MyDrive/session05_Rare_variants/1000G_data/aGDS/"
adgs_file_name_1 <- "1kGP_high_coverage_Illumina.chr"
agds_file_name_2 <- ".filtered.SNV_INDEL_SV_phased_panel.gds"
QC_label <- "annotation/info/QC_label"

output_path <- "/content/analysis_dir05/" 
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info","/genecode_comprehensive_exonic_category","/metasvm_pred","/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed","/apc_epigenetics_transcription",
         "/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")

################################################
#           Main Function 
################################################
### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

#### aGDS directory
agds_dir <- paste0(dir.geno,adgs_file_name_1,seq(19,19),agds_file_name_2) 
save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))

#### Annotation dir
Annotation_name_catalog <- data.frame(name=name,dir=dir)
save(Annotation_name_catalog,file=paste0(output_path,"Annotation_name_catalog.Rdata",sep=""))




