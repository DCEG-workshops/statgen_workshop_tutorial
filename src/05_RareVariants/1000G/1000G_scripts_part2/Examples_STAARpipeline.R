###########################################################
#    Run STAARpipeline 
#    Date: 20231105  
###########################################################

rm(list=ls())
gc()

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)

## set directories
input_dir <- "/content/drive/MyDrive/session05_Rare_variants/"
analysis_dir <- "/content/analysis_dir05/" 

## gds file
agds_dir <- get(load(paste0(analysis_dir, "agds_dir.Rdata")))
## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load(paste0(analysis_dir, "Annotation_name_catalog.Rdata")))

## channel name of the QC label in the GDS/aGDS file
QC_label <- "annotation/info/QC_label"
## variant type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## genes info
genes <- genes_info
## phenotype
pheno_cov <- read.table(paste0(input_dir, "1000G_data/integrated_call_samples_v3.20130502.ALL.panel"),header=TRUE)
phenotype.id <- as.vector(pheno_cov$sample)
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
## Annotation name (variants info)
Annotation_name_info <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## known loci: rs7412 (APOE), rs429358 (APOE), rs35136575 (APOE), rs12151108 (LDLR), rs688 (LDLR), rs6511720 (LDLR)
known_loci <- get(load(paste0(input_dir, "Reproducible_files/known_loci_info.Rdata")))


##########################################
#       load sparse GRM
##########################################

sgrm <- get(load("/content/drive/MyDrive/session05_Rare_variants/1000G_data/PC_sGRM/output.sparseGRM.sGRM.RData"))
sample_id <- unlist(lapply(strsplit(colnames(sgrm),"_"),`[[`,2))

colnames(sgrm) <- sample_id
rownames(sgrm) <- sample_id


### phenotype
pheno <- get(load(paste0(analysis_dir, "phenotype_LDLR_plof_ds_APOE_enhancer_DHS.Rdata")))
### PCs
PCs <- read.table(paste0(input_dir, "1000G_data/PC_sGRM/output.sparseGRM.score"),header=FALSE)
colnames(PCs)[2:12] <- c("id","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

pheno <- dplyr::left_join(pheno,PCs[,2:12],by=c("sample"="id"))

### fit null model
obj_nullmodel <- fit_null_glmmkin(Y~gender+super_pop+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
data = pheno, kins = sgrm, id = "sample", use_sparse = TRUE, family = gaussian(link = "identity"), verbose=T)


#####################################################
#        Gene-Centric Coding: LDLR, plof_ds 
#####################################################

### run coding mask of LDLR
gene_name <- "LDLR"

## genotype: chr 
chr <- 19
#gds.path <- agds_dir[chr]
gds.path <- paste0(input_dir, "1000G_data/aGDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds")
genofile <- seqOpen(gds.path)

results_coding <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
								QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

results_coding

## Conditional Analysis
category <- "plof_ds"								
results_coding_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,category=category,
                                known_loci=known_loci,rare_maf_cutoff=0.01,rv_num_cutoff=2,
								QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
								
results_coding_cond

## variants info
results_coding_info <- Gene_Centric_Coding_Info(category=category,chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,
										QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
										Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name_info)
		
dim(results_coding_info)
# [1] 17 33						
head(results_coding_info)

seqClose(genofile)
								
########################################################
#        Gene-Centric Noncoding: APOE, enhancer_DHS 
########################################################

### run noncoding mask of APOE
gene_name <- "APOE"

## genotype: chr 
chr <- 19
#gds.path <- agds_dir[chr]
gds.path <- paste0(input_dir, "1000G_data/aGDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds")
genofile <- seqOpen(gds.path)

results_noncoding <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
								QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

results_noncoding

## Conditional Analysis
category <- "enhancer_DHS"								
results_noncoding_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,category=category,
                                known_loci=known_loci,rare_maf_cutoff=0.01,rv_num_cutoff=2,
								QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

results_noncoding_cond
								
## variants info
results_noncoding_info <- Gene_Centric_Noncoding_Info(category=category,chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,
										QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
										Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name_info)

dim(results_noncoding_info)
# [1] 50 33						
head(results_noncoding_info)
								
seqClose(genofile)

########################################################
#        Sliding Window: 44,721,001-44,908,000 
########################################################

## genotype: chr 
chr <- 19
#gds.path <- agds_dir[chr]
gds.path <- paste0(input_dir, "1000G_data/aGDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds")
genofile <- seqOpen(gds.path)

start_loc_sub <- 44721001
end_loc_sub <- 44908000

a <- Sys.time()
results_sliding_window <- try(Sliding_Window(chr=chr, start_loc=start_loc_sub, end_loc=end_loc_sub, sliding_window_length = 2000, genofile=genofile, obj_nullmodel=obj_nullmodel, 
						type="multiple",QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
						Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
b <- Sys.time()
b - a

seqClose(genofile)

## alpha level
alpha <- 0.05/dim(results_sliding_window)[1]	
alpha
# [1] 0.0002688172
			
results_sig <- results_sliding_window[as.numeric(results_sliding_window[,90])<alpha,]

## conditional analysis
results_sig_cond <- c()
if(length(results_sig)!=0)
{
	for(kk in 1:dim(results_sig)[1])
	{
		chr <- as.numeric(results_sig[kk,1])
		start_loc <- as.numeric(results_sig[kk,2])
		end_loc <- as.numeric(results_sig[kk,3])

		#gds.path <- agds_dir[chr]
		gds.path <- paste0(input_dir, "1000G_data/aGDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds")
		genofile <- seqOpen(gds.path)
	
		res_cond <- Sliding_Window_cond(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,
				                                start_loc=start_loc,end_loc=end_loc,known_loci=known_loci,
				                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
				                                Annotation_name_catalog=Annotation_name_catalog,Annotation_dir=Annotation_dir,
				                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
		results_sig_cond <- rbind(results_sig_cond,res_cond)

		seqClose(genofile)
	}
}

## unconditional analysis results
results_sig[,c(1:4,90)]
## conditional analysis results
results_sig_cond[,c(1:4,90)]

## variants info
#gds.path <- agds_dir[chr]
gds.path <- paste0(input_dir, "1000G_data/aGDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds")
genofile <- seqOpen(gds.path)

start_loc <- 44900001
end_loc <- 44902000

results_sliding_window_info <- Sliding_Window_Info(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,start_loc=start_loc,end_loc=end_loc,known_loci=known_loci,
	QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
	Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name_info)

dim(results_sliding_window_info)
# [1] 32 32				
head(results_sliding_window_info)

seqClose(genofile)