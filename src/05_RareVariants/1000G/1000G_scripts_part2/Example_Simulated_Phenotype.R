###########################################################
#    Simulated Phenotype
#    Date: 20231106  
###########################################################

rm(list=ls())
gc()

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(STAAR)
library(STAARpipeline)

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


#######################################################
#             Chr 19, LDLR, plof_ds
#######################################################

## gene name
gene_name <- "LDLR"

## genotype: chr 
chr <- 19
#gds.path <- agds_dir[chr]
gds.path <- paste0(input_dir, "1000G_data/aGDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds")
genofile <- seqOpen(gds.path)

############ genotype 
## get SNV id, position, REF, ALT (whole genome)
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant")
{
	SNVlist <- filter == "PASS"
}

if(variant_type=="SNV")
{
	SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel")
{
	SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

position <- as.numeric(seqGetData(genofile, "position"))
variant.id <- seqGetData(genofile, "variant.id")

rm(filter)
gc()


### Gene
kk <- which(genes[,1]==gene_name)

sub_start_loc <- genes[kk,3]
sub_end_loc <- genes[kk,4]

is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
variant.id.gene <- variant.id[is.in]

seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

## plof_ds
## Gencode_Exonic
GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
## Gencode
GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
## Meta.SVM.Pred
MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

variant.id.gene <- seqGetData(genofile, "variant.id")
lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
variant.id.gene <- variant.id.gene[lof.in.plof]

seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

## genotype id
id.genotype <- seqGetData(genofile,"sample.id")
# id.genotype.match <- rep(0,length(id.genotype))

id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
phenotype.id.merge <- data.frame(phenotype.id)
phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
id.genotype.match <- phenotype.id.merge$index

## Genotype
Geno <- seqGetData(genofile, "$dosage")
Geno <- Geno[id.genotype.match,,drop=FALSE]

## impute missing
if(!is.null(dim(Geno)))
{
	if(dim(Geno)[2]>0)
	{
		if(geno_missing_imputation=="mean")
		{
			Geno <- matrix_flip_mean(Geno)$Geno
		}
		if(geno_missing_imputation=="minor")
		{
			Geno <- matrix_flip_minor(Geno)$Geno
		}
	}
}

Geno_LDLR <- Geno
AF_LDLR <- apply(Geno_LDLR,2,mean)/2
MAF_LDLR <- pmin(AF_LDLR,1-AF_LDLR)

############ LDLR effect
c0_LDLR <- 0.25
beta_LDLR <- -c0_LDLR * log10(MAF_LDLR)

LDLR_Effect <- Geno_LDLR %*% beta_LDLR

seqClose(genofile)


#######################################################
#             Chr 19, APOE, enhancer_DHS
#######################################################

## gene name
gene_name <- "APOE"

## genotype: chr 
chr <- 19
#gds.path <- agds_dir[chr]
gds.path <- paste0(input_dir, "1000G_data/aGDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds")
genofile <- seqOpen(gds.path)

############ genotype 
## Enhancer
varid <- seqGetData(genofile, "variant.id")

# Now extract the GeneHancer with rOCRs Signal Overlay
genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
genehancer <- genehancerAnno!=""

rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
rOCRs <- rOCRsAnno!=""
rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])

# Variants that covered by whole GeneHancer without rOCRs overlap
genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
enhancervpos <- as.numeric(seqGetData(genofile,"position"))
enhancervref <- as.character(seqGetData(genofile,"$ref"))
enhancervalt <- as.character(seqGetData(genofile,"$alt"))
dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant")
{
	SNVlist <- filter == "PASS"
}

if(variant_type=="SNV")
{
	SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel")
{
	SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV <- variant.id[SNVlist]

dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

seqResetFilter(genofile)

rm(dfHancerrOCRsVarGene)
gc()

### Gene
is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
variant.is.in <- variant.id.SNV[is.in]

seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

## genotype id
id.genotype <- seqGetData(genofile,"sample.id")
# id.genotype.match <- rep(0,length(id.genotype))

id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
phenotype.id.merge <- data.frame(phenotype.id)
phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
id.genotype.match <- phenotype.id.merge$index

## Genotype
Geno <- seqGetData(genofile, "$dosage")
Geno <- Geno[id.genotype.match,,drop=FALSE]

## impute missing
if(!is.null(dim(Geno)))
{
	if(dim(Geno)[2]>0)
	{
		if(geno_missing_imputation=="mean")
		{
			Geno <- matrix_flip_mean(Geno)$Geno
		}
		if(geno_missing_imputation=="minor")
		{
			Geno <- matrix_flip_minor(Geno)$Geno
		}
	}
}

Geno_APOE <- Geno
AF_APOE <- apply(Geno_APOE,2,mean)/2
MAF_APOE <- pmin(AF_APOE,1-AF_APOE)

Geno_APOE <- Geno_APOE[,MAF_APOE<0.01]
MAF_APOE <- MAF_APOE[MAF_APOE<0.01]


############ APOE effect
c0 <- 0.22
beta_APOE <- -c0 * log10(MAF_APOE)

set.seed(6)
# Generate two covariates
beta_dir_APOE <- 2*(runif(length(beta_APOE),0,1)<0.5) - 1
beta_APOE <- beta_APOE*beta_dir_APOE

APOE_Effect <- Geno_APOE %*% beta_APOE

seqClose(genofile)


########################################################
#       Simulated Phenotype
########################################################

alpha0 <- 0
alpha1 <- 0.5

set.seed(666)
# Generate two covariates
N <- length(phenotype.id)
sex <- pheno_cov$gender=="female"

# Generate error distributions
eps <- rnorm(N)

Y <- alpha0 + alpha1 * sex + LDLR_Effect + APOE_Effect + eps

pheno <- cbind(pheno_cov,Y)

save(pheno,file=paste0(analysis_dir, "phenotype_LDLR_plof_ds_APOE_enhancer_DHS.Rdata"))
