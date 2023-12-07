#!/bin/bash
#
# this file is myjob.sh
#
#SBATCH --job-name Task2-toAlignFastaToReference
#SBATCH --mail-type BEGIN,END
#
# load minimap2/2.26 and samtools/1.17
module load minimap2/2.26
module load samtools/1.17

## Author: Oscar Florez-Vargas

##########################
## == Set Parameters == ##

genomicRegion="hg38_chr1_109655000_109742000"

## Set path to where BAM files are being stored
PATH_TO_FILES=/data/Prokunina_Group/GenomeAssemblies/HPRC_workshop/retrievedRegions
PATH_TO_REFERENCE=/data/Prokunina_Group/GenomeAssemblies/referenceData

##########################


# Create the "bam_files" folder
if [ ! -d $PATH_TO_FILES/$genomicRegion/bam_files ]; then
  mkdir -p $PATH_TO_FILES/$genomicRegion/bam_files
  chmod 775 $PATH_TO_FILES/$genomicRegion/bam_files;
fi


# Step 1: Concatenate all FASTA files in into input.fa
cat $PATH_TO_FILES/$genomicRegion/fasta_files/*.fa > $PATH_TO_FILES/$genomicRegion/bam_files/HPRC_input.fa


# Step 2: Align the FASTA sequences to reference human genome build hg38
minimap2 \
	-Y -t 6 -a \
	$PATH_TO_REFERENCE/genome_hg38.fa \
	$PATH_TO_FILES/$genomicRegion/bam_files/HPRC_input.fa | \
	samtools sort -@ 6 -o $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.bam

# Step 3: Index the BAM file
samtools index $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.bam

# Step 4: remove the temporal FASTA file
rm $OUTDIR/$genomicRegion/bam_files/HPRC_input.fa

