#!/bin/bash
#
# this file is myjob.sh
#
#SBATCH --job-name Task2-toAlignFastaToReference
#SBATCH --mail-type BEGIN,END
#
# load minimap2/2.26 and samtools/1.17 (uncomment the next two lines if it's run on Biowulf)
#module load minimap2/2.26
#module load samtools/1.17

## Author: Oscar Florez-Vargas

##########################
## == Set Parameters == ##

genomicRegion="hg38_chr1_109655000_109742000"

## Set path to where BAM files are being stored
PATH_TO_FILES=${analysis_dir}/retrievedRegions
PATH_TO_REFERENCE=${input_dir}

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
	-Y -t1 -a -K5M -I2G \
	$PATH_TO_REFERENCE/hg38_chr1.fa \
	$PATH_TO_FILES/$genomicRegion/bam_files/HPRC_input.fa | \
	samtools sort -@ 6 -o $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.bam

# Step 3: Index the BAM file
samtools index $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.bam

# Step 4: remove the temporal FASTA file
rm $PATH_TO_FILES/$genomicRegion/bam_files/HPRC_input.fa

