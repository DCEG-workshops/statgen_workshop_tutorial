#!/bin/bash
#SBATCH -J 1kGM
#SBATCH -p test
#SBATCH --time=0-08:00
#SBATCH --mem=40000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-type=NONE

INPUT_PATH=/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/BED
OUTPUT_PATH=/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/BED

for i in {1..22}; do echo $INPUT_PATH/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel >> $INPUT_PATH/mergeBED.list; done;
/n/home05/zilinli/plink --merge-list mergeBED.list --make-bed --out ${OUTPUT_PATH}/genome