#!/bin/bash
#SBATCH -J 1kGP
#SBATCH -p test
#SBATCH --time=0-08:00
#SBATCH --mem=40000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-type=NONE

INPUT_PATH=/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/BED
OUTPUT_PATH=/n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/BED

/n/home05/zilinli/plink --bfile $INPUT_PATH/genome --indep-pairwise 50 5 0.1 --out ${OUTPUT_PATH}/genome.prunedlist;
/n/home05/zilinli/plink --bfile ${OUTPUT_PATH}/genome --extract genome.prunedlist.prune.in --make-bed --out ${OUTPUT_PATH}/genome_pruned