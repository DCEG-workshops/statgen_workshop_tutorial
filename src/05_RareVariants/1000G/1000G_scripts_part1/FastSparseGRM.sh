#!/bin/bash
#SBATCH -J 1kGP
#SBATCH -p test
#SBATCH --time=0-08:00
#SBATCH --mem=40000
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mail-type=NONE

module purge
module load gcc/9.5.0-fasrc01 openmpi/4.1.4-fasrc01
module load intel-mkl/23.0.0-fasrc01
module load cmake/3.25.2-fasrc01

export R_LIBS_USER=$HOME/R-4.1.0-MKL
echo $R_LIBS_USER

/n/home05/zilinli/R-4.1.0/bin/R CMD BATCH --vanilla '--args --prefix.in /n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/BED/genome_pruned --prefix.out /n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/PC_GRM/output.sparseGRM  --KING.executable ./king --degree 4 --num_threads 4 --tempDir /n/holystore01/LABS/xlin/Lab/xihao_zilin/1000G/PC_GRM' runPipeline_wrapper.R runPipeline.Rout