#!/bin/bash
#SBATCH -J 1kG
#SBATCH -p shared,xlin,xlin-lab
#SBATCH --time=0-72:00
#SBATCH --array=1-22 --mem=60000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-type=NONE

module purge
module load gcc/9.5.0-fasrc01 openmpi/4.1.4-fasrc01
module load intel-mkl/23.0.0-fasrc01
module load cmake/3.25.2-fasrc01

export R_LIBS_USER=$HOME/R-4.1.0-MKL
echo $R_LIBS_USER

/n/home05/zilinli/R-4.1.0/bin/Rscript --slave --no-restore --no-save gds2agds.R ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout