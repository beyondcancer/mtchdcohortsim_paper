#!/bin/bash
#SBATCH --job-name=HR_MATCHED_SIM_SS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lsh1804992@lshtm.ac.uk
#SBATCH --ntasks=1
#SBATCH --mem=200gb
#SBATCH --time=120:05:00
#SBATCH --array=1-40

module load stata
stata -b do jobhr${SLURM_ARRAY_TASK_ID}

