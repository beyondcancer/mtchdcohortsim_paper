#!/bin/bash
#SBATCH --job-name=TEST_HR_MATCHED_SIM_SS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lsh1804992@lshtm.ac.uk
#SBATCH --ntasks=1
#SBATCH --mem=200gb
#SBATCH --time=120:05:00
#SBATCH --array=1-600

module load stata
stata -b do install
stata -b do job_${SLURM_ARRAY_TASK_ID}

