#!/bin/bash
#SBATCH --job-name=daemon_job_test    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cuffaro.m@ufl.edu # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=daemon_test_%j.log   # Standard output and error log
pwd; hostname; date

module load libcurl/7.49.1 gcc/12.2.0 openmpi/4.1.5 julia

echo "Running some tests!?! webhook updated∑!!!"

srun --mpi=$HPC_PMIX ./job.sh

date


