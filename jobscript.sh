#!/bin/bash
#SBATCH --job-name=daemon_job_test    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cuffaro.m@ufl.edu # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=daemon_test_%j.log   # Standard output and error log
pwd; hostname; date

module load gcc/12.2.0 openmpi/4.1.5 julia

echo "Running some docss!?! webhook updated!! ahh!! ah!!!"

julia --project=docs/ -e 'using Pkg; Pkg.status()'
julia --project=docs/ -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
# julia --project -e 'using Pkg; Pkg.status; Pkg.test()' > log_test.md
# ./job.sh

# mpiexec -np 1 julia --project -e test/runtests.jl

date


