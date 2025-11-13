#!/bin/bash -l

#SBATCH -p med                                                                  # setting medium priority
#SBATCH --job-name=snakemake                                                    # setting name of job
#SBATCH -c 1                                                                    # number of cores
#SBATCH --mem=8G                                                              # setting the memory per cpu
#SBATCH -t 5-0:00:00                                                            # setting the max time
#SBATCH -D /home/kslauck/projects/nestbox-temp-monitoring                                     # setting home directory
#SBATCH -e /home/kslauck/projects/nestbox-temp-monitoring/slurm_log/sterror_%j.txt            # setting standard error output
#SBATCH -o /home/kslauck/projects/nestbox-temp-monitoring/slurm_log/stdoutput_%j.txt          # setting standard output
#SBATCH --mail-type=BEGIN                                                       # mail alerts at beginning and end of job
#SBATCH --mail-type=END
#SBATCH --mail-user=kslauck@ucdavis.edu                                         # send mail here

# Standard output file contents
echo "## Script:"
cat Snakefile
echo "## Output:"

# initialize mamba
. ~/miniconda3/etc/profile.d/conda.sh

# activate mamba environment
conda activate base

# fail on weird errors
set -e
set -x

# clean up snakemake
snakemake --unlock

# Run script
snakemake -j 500 --slurm --profile .snakemake/slurm_profile/

  # Print out values of the current jobs SLURM environment variables
  env | grep SLURM

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
