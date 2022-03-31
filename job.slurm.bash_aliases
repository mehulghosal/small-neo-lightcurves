#!/bin/bash
SBATCH --job-name=lightcurves_a
SBATCH --partition=
## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=0-01:00:00 ## time format is DD-HH:MM:SS


#SBATCH --cpus-per-task=1
#SBATCH --mem=6400 ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    

#SBATCH --error=hello-%A.err ## %A - filled with jobid
#SBATCH --output=hello-%A.out ## %A - filled with jobid

## Useful for remote notification
##SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
##SBATCH --mail-user=user@test.org

## All options and environment variables found on schedMD site: http://slurm.schedmd.com/sbatch.html

module purge
./build.sh
./hello ${RANDOM}