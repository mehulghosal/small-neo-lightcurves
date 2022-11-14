#!/bin/bash
#SBATCH --job-name=star_analysis
#SBATCH --partition=shared
## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=00-08:00:00 ## time format is DD-HH:MM:SS

#SBATCH --cpus-per-task=1
#SBATCH --mem=6G ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    

#SBATCH --error="./output/$1.err"
#SBATCH --output="./errors/$1.out"

## Useful for remote notification
##SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80

#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=mehul@hawaii.edu

## All options and environment variables found on schedMD site: http://slurm.schedmd.com/sbatch.html

##python3 ./my_python_script $1 $2 $3 $4 $5

python3 ./new_magic.py $1
##echo $1
