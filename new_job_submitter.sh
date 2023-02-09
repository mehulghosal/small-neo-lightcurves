#!/bin/bash

# READ MY CSV HERE

# For i in row:
# sbatch ./my_python_executor.slurm $arg1 $arg2 $arg3 $arg4 $arg5

input="fits_files_all_chips.txt"
while read -r line
do
	echo $line
	sbatch ./new_magic.slurm $line
done <fits_files.txt
