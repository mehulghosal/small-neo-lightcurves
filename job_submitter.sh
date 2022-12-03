#!/bin/bash

# READ MY CSV HERE

# For i in row:
# sbatch ./my_python_executor.slurm $arg1 $arg2 $arg3 $arg4 $arg5

input="ast_files.txt"
while read -r line
do
	echo $line
	sbatch ./ast_magic.slurm $line
done <ast_files.txt
