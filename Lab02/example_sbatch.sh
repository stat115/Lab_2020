#!/bin/bash
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 15 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=1G # Memory in GB (see also --mem-per-cpu)
#SBATCH -o output_%j.out # Standard out goes to this file
#SBATCH -e error_%j.err # Standard err goes to this file

# LOAD_MODULES
module load bwa/0.7.15-fasrc02


# YOUR_COMMANDS_HERE
echo JOB_FINISHED
