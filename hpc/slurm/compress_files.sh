#!/bin/bash -l
# vim: syntax=sh tabstop=2 expandtab
# coding: utf-8

#SBATCH --mail-type ALL 
#SBATCH --mail-user firstname.surname@unil.ch 

#SBATCH --chdir /scratch/<your_username>/
#SBATCH --job-name xz
#SBATCH --output=xz_%A_%a.slurm.out
#SBATCH --error=xz_%A_%a.slurm.err

#SBATCH --partition cpu
#SBATCH --ntasks 1

#SBATCH --cpus-per-task 8
#SBATCH --mem 32G 
#SBATCH --time 04:00:00 
#SBATCH --export NONE

#SBATCH --array=0-99

#module load gcc/9.3.0 python/3.8.8

FILES=(/path_to_files/*)

# How many elements per job?
CHUNK_SIZE=20
# Start index in list
IDX_START=$(( $SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))

echo "$SLURM_ARRAY_TASK_ID"

for f in ${FILES[@]:$IDX_START:$CHUNK_SIZE}; do
  echo "xz --compress -9 --threads $SLURM_CPUS_PER_TASK $f"
done
