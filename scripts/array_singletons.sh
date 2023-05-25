#!/bin/bash 

#SBATCH --account=XXXXX
#SBATCH --mail-type=NONE
#SBATCH --array=11-2600%6
#SBATCH -e s_jobs/mysing_%A_%a.err
#SBATCH -o s_jobs/mysing_%A_%a.out

readarray list < ../1000G_singletons.txt

cramfile=${list[$SLURM_ARRAY_TASK_ID-1]}
module load SAMtools/1.16.1-GCC-11.3.0

./depths_singleton.sh $cramfile


