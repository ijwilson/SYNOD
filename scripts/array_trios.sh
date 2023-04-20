#!/bin/bash 

#SBATCH --account=hcexsukb
#SBATCH --mail-type=NONE
#SBATCH --array=1-2%2
#SBATCH -e jobs/myjob_%A_%a.err
#SBATCH -o jobs/myjob_%A_%a.out

readarray list < remaining_trios3.txt

cramfiles=${list[$SLURM_ARRAY_TASK_ID-1]}

module load SAMtools

./depths_trio.sh $cramfiles $SLURM_ARRAY_TASK_ID




