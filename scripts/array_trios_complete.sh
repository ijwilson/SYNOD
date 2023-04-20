#!/bin/bash 

#SBATCH --account=hcexsukb
#SBATCH --mail-type=NONE
#SBATCH --array=99,192,197,206,219,220,279,286,295,308,309,400,403,417,537,538,542,543,558,595%4
#SBATCH -e jobs/myjobc_%A_%a.err
#SBATCH -o jobs/myjobc_%A_%a.out

readarray list < ../data/1000G_trios.txt

cramfiles=${list[$SLURM_ARRAY_TASK_ID-1]}

module load SAMtools

./depths_trio.sh $cramfiles




