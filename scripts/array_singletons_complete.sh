#!/bin/bash 


#SBATCH --account=hcexsukb
#SBATCH --mail-type=NONE
#SBATCH --array=78,408,987,988,1033
#SBATCH -e s_jobs/mysingcomplete_%A_%a.err
#SBATCH -o s_jobs/mysingcomplete_%A_%a.out

readarray list < ../data/1000G_singletons.txt

cramfile=${list[$SLURM_ARRAY_TASK_ID-1]}
module load SAMtools/1.16.1-GCC-11.3.0


./depths_singleton.sh $cramfile


