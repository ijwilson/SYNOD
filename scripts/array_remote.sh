#!/bin/bash 

#SBATCH --partition=long
#SBATCH --account=hcexsukb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.wilson@newcastle.ac.uk

echo "read $# arguments"

if [$# -ne 1];
then 
	echo "We need an argument with a list of bam file names or locations"
fi


readarray list < $1
cramfile=${list[$SLURM_ARRAY_TASK_ID-1]}
ff=$(basename -- $cramfile)
filestem="${ff%.*}"
bamfilename=$filestem".bam"
results_dir=/nobackup/nijw/depth
depthfile=$results_dir/$filestem".depth" 
echo "Going to convert $cramfile to $bamfilename "
echo "and then get the depth and put results into $depthfile "

module load SAMtools/1.16.1-GCC-11.3.0
TEMPD=$(mktemp -d -p /scratch/nijw )

reffile=/nobackup/nijw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

samtools view -b -T $reffile \
              --threads 4 \
              -o $TEMPD/$bamfilename $cramfile

samtools index $TEMPD/$bamfilename

ls -l $TEMPD/$bamfilename

region_file=/nobackup/nijw/reference/gene_regions38.bed

sambamba depth region --regions=$region_file \
                      --nthreads 6 \
                       $TEMPD/$bamfilename \
  > $depthfile 


