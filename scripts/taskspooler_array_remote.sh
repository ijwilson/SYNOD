#!/bin/bash 

echo "read $# arguments"

if (($# != 1))
then 
	echo "We need an argument with a list of bam file names or locations"
	exit 1
fi

results_dir=/home/nijw/data/depth
mkdir -p $results_dir
reffile=/mnt/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
region_file=/home/nijw/BB/del-lite/working/gene_regions38.bed

tsp -S 3

readarray list < $1

for cramfile in ${list[@]};do
	ff=$(basename -- $cramfile)
	filestem="${ff%.*}"
	bamfilename=$filestem".bam"
	depthfile=$results_dir/$filestem".depth" 
	echo "Going to convert $cramfile to $bamfilename "
	echo "and then get the depth and put results into $depthfile "

	TEMPD=$(mktemp -d -p /mnt/data/tmp )

	tsp -L $filestem.dl  samtools view -b -T $reffile \
 	--threads 4 -o $TEMPD/$bamfilename $cramfile 

	tsp -L $filestem.index -d samtools index $TEMPD/$bamfilename

	tsp -L $filestem.depth -d sambamba depth region --regions=$region_file -o $depthfile --nthreads 4 $TEMPD/$bamfilename 

	tsp -L $filestem.rm -d rm -r $TEMPD
done
