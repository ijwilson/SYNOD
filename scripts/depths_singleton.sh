#!/bin/bash 

working_dir=${HOME}/git/delocalise
results_dir=singleton_results
scratch_dir="$(mktemp -d)"

trap 'rm -rf -- "$scratchdir"' EXIT

set -e      ## script finishes on error 
cramfile=$1 ## cramfile given by command line argument
samtools quickcheck $1  ## check the file

reffile=${working_dir}/data/GRCh38_full_analysis_set_plus_decoy_hla.fa ## downloaded copy of reference
region_file=${working_dir}/data/100_regions.bed

base=$(basename $cramfile .final.cram)
bamfile=${scratch_dir}/${base}.bam
delfiles=("$bamfile"  "${bamfile}.bai" "${base}.final.cram.crai" )  

samtools view \
	    -b -T $reffile \
	    --region-file $region_file \
	    --threads 6 \
	    -o $bamfile \
	    $cramfile  
samtools index ${bamfile}

sambamba depth region --regions=$region_file --nthreads 6 \
                 ${bamfile} > ${working_dir}/${results_dir}/${base}.depth 

rm ${delfiles[@]}
