#!/bin/bash 

working_dir=${HOME}/git/delocalise
scratch_dir="$(mktemp -d /tmp/BAMfiles_XXX)"
echo "files going in ${scratch_dir}"

trap 'rm -rf -- "$scratch_dir" ' EXIT   ## if exits with error, remove temp files

cramfiles=( "$@" )

set -e
samtools quickcheck $@      ## check to see if files exist and are cram files

reffile=${working_dir}/data/GRCh38_full_analysis_set_plus_decoy_hla.fa
region_file=${working_dir}/data/100_regions.bed
bamfiles=()
delfiles=()   ## keeps track of files to delete at the end

mkdir -p ${working_dir}/test_results/

for cramfile in ${cramfiles[@]}
do
    base=$(basename $cramfile .final.cram)
    echo "converting" $base
    bamfile=${scratch_dir}/${base}.bam
    bamfiles+=("$bamfile")
    delfiles+=("$bamfile"  "${bamfile}.bai" "${base}.final.cram.crai" )  

    samtools view \
	    -b -T $reffile \
	    --region-file $region_file \
	    --threads 6 \
	    -o $bamfile \
	    $cramfile  
    samtools index ${bamfile}
done

base=$(basename ${cramfiles[0]} .final.cram)

sambamba depth region --regions=$region_file --nthreads 6 \
                ${bamfiles[@]} > ${working_dir}/test_results/${base}.depth 

rm ${delfiles[@]} 
