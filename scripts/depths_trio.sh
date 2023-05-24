#!/bin/bash 

working_dir=${HOME}/git/SYNOD

set -e
samtools quickcheck $1 $2 $3 

cramfiles=($1 $2 $3)

reffile=${HOME}/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
region_file=${working_dir}/data/100_regions.bed
scratch_dir=/tmp/bam
bamfiles=()
delfiles=()   ## keeps track of files to delete at the end

mkdir -p $scratch_dir

for cramfile in ${cramfiles[@]}
do
    base=$(basename $cramfile .final.cram)
    echo $base
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
                  ${bamfiles[@]} > ${working_dir}/results/${base}.depth 

rm ${delfiles[@]} 
