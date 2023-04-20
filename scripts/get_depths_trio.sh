#!/bin/bash 

working_dir=${HOME}/git/delocalise

cramfiles=("http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram" \
              "http://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989341/NA12891.final.cram" \
              "http://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989342/NA12892.final.cram") 


reffile=${working_dir}/data/GRCh38_full_analysis_set_plus_decoy_hla.fa
region_file=${working_dir}/data/100_regions.bed
bamfiles=()
delfiles=()   ## keeps track of files to delete at the end

for cramfile in ${cramfiles[@]}
do
    base=$(basename $cramfile .final.cram)
    echo $base
    bamfile=${base}.bam
    bamfiles+=("$bamfile")
    delfiles+=("$bamfile"  "${bamfile}.bai" "${cramfile}.crai" )  

    samtools view \
	    -b -T $reffile \
	    --region-file $region_file \
	    --threads 6 \
	    -o $bamfile \
	    $cramfile  
    samtools index ${bamfile}
done

sambamba depth region --regions=$region_file --nthreads 6 \
                  ${bamfiles[@]} > ${working_dir}/results/trio1.depth 

rm ${delfiles[@]} 
