#!/bin/bash 
set -e
# Check for the required arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 [-r region_file] [-f reference_file] cram_1 ... cram_k"
  exit 1
fi

reference_file="${HOME}/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
region_file="data/100_regions.bed"

# Read in the required arguments
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -r|--region-file)
    region="$2"
    shift 2 # past argument past value
    ;;
    -f|--reference-file)
    ref="$2"
    shift 2 #past argument past value
    ;;
    *)
    echo "all options arguments read"
    break
    ;;
esac
done

cramfiles=( "$@" )

# A couple of variables to get paths correct
SCRIPT_DIR_REL=$(dirname ${BASH_SOURCE[0]})
SCRIPT_DIR_ABS=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
CURRENT_DIR=$(pwd)
# I first need a directory to work with.  By default I use the 
# current directory.  Change if you need
working_dir=${CURRENT_DIR}
results_dir=${working_dir}/results
mkdir -p $results_dir
# The directory for temporary files is /tmp but you may need to change for HPCs
TMP_DIR=/tmp
# create a directory for putting working files
scratch_dir="$(mktemp -d ${TMP_DIR}/BAMfiles_XXX)"
echo "temporary files going in ${scratch_dir}"

trap 'rm -rf -- "$scratch_dir" ' EXIT   ## if exits with error, remove temp files

samtools quickcheck ${cramfiles[@]}      ## check to see if files exist and are cram files

bamfiles=()
delfiles=()   ## keeps track of files to delete at the end

for cramfile in ${cramfiles[@]}
do
    base=$(basename $cramfile .final.cram)
    echo "converting" $base
    bamfile=${scratch_dir}/${base}.bam
    bamfiles+=("$bamfile")
    delfiles+=("$bamfile"  "${bamfile}.bai" "${base}.final.cram.crai" )  

    samtools view \
	    -b -T $reference_file \
	    --region-file $region_file \
	    --threads 6 \
	    -o $bamfile \
	    $cramfile  
    samtools index ${bamfile}
done

base=$(basename ${cramfiles[0]} .final.cram)

sambamba depth region --regions=$region_file --nthreads 6 \
                ${bamfiles[@]} > ${results_dir}/${base}.depth 

rm ${delfiles[@]} 
