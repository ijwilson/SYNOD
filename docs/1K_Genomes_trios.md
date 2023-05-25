# Getting trios from the 1000 Genomes

R code to get trios is in the file 
[scripts/find_all_trios.r](../scripts/find_all_trios.r).  This produces
a file `1000G_trios.txt` which contains the filenames and paths for CRAM
files for each trio, one trio per line in order: child, father, mother.

## Analysing Trios

The bash script [scripts/depth_inds.sh](scripts/depth_inds.sh) takes as inputs
paths to any number of CRAM files on the command line.  For analysing trios
we take as input a line from the "1000G_trios.txt" file.  

This script takes two extra options, the path to the sequence file used for assembly, and a bed file containing the regions over which to calculate the depth  

    Usage: depth_inds.sh [-r region_file] [-f reference_file] cram_1 ... cram_k


```bash
scripts/depth_inds.sh\
 -r ~/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa\
 -f data/100_regions.bed\
 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3988761/HG00405.final.cram \
 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241665/HG00403.final.cram \
 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241666/HG00404.final.cram

```

Analysis of trios is done in the file [scripts/analyse_all_trios.r](../scripts/analyse_all_trios.r).