
# Example Chromosome 11 for NA20752 from 1000G

To perform this example we need to create

* BAM files
* A set of genes to check.  See the [scripts/get_regions.r](../scripts/get_regions.r) file to get `wide_autosomal.bad`.
* [sambamba](https://lomereiter.github.io/sambamba/) and [samtools](https://samtools.github.io/bcftools/howtos/roh-calling.html) installed on a local system.

The instructions depend on using a linux system, or a windows linux subsystem [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

Note that under ubuntu both samtools and sambamba (and igv) can be installed using\
`sudo apt install samtools sambamba igv`.

## Setting up BAM files

Some BAM files are available to download files from the 1000 genomes,
but the main alignment files available are CRAM files, which are
compressed BAM files.

The latest versions of `sambamba` do not have support for CRAM files
and the files must be converted into bam files first.

We can proceed by downloading the cram files from the 1000 genomes
web site and then converting to bam.  To save space it is better to
download and convert to bam in one step and **samtools** has the ability
to work on remote files.  This is explored more fully in the file
[1000_genomes](docs/1000_genomes).

It is also preferable to have a local sequence file used to produce the cram files.  For the 1000 genomes data this is
[GRCh38_full_analysis_set_plus_decoy_hla.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa).
We also need the index file
[GRCh38_full_analysis_set_plus_decoy_hla.fa.fai](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai).

We assume that all this is being done in the directory `${HOME}/SYNOD`

### Download reference genome

```bash
cd ~/SYNOD
mkdir reference
cd reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
```

The script below gets depths for a subset of the whole ensembl version 96
for a sample from the 1000 genomes.

```bash
#!/bin/bash

# prepare bed file for chr11
grep chr11 ~/SYNOD/data/wide_autosomal38.bed > chr11.bed

# convert to bam file
samtools view -b -T ~/SYNOD/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa\
         --threads 4 -o NA20752_chr11.bam http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239832/NA20752.final.cram chr11
samtools index NA20752_chr11.bam
sambamba depth region --regions=chr11.bed NA20752_chr11.bam > NA20752_chr11.depth
```

Depending on the machine, this is likely to take about 15 minutes to run, and requires about
1.6 GB for the temporary BAM file.

The head of the output file is shown below

```text
# chrom chromStart      chromEnd        F3      F4      readCount       meanCoverage    sampleName
chr11   112967  125927  ENSG00000254468 AC069287.1      1181    13.3532 NA20752
chr11   127204  139612  ENSG00000230724 LINC01001       2102    24.717  NA20752
chr11   129279  186136  ENSG00000255229 AC069287.3      7750    20.1147 NA20752
chr11   167784  207428  ENSG00000177951 BET1L   8243    30.2585 NA20752
chr11   215030  236931  ENSG00000142082 SIRT3   5159    34.7775 NA20752
chr11   236966  252984  ENSG00000185627 PSMD13  3709    34.0194 NA20752
chr11   369499  382117  ENSG00000182272 B4GALNT4        3263    37.3986 NA20752
chr11   392614  404908  ENSG00000184363 PKP3    2848    32.2268 NA20752
chr11   405716  417455  ENSG00000185187 SIGIRR  2865    33.7015 NA20752
```

## Postprocessing

Postprocessing of these results can be done using R.  The method is really simple.  Divide the
mean coverage by the median mean coverage for the individual.

```r
aa <- read.table("NA20752_chr11.depth",
    col.names = c("chrom", "start", "end", "symbol","ensemblid","readCount", "meanCoverage", "ID"))
aa$relative_depth <- 2*aa$meanCoverage/median(aa$meanCoverage)
aa$cn <- cut(aa$relative_depth, breaks=c(0, 0.4,1.5,2.5,3.5,100), labels=0:4)

table(aa$cn)
```
We used the bioconductor package to plot the copy numbers

```r
library(karyoploteR)
g <- GRanges(aa)
kp <- plotKaryotype(genome="hg38",chromosomes=c("chr11"), plot.type = 1)
g$col <- c("red","magenta","lightgrey","cyan","green")[as.integer(g$cn)]
kpPlotRegions(kp, g[g$cn != "2"], col=g$col[g$cn != 2])
kpPoints(kp, data = g, y = g$relative_depth, ymax=5, ymin=0)
kpAxis(kp, ymin=0, ymax=5, numticks=6)
```

![Relative Depths in Chromosome 11](chr11.png "chromosome 11 depth")