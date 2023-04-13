---
output:
  pdf_document: default
  html_document: default
  word_document: default
---
# Calling CNVs in 1000 Genomes data

The 1000 Genomes Project was an international research effort that aimed to sequence the genomes of at
least 1,000 individuals from around the world to create a comprehensive catalogue of human genetic
variation. The project was initiated in 2008 and completed in 2015, and involved scientists from over
75 institutions in more than 20 countries.

The primary goal of the project was to create a reference database of genetic variation that could be
used to improve our understanding of the genetic basis of disease, as well as to aid in the
development of new diagnostic and therapeutic tools. By sequencing the genomes of individuals from
diverse populations, the project aimed to capture a broad range of genetic variation, including rare
variants that are difficult to detect with smaller sample sizes.  Here we show how these data can be
used to look at the frequency of a set of potential disease genes.  A advantage of the
1000 genomes data is that a number of the individuals have both parents sequenced.

The 1000 Genomes Project generated a vast amount of data, including more than 84 million genetic
variants across the genomes of 2,504 individuals from 26 populations around the world. This data has
been made publicly available and has been widely used by researchers in a variety of fields.

You can find more detailed information about the populations included in the 1000 Genomes Project, as
well as access to the project data, at the official website: <https://www.internationalgenome.org>.

In this document I describe how to find copy number variants in teh freely available 1000 genomes data.

## CRAM Files

The CRAM file format is a compressed file format for storing
Next-Generation Sequencing (NGS) read alignments. It was developed as a more efficient alternative to the popular BAM file format, which can become unwieldy for large-scale genomic datasets.

CRAM files use a combination of reference-based compression and lossless coding to achieve
significant reductions in file size, typically around 50-60% compared to BAM files. This can be
especially beneficial for large-scale genomic projects, where storage and transfer of large datasets
can be a challenge. They also allow for efficient retrieval of specific subsets of data, which can be helpful for analyses -- like this -- that focus on specific genomic regions.

### Accessing the 1000 genomes CRAM files

On the data section of the 1000 genomes home page, you can find links to the CRAM files for each of the 2,504 samples in
the project. You can also download metadata and other information about the project and individual samples from this
page.  [Data reuse policy](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20200526_1000G_2504plus698_high_cov_data_reuse_README.txt).

A list of all the files can be obtained from the *Available data* tab by clicking "download the list"
button at
<https://www.internationalgenome.org/data-portal/data-collection/30x-grch38>

This file is available [here](/data/all_1000G_samples.tsv).  The first 4 lines of this file are

|url|md5|Data collection| Data type  |Analysis group|Sample  |Population |Data reuse policy|
|---|---|---------------|------------|--------------|--------|-----------|------------------|
|<ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram>|923...|1000 Genomes 30x on GRCh38|alignment|PCR-free high coverage|NA12718|Utah residents (CEPH) with Northern and Western European ancestry|as above
|<ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239482/NA12775.final.cram>|619...|1000 Genomes 30x on GRCh38|alignment|PCR-free high coverage|NA12775|Utah residents (CEPH) with Northern and Western European ancestry|as above
|<ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239487/NA12842.final.cram>|174...|1000 Genomes 30x on GRCh38|alignment|PCR-free high coverage|NA12842|Utah residents (CEPH) with Northern and Western European ancestry|as above

## Converting CRAM to BAM files

Converting CRAM files to BAM files is necessary when working with sambamba which no longer supports the CRAM format. It
can also be helpful for data visualization and exploration, as the BAM format is widely supported by many genome browsers
and visualization tools.  To convert a CRAM file to a BAM file, you can use a tool such as Samtools or Picard.  Our examples
use samtools.

Samtools requires a reference sequence to be specified when converting the file format. The reference sequence can be
provided as a separate file or as a URL link to a reference genome database.  You also need Samtools installed on your
computer. If it's not already installed, you can download it from the [Samtools website](http://www.htslib.org/download/),
or it is available as a package on ubuntu (`apt install samtools`)

Open a terminal window and navigate to the directory where your CRAM file is located.

Use the Samtools "view" command to convert the CRAM file to a BAM file. The basic syntax for the command is as follows:

```bash
samtools view --bam --reference reference.fa --output output.bam input.cram
```

In this command, `reference.fa` is the name of the reference sequence file, and `input.cram` and `output.bam` are the
input and output file names, respectively.  The index for the reference file must also be available.

The reference sequence is used as a guide to reconstruct the original sequences of the reads and to ensure that the
alignments in the resulting BAM file are accurate. The reference sequence can be in FASTA format, and it should be the
same reference sequence that was used to create the original CRAM file.  The reference sequence used for the 1000
genomes data is at <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa>, along with its index GRCh38_full_analysis_set_plus_decoy_hla.fa.fai.

Depending on the size of your input file and the resources available on your computer, the conversion process may take some time to complete.

### Converting a subset of regions from the remote file

The following command converts regions given in the bed file `data/100_regions.bed`
from the CRAM file on the ebi web site into a local bam file.  It uses an online version of the reference file and the
reference file index.

```bash
samtools view \
   --reference http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   --threads 4 \
   --region-file ../data/100_regions.bed \
   --output HG00121.bam \
   http://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240195/HG00121.final.cram 
```

This can be slow, and is unreliable.  It can be made quicker and more stable by using a local copy of the reference and
reference index file, here downloaded into the data directory.

```bash
samtools view \
   --reference ../data/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   --threads 4 \
   --region-file ../data/100_regions.bed \
   --output HG00121.bam \
   http://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240195/HG00121.final.cram 
```

To work with this BAM file, we need an index which can be done using *samtools*.

```bash
samtools index HG00121.bam 
```

### Getting depths using sambamba

Sambamba is used to get the depths.  We only downloaded and transformed the bits of the original CRAM file in the file
100_regions.bed.  We get the mean depths over those regions.

```bash
mkdir -p ../test_results
sambamba depth region --regions=../data/100_regions.bed \
                      --output-filename ../test_results/HG00121.depth HG00121.bam 
```

### A complete script

The script in `scripts/depth_inds.sh`, download and transforms CRAM files given on the command line and get the depths, doing a quick test of each CRAM file, dealing with temporary files, and removing files on error.  The output is named after the first file.

```bash
./depth_inds.sh \
  http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram \
  http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239482/NA12775.final.cram \
  http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239487/NA12842.final.cram \
  http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239494/NA18536.final.cram
```

## Post-process in R

The script below analyses a single depth file as a command line argument.  It
depends on R being installed.  The script is in `/R/simple_analysis.r`.

```r
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
a <- read.table(args[1], skip = 1, sep = "\t",
           col.names = c("chrom", "start", "end", "geneid", "symbol",
                        "readCount", "meanCoverage", "sampleName"))

sp <- split(a, a$sampleName)
rel <- sapply(sp, function(xx) 2 * xx$meanCoverage / median(xx$meanCoverage))
cn <- apply(rel, 2,
            cut, breaks = c(-2, 0.4, 1.5, 2.5, 3.5, 500), labels = c(0:4))
rownames(rel) <- rownames(cn) <- sp[[1]]$symbol

cn["NPHP1", ]
```

Running this with the output from our depth_inds command:

```bash
R/simple_analysis.r test_results/NA12718.depth
```

gives results

```text
NA12718 NA12775 NA12842 NA18536 
   "2"     "2"     "2"     "2" 
```

NPHP1 has copy number 2 for all four files.
