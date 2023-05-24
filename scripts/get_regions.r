#!/usr/bin/env Rscript
# We need the rtracklayer library for R
# this is a Bioconductor package.  To install this you will need The R BiocManager package
# install.packages("BiocManager")
# BiocManager::install("rtracklayer")
suppressPackageStartupMessages(library(rtracklayer))
## GRCh38
working_dir <- file.path("c:/data", "git/SYNOD")
working_dir <- file.path("~", "git/SYNOD")

# The file below is 42MB so may take a while to download
tr <- import(file.path(
    "https://ftp.ensembl.org/pub/release-96/gtf",
    "homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz")
)
seqlevelsStyle(tr) <- "UCSC"
tr <- tr[tr$type == "gene" & seqnames(tr) %in% 
                               paste0("chr", c(1:22, "X", "Y", "MT"))]

df <- data.frame(seqnames(tr), 
                 start(tr), end(tr), id = tr$gene_id, gene_name=tr$gene_name)
write.table(df, file = file.path(working_dir, "data/gene_regions38.bed"),
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

## now try to get wide autosomal data
tr2 <- tr[width(tr) >= 8000 & seqnames(tr) %in% paste0("chr", 1:22)]
df2 <- data.frame(
    seqnames(tr2),
    start(tr2),
    end(tr2),
    id = tr2$gene_id,
    gene_name = tr2$gene_name
)
write.table(df2, file = file.path(working_dir, "data/wide_autosomal38.bed"),
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Now for 100 genes including NPHP1 & SLC66A1.  Silly way to do it!
# is in the git directory anyway
if (FALSE) {
  repeat {
    df3 <- df2[sort(sample(nrow(df2), 100)), ]
    if ("NPHP1" %in% df3$gene_name & "PQLC2" %in% df3$gene_name ) break
  }
  write.table(df3, file = file.path(working_dir, "data/100_regions.bed"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  system(paste0("chmod uga-w ", file.path(working_dir, "data/100_regions.bed")))
}

# plot

if (FALSE) {
   df <- read.table("data/gene_regions38.bed")
  colnames(df) <- c("chrom","start","end","gene")
  library(ggplot2)
  df$chrom <- factor(df$chrom, levels=paste0("chr",c(1:22,"X","Y")))
  ggplot(df, aes(y = end-start+1, x=chrom, fill=chrom)) + geom_violin() + scal$
  ggsave(file = "violin.png")
}