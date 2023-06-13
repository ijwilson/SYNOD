# analyse trios
library(data.table)

trios_dir <- "~/rocket/git/delocalise/results_trios"
trio_filenames <- list.files(trios_dir, "*.depth", full.names = TRUE)
## all the files available.  which are wanted?
probands_present <- gsub(".depth", "", basename(trio_filenames))
## read the original trios files
trios <- fread("data/1000G_trios.txt", header = FALSE)
## and get the trios
trio_names <- data.frame(
  probands = gsub(".final.cram", "", basename(trios$V1)),
  fathers  = gsub(".final.cram", "", basename(trios$V2)),
  mothers = gsub(".final.cram", "", basename(trios$V3))
)
cat(length(trio_filenames), nrow(trios), nrow(trio_names), "\n")
## reorder the rows of trio_names to match the trio_filenames
m <- match(trio_names$proband, probands_present)
trio_filenames <- trio_filenames[m]
cat(length(trio_filenames), nrow(trios), nrow(trio_names), "\n")

raw <- lapply(trio_filenames, fread)
## extract position information
info <- split(raw[[1]], raw[[1]]$sampleName)[[1]][, 1:5]
colnames(info) <- c("chrom", "start", "end", "ensemblid", "symbol")

getrel <- function(ii) {
  dat <- raw[[ii]]
  trio <- trio_names[ii, ]
  inds <- split(dat, factor(dat$sampleName, levels = trio))
  median_depth <- sapply(inds, function(x) median(x$meanCoverage))
  return(
    sapply(1:3, function(ii) 2*inds[[ii]]$meanCoverage / median_depth[ii])
  )
}

reldepth <- lapply(seq_len(nrow(trio_names)), getrel)
cn <- lapply(reldepth, function(x) 
  matrix(
    cut(x, breaks= c(-0.01, 0.4, 1.5, 2.5, 3.5, 1000),
    labels = c(0:4)),
  ncol = 3)
)
trio_string <- sapply(cn, function(x) apply(x, 1, paste, collapse = "-"))
rownames(trio_string) <- info$symbol
colnames(trio_string) <- trio_names$proband

table(trio_string["NPHP1", ])
table(trio_string["PQLC2", ])
