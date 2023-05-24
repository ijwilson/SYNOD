# analyse trios
library(data.table) 
trio_filenames <- list.files("results_trios", "*.depth", full.names = TRUE)

probands <- gsub(".depth", "", basename(trio_filenames))
## read the original trios files
trios <- fread("data/1000G_trios.txt", header = FALSE) 

trio_names <- data.frame(
  probands = gsub(".final.cram", "", basename(trios$V1)),
  fathers  = gsub(".final.cram", "", basename(trios$V2)),
  mothers =gsub(".final.cram", "", basename(trios$V3))
)

## reorder the rows of trio_names to match the trio_filenames
m <- match(probands, trio_names$proband)
trio_names <- trio_names[m, ]

table(trio_names$probands %in% probands)  ## so all present

raw <- lapply(trio_filenames, fread)
## extract position information
info <- split(raw[[1]], raw[[1]]$sampleName)[[1]][,1:5]
colnames(info) <- c("chrom", "start", "end", "ensemblid","symbol")

getrel <- function(ii) {
  dat <- raw[[ii]]
  trio <- trio_names[ii, ]
  inds <- split(dat, factor(dat$sampleName, levels = trio))
  median_depth <- sapply(inds, function(x) median(x$meanCoverage))
  return(
    sapply(1:3, function(ii) 2*inds[[ii]]$meanCoverage/median_depth[ii])
  )
}


reldepth <- lapply(seq_len(nrow(trio_names)), getrel)
cn <- lapply(reldepth, function(x) 
  matrix(cut(x, breaks= c(-0.01,0.4,1.5,2.5,3.5,1000), labels = c(0:4)), ncol = 3)
)
trio_string <- sapply(cn, function(x) apply(x, 1, paste, collapse = "-"))
rownames(trio_string) <- info$symbol
colnames(trio_string) <- trio_names$proband

table(trio_string["NPHP1",])
table(trio_string["PQLC2",])
