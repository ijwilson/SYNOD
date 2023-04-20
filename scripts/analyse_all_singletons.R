# analyse singletons, one individual per file
library(data.table)
singleton_filenames <- list.files("results_singletons/", "*.depth", full.names = TRUE)

probands <- gsub(".depth", "", basename(singleton_filenames))
# check that the results we have match the original set of files to analyse
singletons <- fread("data/1000G_singletons.txt", header = FALSE) 
expected_probands <- gsub(".final.cram", "", basename(singletons$V1))

table(expected_probands %in% probands)  ## are all present

raw <- lapply(files, fread)  ## read raw data
## get the relative depth for each individual
rel <- sapply(raw, function(x) {
  m <- median(x$meanCoverage); return(2*x$meanCoverage/m);})

samples <-  sapply(raw, function(x) x$sampleName[1])
colnames(rel) <- samples
cn <- apply(rel,2, cut, c(-1, 0.4, 1.5, 2.5, 3.5, 1000), labels=0:4)

pos <- raw[[1]][, 1:5]
colnames(pos) <- c("chrom", "start", "end", "ensemblid", "symbol")

results <- cbind(pos, cn)

table(cn[results$symbol=="NPHP1", ])
table(cn[results$symbol=="PQLC2", ])



