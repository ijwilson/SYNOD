# analyse singletons, one individual per file
library(data.table)

singletons_dir <- "~/rocket/git/delocalise/results_singletons"

singleton_filenames <- list.files(singletons_dir, "*.depth", full.names = TRUE)

probands <- gsub(".depth", "", basename(singleton_filenames))
# check that the results we have match the original set of files to analyse
singletons <- fread("data/1000G_singletons.txt", header = FALSE) 
expected_probands <- gsub(".final.cram", "", basename(singletons$V1))

table(expected_probands %in% probands)  ## are all present
singletons[!(expected_probands %in% probands)]
raw <- lapply(singleton_filenames, fread)  ## read raw data
## get the relative depth for each individual
rel <- sapply(raw, function(x) {
  m <- median(x$meanCoverage); return(2*x$meanCoverage/m);})

samples <-  sapply(raw, function(x) x$sampleName[1])
colnames(rel) <- samples
cn <- apply(rel,2, cut, c(-1, 0.4, 1.5, 2.5, 3.5, 1000), labels=0:4)

pos <- raw[[1]][, 1:5]
colnames(pos) <- c("chrom", "start", "end", "ensemblid", "symbol")

results <- cbind(pos, cn)

table(cn[results$symbol == "NPHP1", ])
table(cn[results$symbol == "PQLC2", ])

ebi_dir <- "http://ftp.1000genomes.ebi.ac.uk/"
ftp_dir <- "vol1/ftp/data_collections/1000G_2504_high_coverage"

ped_address <- file.path(ebi_dir, ftp_dir,
                               "20130606_g1k_3202_samples_ped_population.txt")

ped <- read.table(ped_address, header = TRUE)
colnames(ped)[1:4] <-
              c("Family_ID", "Individual_ID", "Paternal_ID", "Maternal_ID")
