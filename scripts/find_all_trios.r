# find trios

library(data.table)
#'
#' The files for our analyses can be foundin the directory
#' http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/
#' on the 1000 genomes data site.
#'
#'
#' The README for this directory is below.  Information about reuse and a
#' pipeline description is also in this directory
#'
#'
#' 25th May 2020
#'
#' 2504 phase three panel plus additional 698 related samples, totaling 3202 samples # nolint
#' ================================================================================= # nolint

#' As described in the READMEs in these directories, NYGC released 30x coverage
#' seqeunce data on the 2504 samples in the 1000 Genomes Project's phase three
#' panel. This data was released in 2019 and is listed in our index files and
#' available in ENA as study ERP114329.

#' Following release of the 2504 unrelated samples in the phase three panel,
#' NYGC have, in addition, sequenced a further 698 samples from the collections
#' created #' for the 1000 Genomes Project. This data is also at 30x coverage
#' and was generated in the same way at NYGC. The set of 698 samples is composed
#' of samples that are related to samples in the main 2504 panel, mainly
#' completing trios (sets of parents and their children). The additional 698
#' have been added here and are available in ENA under study ERP120144. The data
#' for the additional 698 related samples was released in 2020.

#' The total number of samples sequenced to 30x coverage by NYGC is 3202.

#' For further details, including data reuse and accurate citation and
#' attribution of these data sets, please see the accompanying documents in
#' these directories. If you have any questions, please
#' contact info@1000genomes.org.

#+ run
ebi_dir <- "http://ftp.1000genomes.ebi.ac.uk/"
ftp_dir <- "vol1/ftp/data_collections/1000G_2504_high_coverage"

highcoverage_index <- file.path(ebi_dir, ftp_dir,
                                "1000G_2504_high_coverage.sequence.index")
related_index_address <- file.path(ebi_dir, ftp_dir,
                               "1000G_698_related_high_coverage.sequence.index")
ped_address <- file.path(ebi_dir, ftp_dir,
                               "20130606_g1k_3202_samples_ped_population.txt")

ped <- read.table(ped_address, header = TRUE)
colnames(ped)[1:4] <-
              c("Family_ID", "Individual_ID", "Paternal_ID", "Maternal_ID")

samples <- read.table(highcoverage_index, sep = "\t")
relatives <- read.table(related_index_address, sep = "\t")
relatives$V12 <- NULL
col_names <- scan(highcoverage_index, skip = 23, n = 22, what = character())
col_names[1] <- "file_path"
colnames(samples) <- colnames(relatives) <- col_names
samples$set <- "base"
relatives$set <- "related"
all <- rbind(samples, relatives)
# tidy columns
all$STUDY_NAME <- all$PAIRED_FASTQ <- all$ANALYSIS_GROUP <- all$RUN_NAME <-
all$CENTER_NAME <- all$SUBMISSION_ID <- all$INSTRUMENT_PLATFORM <-
all$INSTRUMENT_MODEL <- all$INSERT_SIZE <- all$MD5SUM <- all$SUBMISSION_DATE <-
all$RUN_ID <- NULL
file_locations <- all$file_path
names(file_locations) <- all$SAMPLE_NAME

trios <- ped[
    ped$Individual_ID %in% all$SAMPLE_NAME &
    ped$Paternal_ID %in% all$SAMPLE_NAME &
    ped$Maternal_ID %in% all$SAMPLE_NAME, ]

head(trios)

trio_files <- cbind(file_locations[trios$Individual_ID],
                    file_locations[trios$Paternal_ID],
                    file_locations[trios$Maternal_ID])

sum(is.na(trio_files))

nrow(trio_files)

write.table(trio_files, file = "1000G_trios.txt", row.names = FALSE,
            col.names = FALSE, sep = " ", quote = FALSE)

singletons <- ped[!ped$Individual_ID %in% trios$Individual_ID, ]
nrow(singletons)
singletons <- singletons[singletons$Individual_ID %in% all$SAMPLE_NAME, ]

head(singletons)
nrow(singletons)

cat(file_locations[singletons$Individual_ID],
     file = "100G_singletons.txt", sep = "\n")
