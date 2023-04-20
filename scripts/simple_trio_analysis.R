
library(data.table)

d <- fread("~/git/delocalise/results/trio1.depth")

# get a vector with the sample names in the order child, father, mother
# that way when we split them up, the output will be in that order
cfm_trio <- c("NA12878", "NA12891","NA12892")
# split the data by sampleName
inds <- split(d, factor(d$sampleName, levels=cfm_trio))
names(inds)  ## this is in child-mother-father order as expected

median_depth <- sapply(inds, function(x) median(x$meanCoverage))

reldepth <- sapply(1:3, function(ii) 2*inds[[ii]]$meanCoverage/median_depth[ii] )
cn <- apply(reldepth, 2,  function(x) cut(x, c(-0.01,0.4,1.5,2.5,3.5,100), labels = c(0:4)))

trio <- apply(cn, 1, paste, collapse="-")
tr_depth <- apply(reldepth, 2, paste, sep=";")

table(trio)

res <- data.frame(inds[[1]][,c(1,2,3,4,5)], trio, tr_depth)
colnames(res) <- c("chrom", "start", "end", "ensemblid", "symbol", "trio", "tr_depth")
rownames(res) <-  res$symbol

res["NPHP1", ]
