
library(GenomicRanges)
a <- GRanges(read.csv("results/NA20752.calls", sep=","))
head(a)

gr <- GRanges(a)

library(karyoploteR)


kp <- plotKaryotype(genome = "hg38", chromosomes=c("autosomal"))

kpPlotRegions(kp, gr[gr$cn==0], col="red")
kpPlotRegions(kp, gr[gr$cn==1], col="magenta")
kpPlotRegions(kp, gr[gr$cn==2], col="black")
kpPlotRegions(kp, gr[gr$cn==3], col="green")
kpPlotRegions(kp, gr[gr$cn==4], col="blue")


kp <- plotKaryotype(genome = "hg38", zoom="chr1:140000000-150000000")

kpPlotRegions(kp, gr[gr$cn==0], col="red")
kpPlotRegions(kp, gr[gr$cn==1], col="magenta")
kpPlotRegions(kp, gr[gr$cn==2], col="black")
kpPlotRegions(kp, gr[gr$cn==3], col="green")
kpPlotRegions(kp, gr[gr$cn==4], col="blue")