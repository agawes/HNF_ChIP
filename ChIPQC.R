#module load R/3.2.5

library(ChIPQC)
samples = read.csv("HNF_ChIP.sample_sheet.csv")

exampleExp = ChIPQC(samples,annotaiton="hg19")

ChIPQCreport(exampleExp)

pdf("Reads_in_GenomicIntervals.pdf")
plotRegi(exampleExp,facetBy=c("Tissue","Factor"))
dev.off()