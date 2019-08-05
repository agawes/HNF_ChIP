# module load R/3.3.1
# R

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ReactomePA)
require(lattice)

peak_dir="/well/mccarthy/users/agata/StemBANCC/HNF4_1_ChIP_Singapore/aquas_flt_narrowPeak/"
peak_files=list.files(path=peak_dir, pattern=".narrowPeak.gz$",recursive=TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

for (f in peak_files){
	name = gsub(".filt.narrowPeak.gz","",f)
	print(name)
	peaks <- readPeakFile(gzfile(f), header=F)

	pdf(paste0(name,".covplot.pdf"))
	print(covplot(peaks, weightCol="V5"))
	dev.off()

	### plot peak distances from their nearest TSS:
	tagMatrix <- getTagMatrix(peaks, windows=promoter)
	pdf(paste0(name,".TSS_distance.pdf"))
	print(plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency"))
	dev.off()

	## peak to annotation piechart
	peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
	pdf(paste0(name,".peakAnno_pie.pdf"))
	print(plotAnnoPie(peakAnno))
	dev.off()

	### top 10 enriched terms from Reactome
	genes <- seq2gene(peaks, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
	pathways <- enrichPathway(genes)

	pdf(paste0(name,".Reactome.pdf"))
	print(barplot(pathways, showCategory=10))
	dev.off()
}

peak_files_list=list()
for (f in peak_files){
	name = gsub(".filt.narrowPeak.gz","",f)
	print(name)
	peak_files_list[[name]] <- readPeakFile(gzfile(f), header=F)
}
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peak_files_list, getTagMatrix, windows=promoter)

pdf("avgProfTSS.allSamples.pdf",height=12)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()


pdf("peakAnno.allSamples.pdf")
peakAnnoList <- lapply(peak_files_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
dev.off()

pdf("DistToTSS.allSamples.pdf")
plotDistToTSS(peakAnnoList)
dev.off()

library(clusterProfiler)
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

 pdf("KEGG_enrichment.allSamples.pdf", width=20, height=10)
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()



###### repeat all for IDR peaks
peak_dir="/well/mccarthy/users/agata/StemBANCC/HNF4_1_ChIP_Singapore/aquas_IDR_narrowPeak/"
peak_files=list.files(path=peak_dir, pattern=".narrowPeak.gz$",recursive=TRUE)

### filter to files with >100 IDR peaks; otherwise throws errors...
peak_files=peak_files[c(2, 7, 9, 10, 11, 13:15)]

for (f in peak_files){
	name = gsub("_rep1-pr.IDR0.05.filt.narrowPeak.gz","",f)
	print(name)
	peaks <- readPeakFile(gzfile(f), header=F)

	pdf(paste0(name,".covplot.pdf"))
	print(covplot(peaks, weightCol="V5"))
	dev.off()

	### plot peak distances from their nearest TSS:
	tagMatrix <- getTagMatrix(peaks, windows=promoter)
	pdf(paste0(name,".TSS_distance.pdf"))
	print(plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency"))
	dev.off()

	## peak to annotation piechart
	peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
	pdf(paste0(name,".peakAnno_pie.pdf"))
	print(plotAnnoPie(peakAnno))
	dev.off()

	### top 10 enriched terms from Reactome
	genes <- seq2gene(peaks, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
	pathways <- enrichPathway(genes)

	pdf(paste0(name,".Reactome.pdf"))
	print(barplot(pathways, showCategory=10))
	dev.off()
}

peak_files_list=list()

for (f in peak_files){
	name = gsub("_rep1-pr.IDR0.05.filt.narrowPeak.gz","",f)
	print(name)
	peak_files_list[[name]] <- readPeakFile(gzfile(f), header=F)
}
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peak_files_list, getTagMatrix, windows=promoter)

pdf("avgProfTSS.allSamples.pdf")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()


pdf("peakAnno.allSamples.pdf")
peakAnnoList <- lapply(peak_files_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
dev.off()

pdf("DistToTSS.allSamples.pdf")
plotDistToTSS(peakAnnoList)
dev.off()

library(clusterProfiler)
pdf("KEGG_enrichment.allSamples.pdf")
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()


## get % of reads in <1kb Promoter
data.frame(unlist(lapply(peakAnnoList, function(x) x@annoStat[3,2])))
