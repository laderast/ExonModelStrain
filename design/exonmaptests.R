out2 <- out[-indnas,]
dim(out)
dim(out2)
out2[1:5,]
ls()
candidates <- merge(unbalanced, out2, by.x=1, by.y=2)
dim(candidates)
write.table(candidates, "candidate-exons-unbalanced.txt", quote=F, sep = "\t", row.name=F)
savehistory('ranksortout.R')
PlotExonResultsCI
genes[1:5,]
genes[1:5]
dim(genes)
length(genes)
trxs <- gene.to.transcript(gen)
trx
trxs
exons <- transcript.to.exon(trxs[1])
exons
exon.details(exons)
getExonLocationsTranscript <- function(trxid, orderByStrand=FALSE){ exons <- transcript.to.exon(trxs) exoninfo <- exon.details(exons) start <- exoninfo$seq_region_start strand <- exoninfo$seq_region_strand result <- data.frame(exons, start, strand) colnames(result) <- c("EnsemblExonID", "ExonStart", "Strand") if (orderByStrand) { if (result$Strand[1] == -1) { result <- result[order(result$ExonStart, decreasing = TRUE),] } } result }
getExonLocationsTranscript(trxs[1])
> + + getExonLocationsTranscript(trxs[1])
> }+ +
{} getExonLocationsTranscript(trxs[1]) > + + getExonLocationsTranscript(trxs[1])
> + + getExonLocationsTranscript(trxs[1]) > + + )
> + + getExonLocationsTranscript(trxs[1]) > + + getExonLocationsTranscript(trxs[1])
library(exonmap)
exonmap:::.xmap.get.db.name()
xmapConnect()
exonmap:::.xmap.get.db.name()
exonmap:::.xmap.internals
