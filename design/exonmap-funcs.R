library(exonmap)

xmapConnect()

con <- exonmap:::.xmap.get.con("core")
test <- exonmap:::.xmap.call(NULL,"xmap_AllGenes",con=con)

genes <- test$stable_id

gen <- genes[2]

trxs <- gene.to.transcript(gen)

trx <- trxs[1]

##for(j in 1:length(trxs)){




getAllGenes <- function(subset="core"){
  con <- exonmap:::.xmap.get.con(subset)
  test <- exonmap:::.xmap.call(NULL, "xmap_AllGenes", con=con)
  genes <- test$stable_id
  genes
}

getAllTranscripts <- function(subset="core") {
  genes <- getAllGenes(subset=subset)
  trxs <- gene.to.transcript(genes)
  trxs
}

getAllPredictionTranscripts <- function(subset="core") {

  con <- exonmap:::.xmap.get.con(subset)
  test <- exonmap:::.xmap.call(NULL, "xmap_predictionTranscripts", con=con)
  genes <- test$display_label
  genes <- as.numeric(sub("GENSCAN0+", "",genes, perl=TRUE))
  genes
}
  
getTranscriptInfo <- function(trxid, subset="core"){

if(is.numeric(trxid)){
  subset <- "prediction" }
  
exons <- transcript.to.exon(trxid, subset=subset)

exonvec <- NULL
probesetvec <- NULL

for(j in 1:length(exons)){
   res <- exon.to.probeset(exons[j], subset = subset)
   #print(res)
   if(length(res) != 0){
   exonvec <- c(exonvec,rep(exons[j], length(res)))
   probesetvec <- c(probesetvec, res)
   }
   else{ exonvec <- c(exonvec, exons[j])
         probesetvec <- c(probesetvec, NA)
       }
   }
   if(any(subset == c("core","est"))){
   gen <- transcript.to.gene(trxid)  
   }
   else{
     pad <- 11 - nchar(trxid)
     pad <- paste(rep(0,pad), collapse="")
     trxid <- paste("GENSCAN", pad, trxid, sep="")
     gen <- trxid
     
   }
   genvec <- rep(gen, length(probesetvec))
   tranvec <- rep(trxid, length(probesetvec))

   #print(genvec)
   #print(tranvec)
   #print(exonvec)
   #print(probesetvec)

   out <- data.frame(tranvec, exonvec, probesetvec)
   colnames(out) <- c("transid", "EnsemblExonID", "ProbesetID") 

   out
   
}


getExonLocationsTranscript <- function(trxid, orderByStrand=FALSE, subset="core"){
   exons <- transcript.to.exon(trxs, subset=subset)
   exoninfo <- exon.details(exons, subset=subset)
   start <- exoninfo$seq_region_start
   strand <- exoninfo$seq_region_strand
   result <- data.frame(exons, start, strand)
   colnames(result) <- c("EnsemblExonID", "ExonStart", "Strand")

   if (orderByStrand) {
    if (result$Strand[1] == -1) {
        result <- result[order(result$ExonStart, decreasing = TRUE),]
       }
    }
   
   result
 }


annotateResults <- function(resultsframe) {

       testid <- resultsframe[1,1]


       #print(testid)
	if(length(grep("ENSMUSG",testid))==1){
          genes <- as.character(resultsframe[,1])   
          outframe <- gene.details(as.character(genes)
          cols <- c("display_label", "name", "seq_region_start", "seq_region_end", "biotype")
          outframe[,cols]
          dbres <- data.frame(genes, outframe)
          colnames(dbres) <- c("EnsemblGeneID", "sym", "chr", "TrStart", "TrEnd", "Biotype") 
          #sql <- paste("select DISTINCT EnsemblGeneID, sym, chr, TrStart, TrEnd from probe2ensembl;")
	 }
      if(length(grep("ENSMUST",testid))==1){
          trxs <- as.character(resultsframe[,1])
          outframe <- transcript.details(trxs)
          genes <- gene.to.transcript(trxs)
          cols <- c("display_label", "name", "seq_region_start", "seq_region_end", "biotype")
          outframe <- outframe[,cols]
          dbres <- data.frame(trxs, genes, outframe)
          colnames(dbres) <- c("EnsemblTranscriptID","EnsemblGeneID", "sym", "chr", "TrStart", "TrEnd", "Biotype") 
       #sql <- paste("select DISTINCT EnsemblTranscriptID, EnsemblGeneID, sym, chr, TrStart, TrEnd from probe2ensembl;")
       }
	#print(sql)
 
       #dbres <- dbGetQuery(dbconnobj, sql)
	res <- merge(dbres, resultsframe, by=1)
  #mouseDisconnect(dbconnobj)
  
	 return(res)
}
