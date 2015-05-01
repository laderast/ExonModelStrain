xmap_getAllGenes <- function(){
  #con <- xmapcore:::.xmap.get.con(subset)
  #test <- xmapcore:::.xmap.call(NULL, "xmap_AllGenes", con=con)
  #genes <- test$stable_id
  genes <- all.genes()
  genes
}

xmap_getAllTranscripts <- function() {
  #genes <- getAllGenes(subset=subset)
  trxs <- all.transcripts()
  trxs
}

xmap_getAllPredictionTranscripts <- function() {

  con <- xmapcore:::.xmap.get.con(subset)
  test <- xmapcore:::.xmap.call(NULL, "xmap_predictionTranscripts", con=con)
  genes <- test$display_label
  genes <- as.numeric(sub("GENSCAN0+", "",genes, perl=TRUE))
  genes
}
  
xmap_getTranscriptInfo <- function(trxid, noOverhang, uniqueTargets, remove.na=T, grabProbesetStart=F){

if(is.numeric(trxid)){
  subset <- "prediction" }
  
#exons <- transcript.to.exon(trxid, subset=subset)
exons <- transcript.to.exon(trxid)

exonvec <- NULL
probesetvec <- NULL
startvec <- NULL

  for(j in 1:length(exons)){
   #res <- exon.to.probeset(exons[j], subset = subset)
   res <- exon.to.probeset(exons[j])
   print(res)


   if(length(res) > 0){
	#grab probeset start (first probe position)
	if(grabProbesetStart){
		for(re in res){
		st <- start(probeset.to.hit(as.character(re)))[1]
		#print(st)
		startvec <- c(startvec, st)
		#print(startvec)
		}
	}
   
   
   
   exonvec <- c(exonvec,rep(exons[j], length(res)))
   #print(exonvec)
   probesetvec <- c(probesetvec, res)
   #print(probesetvec)
   }

   #if(!remove.na){
   #   exonvec <- c(exonvec, exons[j])
    #     probesetvec <- c(probesetvec, NA)
    #  }

 }


#print(exonvec)
#print(probesetvec)
#print(startvec)

   #if(any(subset == c("core","est"))){
   #gen <- transcript.to.gene(trxid)  
   #}

   #else{
   #  pad <- 11 - nchar(trxid)
   #  pad <- paste(rep(0,pad), collapse="")
   #  trxid <- paste("GENSCAN", pad, trxid, sep="")
   #  gen <- trxid     
   #}

   #genvec <- rep(gen, length(probesetvec))
   tranvec <- rep(trxid, length(probesetvec))

   #print(genvec)
   #print(tranvec)
   #print(exonvec)
   #print(probesetvec)
   #print(startvec)

if(any(is.null(exonvec), is.null(probesetvec))){
   ##return null data.frame here
  out <- NULL
}

else{
   if(grabProbesetStart){
	out <- data.frame(tranvec, exonvec, probesetvec, startvec)
	colnames(out) <- c("transid", "EnsemblExonID", "ProbesetID", "Start") 
   }
   else{
	out <- data.frame(tranvec, exonvec, probesetvec)
	colnames(out) <- c("transid", "EnsemblExonID", "ProbesetID") 
   
   }
 }

out
   
}


xmap_getExonLocationsTranscript <- function(trId, orderByStrand=FALSE, noOverhang, uniqueTargets){
   
   #exons <- transcript.to.exon(trId, subset=subset)
   exons <- transcript.to.exon(trId)
   exoninfo <- exon.details(exons)

   #exoninfo <- exoninfo[exoninfo$transcript_stable_id == trId,]
   #start <- exoninfo$seq_region_start
   #strt <- start(exoninfo)
   #exons <- exoninfo$stable_id
   #strand <- exoninfo$strand
   result <- data.frame(exoninfo$stable_id, start(exoninfo), exoninfo$strand)
   colnames(result) <- c("EnsemblExonID", "ExonStart", "Strand")

   if (orderByStrand) {
    if (result$Strand[1] == -1) {
        result <- result[order(result$ExonStart, decreasing = TRUE),]
       }
    }
   
   result
 }


xmap_getGeneSym <- function(Id) {

#dbconnobj <- mouseConnect()
  
#if(length(grep("G",Id)==1)){
 #info <- gene.details(Id)
#}


#if(length(grep("T",Id)==1)){
 result <- gene.to.symbol(transcript.to.gene(Id))
#}

#result <- info$display_label

result <- as.character(result)

result 

}




xmap_annotateResults <- 
function (resultsframe) 
{
    testid <- resultsframe[1, 1]
    if (length(grep("G", testid)) == 1) {
        genes <- as.character(resultsframe[, 1])
        outframe <- gene.details(as.character(genes))
        cols <- c("stable_id", "name", "seq_region_start", "seq_region_end", 
            "biotype")
        outframe[, cols]
        dbres <- data.frame(genes, as.data.frame(outframe))
        colnames(dbres) <- c("EnsemblGeneID", "sym", "chr", "TrStart", 
            "TrEnd", "Biotype")
    }
    if (length(grep("T", testid)) == 1) {
        trxs <- as.character(resultsframe[, 1])
        outframe <- transcript.details(trxs)

	  #print(outframe)

	  starts <- start(outframe)
	  ends <- end(outframe)
	  biotype <- outframe$biotype
	  status <- outframe$status
	  chrs <- outframe$space
        trxs <- outframe$stable_id
	 

        genes <- NULL
        syms <- NULL

        #cols <- c("space", "biotype", "status", "start", "end")
       
        outframe <- data.frame(chrs, biotype, status, starts, ends)
	  #class(outframe)
	  #colnames(outframe) <- cols

        for (tr in trxs) {
		print(tr)
            g <- transcript.to.gene(tr)
		ss <- gene.to.symbol(g)
            genes <- c(genes, g)
		syms <- c(syms, ss)

        }

        dbres <- data.frame(trxs, genes, syms, outframe)
        colnames(dbres) <- c("EnsemblTranscriptID", "EnsemblGeneSymbol", 
            "Symbol", "Chromosome", "Biotype", "Status", "TrStart", 
            "TrEnd")
    }
    res <- merge(dbres, resultsframe, by = 1)
    return(res)
}