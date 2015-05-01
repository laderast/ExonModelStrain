#Exon <- NULL
#Expression <- NULL
#Strain <- NULL
#Subject <- NULL

dbenv <- new.env(hash=TRUE, parent=emptyenv())

dbInit <- function(){
   dbenv[["connected"]] <- FALSE
   dbenv[["dbPackage"]] <- NULL
   dbenv[["analysisType"]] <- NULL
   dbenv[["getExonLocationsTranscript"]] <- NULL
   dbenv[["getAllGenes"]] <- NULL
   dbenv[["getAllTranscripts"]] <- NULL
   dbenv[["getGeneSym"]] <- NULL
   dbenv[["annotateResults"]] <- NULL
   dbenv[["getTranscriptInfo"]] <- NULL
}


#global functions that are mapped to dbSpecific functions
#by mapConnect()
#note that mapDisconnect() makes these NULL again
getExonLocationsTranscript <- 
	function(trId, orderByStrand=FALSE, subset="core", noOverhang, uniqueTargets) {

	.getExonLocationsTranscript <- dbenv[["getExonLocationsTranscript"]]
	exons <- .getExonLocationsTranscript(trId, orderByStrand=FALSE, noOverhang, uniqueTargets)
	exons
	}

getAllGenes <- function(subset="core") {

	.getAllGenes <- dbenv[["getAllGenes"]]
	genes <- .getAllGenes(subset)
	genes
}
	
getAllTranscripts <- 

	function(subset="core") {
	.getAllTranscripts <- dbenv[["getAllTranscripts"]]
	trxs <- .getAllTranscripts(subset="core")
	trxs
	}
	
getGeneSym <- function(Id) {

	.getGeneSym <- dbenv[["getGeneSym"]]
	sym <- .getGeneSym(Id)
	sym
}

annotateResults <- function(resultsframe) {

	.annotateResults <- dbenv[["annotateResults"]]
	results <- .annotateResults(resultsframe)
	results
}

getTranscriptInfo <- 
	function(trxid, subset="core", noOverhang, uniqueTargets, remove.na=T, grabProbesetStart=FALSE){

	.getTranscriptInfo <- dbenv[["getTranscriptInfo"]]
	result <- .getTranscriptInfo(trxid, noOverhang, uniqueTargets)
	result
	}




mapConnect <- function(dbPackage, dbName="mouse") {

  if(missing(dbPackage)){
    stop("Please specify dbPackage as either 'mouseexonensembl.db' or
          'xmapcore'")
  }
  
  dbenv[["dbPackage"]] <- dbPackage
  
  if(dbPackage == "mouseexonensembl.db"){
       library(mouseexonensembl.db)
       dbenv[["connected"]] <- TRUE
       #dbenv[["dbPackage"]] <- "mouseexonensembl.db"

       #map appropriate db functions for package
       dbenv[["getExonLocationsTranscript"]] <- exonensembl_getExonLocationsTranscript
       dbenv[["getAllGenes"]] <- exonensembl_getAllGenes
       dbenv[["getAllTranscripts"]] <- exonensembl_getAllTranscripts
       dbenv[["getGeneSym"]] <- exonensembl_getGeneSym
       dbenv[["annotateResults"]] <- exonensembl_annotateResults
       dbenv[["getTranscriptInfo"]] <- exonensembl_getTranscriptInfo
   }

   if(dbPackage == "xmapcore") {
       library(xmapcore)
       if(!missing(dbName)){
         xmap.connect(dbName)}

       else{xmap.connect()}
       
       dbenv[["getExonLocationsTranscript"]] <- xmap_getExonLocationsTranscript
       dbenv[["getAllGenes"]] <- xmap_getAllGenes
       dbenv[["getAllTranscripts"]] <- xmap_getAllTranscripts
       dbenv[["getGeneSym"]] <- xmap_getGeneSym
       dbenv[["annotateResults"]] <- xmap_annotateResults
       dbenv[["getTranscriptInfo"]] <- xmap_getTranscriptInfo
       dbenv[["dbPkg"]] <- "xmapcore"
       dbenv[["connected"]] <- TRUE
   }

}


mapDisconnect <- function() {
   if(dbenv[["connected"]] == FALSE){
     stop("Not connected to any database")
   }
   if(dbenv[["dbPackage"]] == "exonmap") {
     xmapDisconnect()
     dbenv[["connected"]] <- FALSE
     dbenv[["dbPackage"]] <- NULL
   }
   if(dbenv[["dbPackage"]] == "mouseexonensembl.db"){
     dbenv[["connected"]] <- FALSE
     dbenv[["dbPackage"]] <- NULL
   }

   dbenv[["getExonLocationsTranscript"]] <- NULL
   dbenv[["getAllGenes"]] <- NULL
   dbenv[["getAllTranscripts"]] <- NULL
   dbenv[["getGeneSym"]] <- NULL
   dbenv[["annotateResults"]] <- NULL
   dbenv[["getTranscriptInfo"]] <- NULL
}



fetchdata <- function(Id, ExpSet, uniqueTargets=TRUE,
                      noOverhang=FALSE, grabProbesetStart=FALSE) {
   require(Biobase)
   if(class(ExpSet)!="ExpressionSet"){
	stop("Expression set required")
    }

   #library(dBPkg,character.only=TRUE)

   samplelength <- length(sampleNames(ExpSet))
   subject <- c(1:samplelength) 
   strain <- ExpSet$Strain
   transdat <- NULL
   info <- NULL
   ind <- NULL

  if(length(grep("T",Id))==1){
     info <- suppressWarnings(getTranscriptInfo(Id, uniqueTargets=uniqueTargets,
	noOverhang=noOverhang, grabProbesetStart=grabProbesetStart))
  }
  if(length(grep("G",Id))==1){
     info <- suppressWarnings(getGeneInfo(Id, uniqueTargets=uniqueTargets, noOverhang=noOverhang))
  }

if(!is.null(info)){
    ind <- intersect(info$ProbesetID, featureNames(ExpSet))
    #print(ind)

    info <- info[which(info$ProbesetID %in% ind),]
    #print(info)
    
    expdat <- exprs(ExpSet)[ind,]

    if(length(ind) ==1){
       expdat <- as.data.frame(expdat)
       expdat <- t(expdat)
       rownames(expdat) <- ind
    }
}

#print(expdat)

if(length(ind) > 0) {

      transvec <- rep(Id, samplelength)
      rowvec <- c(1:nrow(info))
      
        for (rownum in rowvec){
           
           #for each row in dat do the following:
           #extract exon id and pad a vector
           exonid <- info[rownum,"EnsemblExonID"]
           exonvec <- rep(exonid, samplelength)
           #extract probeset id and pad a vector
           psid <- as.character(info[rownum,"ProbesetID"])
           psvec <- rep(psid, samplelength)
	   
	   if(grabProbesetStart){
	   strt <- info[rownum, "Start"]
	   stvec <- rep(strt, samplelength)
	   }
   
	#extract expression values for row
           rowdata <- as.vector(expdat[psid,])
           rowdata <- t(rowdata)
	     if(!is.na(mean(rowdata))){ 
               #put together vectors into a submatrix
	       
	       if(grabProbesetStart){
			rowdatum <- data.frame(transvec, psvec, subject, as.numeric(strain), 
				exonvec, stvec, as.numeric(rowdata))
			colnames(rowdatum) <- c("transid", "probesetid", "Subject", "Strain", 
				"Exon", "Start", "Expression")
			}
	       
	       else{
			rowdatum <- data.frame(transvec,psvec,subject, as.numeric(strain), 
				exonvec, as.numeric(rowdata))
			colnames(rowdatum) <- c("transid","probesetid","Subject", "Strain", 
				"Exon", "Expression")
			}
		
		#append submatrix onto transdat
               transdat <- rbind(transdat,rowdatum)

		} 
		   
         }

}

else {transdat <- NULL}

transdat

}



runModel <- function(transdat, ExpSet, analysisUnit="exon") {

if(class(transdat)!="data.frame"){ 
      stop("Data Frame required")
    }
if(class(ExpSet)!="ExpressionSet"){
	stop("Expression set required")
    }

 #require(Biobase)

Probeset <- factor(transdat$probesetid)
Exon <- factor(transdat$Exon)
Expression <- transdat$Expression
Subject <- factor(transdat$Subject)
Strain <- factor(transdat$Strain)

 #attach(transdat)
 samplenumber <- length(sampleNames(ExpSet))

 ex <- as.character(Exon)
 exontab <- table(ex)
 exonnam <- levels(ex)
 exontab <- as.data.frame(exontab)$Freq
 names(exontab) <- exonnam
 exontab <- exontab/samplenumber
 
 exonlength <- length(exontab)

 #print(exontab)
 status <- NULL

#print(Strain)
#print(Expression)

 means <- tapply(Expression, Strain, mean, simplify=TRUE)

 numPrsetsSingleExon <- NULL

 if(length(levels(Probeset)) ==1 || length(exontab) == 1 )

 {
  status <- "single"
  numPrsetsSingleExon <- nrow(transdat)
  lm1 <- anova(lm(Expression~Strain+Subject%in%Strain))

 } 

 else{
    #run linear model on transdat
    
    if(analysisUnit=="exon"){
    lm1<-anova(lm(Expression~Strain+Exon+Subject%in%Strain + Exon:Strain))
    }
    
    if(analysisUnit=="probeset"){
    lm1<-anova(lm(Expression~Strain+Probeset+Subject%in%Strain + Probeset:Strain))    
    }
 
    if(max(exontab)==1){
       status <- "perexon"
      }
    else{
       status <- "multiple"
     }
   }

 if(status == "single"){
       tab <- print(lm1)
       fields <- c("Strain")
       pvals <- tab[fields,c("Pr(>F)")]
       names(pvals) <- fields
   }

 else {
	#extract anova table from results
      tab <- print(lm1)
      #pvalues of interest
      if(analysisUnit=="exon"){
      fields <- c("Strain", "Exon", "Strain:Exon")}
      
      if(analysisUnit=="probeset"){
      fields <- c("Strain", "Probeset", "Strain:Probeset")
      }
      
      #extract pvalues
      pvals <- tab[fields,c("Pr(>F)")]
      names(pvals) <- fields
   }

 #detach(transdat)
 return(list(pvals=pvals, status=status,
             means=means,
             exonlength=exonlength, 
	numPrsetsSingleExon=numPrsetsSingleExon))

}



findMissingExons <- function(transdat) {

	if(class(transdat)!="data.frame"){stop("Input is not data frame")}

      Id <- as.character(transdat$transid[1])

      if(length(grep("T", Id))==1){ exons <- getExonLocationsTranscript(trId=Id) }
      else{exons <- getExonLocationsGene(geneId=Id)}

      ExInDat <- as.character(unique(transdat$Exon))

	trueexons <- as.character(exons$EnsemblExonID)

	res <- setdiff(exons$Ensembl, ExInDat)

	if (length(res) == 0) {res <- NULL}

	return(res)

}

getProbesetDeltas <- function(transdat, uniqueTargets=T, noOverhang=F){
	if(class(transdat)!="data.frame"){stop("Input is not data frame")}

	B6data <- transdat[transdat$Strain==1,]
	D2data <- transdat[transdat$Strain==2,]
	
	attach(B6data)
      #tapply to get Exon means of Expression values
      B6exonmeans <- tapply(Expression, probesetid, mean)
      detach(B6data)
      
	attach(D2data)
      #tapply to get Exon means of Expression values
      D2exonmeans <- tapply(Expression, probesetid, mean)
      detach(D2data)

	deltas <- B6exonmeans - D2exonmeans

	exonlength <- length(deltas)
	#trueexonlength <- nrow(exons)

      #sort deltas
      deltas <- sort(abs(deltas), decreasing=TRUE)

      #extract max value from sorted deltas
      maxexon <- deltas[1]
      #extract exon id with max value
      maxexonname <- names(deltas[1])

	position <- grep(maxexonname, names(deltas))
      
      positionflag <- NA

      #flag 1 if max exon is at beginning of transcript
      if(position==1){ positionflag <-1 }

      #flag 2 if max exon is at end of transcript
      if(position==length(deltas)){ positionflag <-2 }

      #flag 0 otherwise
      if(position > 1 && position < length(deltas)) { positionflag <- 0 }


	missingexonflag <- 0

	return(list(maxexonname=maxexonname, maxexon = maxexon, exonlength=length(deltas), 
           trueexonlength = length(deltas), missingexonflag=missingexonflag ,position = position, positionflag = positionflag, 
	     strand = NA, deltas = deltas))

}


getExonDeltas <- function(transdat, uniqueTargets=T, noOverhang=F) {

	if(class(transdat)!="data.frame"){stop("Input is not data frame")}

      Id <- as.character(transdat$transid[1])

      if(length(grep("T", Id))==1){ exons <- getExonLocationsTranscript(trId=Id,uniqueTargets=uniqueTargets, noOverhang=noOverhang) }
      else{exons <- getExonLocationsGene(geneId=Id, uniqueTargets=uniqueTargets, noOverhang=noOverhang)}

	#if exon locations reversed (last position first)
	#strand is -1
	strand <- exons$Strand[1]


      #attach(transdat)
      #aggregate all B6 data from transdat
      B6data <- transdat[transdat$Strain==1,]
      #aggregate all D2 data from transdat
      D2data <- transdat[transdat$Strain==2,]

      #detach(transdat)

      attach(B6data)
      #tapply to get Exon means of Expression values
      B6exonmeans <- tapply(Expression, as.character(Exon), mean)
      detach(B6data)

      attach(D2data)
      #tapply to get Exon means of Expression values
      D2exonmeans <- tapply(Expression, as.character(Exon), mean)
      detach(D2data)

      #calculate deltas
      deltas <- B6exonmeans - D2exonmeans

	exonlength <- length(deltas)
	trueexonlength <- nrow(exons)

      #sort deltas
      deltas <- sort(abs(deltas), decreasing=TRUE)

      #extract max value from sorted deltas
      maxexon <- deltas[1]
      #extract exon id with max value
      maxexonname <- names(deltas[1])

	position <- grep(maxexonname, exons$EnsemblExonID)
      
      positionflag <- NA

      #flag 1 if max exon is at beginning of transcript
      if(position==1){ positionflag <-1 }

      #flag 2 if max exon is at end of transcript
      if(position==trueexonlength){ positionflag <-2 }

      #flag 0 otherwise
      if(position > 1 && position < trueexonlength) { positionflag <- 0 }

	if(exonlength!=trueexonlength) {
	    missingexonflag <- 1
	}

	else{ missingexonflag <- 0}

	return(list(maxexonname=maxexonname, maxexon = maxexon, exonlength=exonlength, 
           trueexonlength = trueexonlength, missingexonflag=missingexonflag ,position = position, positionflag = positionflag, 
	     strand = strand, deltas = deltas))

}







####


RunExonModelWorkflow <- function(ExpSet,idlist = NULL,analysisType="transcript", analysisUnit="exon",
                                 uniqueTargets=TRUE,noOverhang=FALSE){

options <- list(uniqueTargets=uniqueTargets, noOverhang=noOverhang)


if(dbenv[["connected"]]==FALSE){
  stop("Please connect to a database using mapConnect()")
}

dbenv[["analysisType"]] <- analysisType

#if(is.null(dBPackage)){ stop("Database Package needs to be specified") }
#else{library(dBPackage, character.only=TRUE)}

if(class(ExpSet)!="ExpressionSet"){stop("input is not ExpressionSet")}

if(!(analysisType %in% c("transcript", "gene"))){ 
   stop("analysis type not valid: should be either 'transcript' or 'gene'")
}

if(is.null(idlist)){

   if(analysisType=="transcript"){
   idlist <- getAllTranscripts()
   }

   if(analysisType=="gene"){
   idlist <- getAllGenes()
   }
}

if(class(idlist)!="character") {
   stop("idlist not valid - must be character vector")
   }

idlist <- unique(idlist)

singles <- NULL
multi <- NULL
notrun <- NULL
i <- 1

for(tt in idlist) {
    tess <- paste("running", i, tt, sep=" ")
    print(tess)
    ##trdat <- fetchdata(tt, ExpSet, dBPkg=dBPackage, uniqueTargets=uniqueTargets,noOverhang=noOverhang)
    trdat <- fetchdata(tt, ExpSet)

 if(!is.null(trdat)){

    res <- runModel(trdat, ExpSet, analysisUnit=analysisUnit)
    if(analysisUnit=="exon"){
    deltres <- getExonDeltas(trdat, uniqueTargets=uniqueTargets, noOverhang=noOverhang)
	}
	
    if(analysisUnit=="probeset"){
    deltres <- getProbesetDeltas(trdat, uniqueTargets=uniqueTargets, noOverhang=noOverhang)
    }
    
    
    if(res$status=="single"){
        singles <- rbind(singles,data.frame(tt,t(res$pvals),t(res$means), res$numPrsetsSingleExon,deltres$maxexonname, deltres$maxexon, 
           deltres$exonlength, deltres$trueexonlength, deltres$missingexonflag, deltres$strand))
        
	  singlepvalsnames <- res$pvalsnames
    
    }

    else{
        if(res$status=="multiple"){ modelflag <- 1}
        else{ modelflag <- 0}

        multipvalsnames <- names(res$pvals)

        multi <- rbind(multi,data.frame(tt,t(res$pvals),t(res$means),modelflag,deltres$maxexonname, 
           deltres$maxexon, deltres$position, deltres$positionflag, deltres$exonlength, 
           deltres$trueexonlength, deltres$missingexonflag, deltres$strand))
    
	}
  }
    
  else{notrun <- c(notrun, tt)}
  i <- i +1

}

if(!is.null(singles)){
colnames(singles) <- c("ID", "pStrain", "Strain1mean", "Strain2mean", "numprobesets", "maxexon", "maxexondelta",
    "numexonsmapped","trueexonnum", "missingexonflag","strand")
rownames(singles) <- NULL
}

if(!is.null(multi)){
colnames(multi) <- c("ID", "pStrain", "pExon", "pExonStrain", "Strain1mean", "Strain2mean", "multprobeflag", "maxexon", "maxexondelta",
    "position", "positionflag", "numexonsmapped","trueexonnum", "missingexonflag","strand") 
rownames(multi) <- NULL
}

return(list(multi=multi,singles=singles, notrun=notrun, options=options))

}


RunQVals <- function(resultframe) {

   require(qvalue)

   if(class(resultframe)!= "data.frame" ) {
	stop("input is not a data.frame")
   }

   pnames <- c("pStrain", "pExon", "pExonStrain")

   pnames <- intersect(colnames(resultframe), pnames)

   if(length(pnames)==0) { stop("input is not a resultframe")}

   ind <- which(colnames(resultframe) %in% pnames)

   #print(pnames)

   pvals <- resultframe[,pnames]
   #print(pvals[1:5,])
   qvals <- NULL


   for(p in pnames) {
        qv <- NULL
        try(qv <- qvalue(pvals[,p]),silent=TRUE)
	  #suppressWarnings(qv <- qv$qvalues)
          #print(qv)
          #print(length(qv))
        
	  if(length(qv) == 1){qv <- rep(NA, nrow(pvals))
                      qv <- as.numeric(qv)
		print(paste("qvalues could not be calculated for ", p, sep = "")) 
		}
          else{qv <- qv$qvalues}
          
          #print(qv)    
	  qvals <- cbind(qvals, qv)
	  }

   qnames <- sub('p', 'q', pnames)
   colnames(qvals) <- qnames

   qvals
   out <- cbind(resultframe[,1],pvals, qvals,resultframe[,-c(1,ind)])

   out

}



MissingExons <- function(id, ExpSet) {

   if(class(id)!="character"){ stop("input is not id")}

   #if(is.null(dBPackage)){ stop("Database Package needs to be specified") }

   #else{library(dBPackage, character.only=TRUE)}

   if(dbenv[["connected"]]==FALSE){
     stop("Please connect to a database using mapConnect()")
    }
   
   if(class(ExpSet)!="ExpressionSet"){stop("input is not ExpressionSet")}

   trdat <- fetchdata(id, ExpSet, dBPkg=dBPackage)
   res <- findMissingExons(trdat)

   return(res)

}



getExonPeaks <- function(transdat, threshold=3, uniqueTargets=T, noOverhang=F) {

	if(class(transdat)!="data.frame"){stop("Input is not data frame")}

      Id <- as.character(transdat$transid[1])

      if(length(grep("T", Id))==1){ exons <- getExonLocationsTranscript(trId=Id, uniqueTargets=uniqueTargets, noOverhang=noOverhang) }
      else{exons <- getExonLocationsGene(trId=Id, uniqueTargets=uniqueTargets, noOverhang=noOverhang)}

	#if exon locations reversed (last position first)
	#strand is -1
	if(rownames(exons)[1] == 1) {strand <- 1}
	else {strand <- -1}

      #attach(transdat)
      #aggregate all B6 data from transdat
      B6data <- transdat[transdat$Strain==1,]
      #aggregate all D2 data from transdat
      #D2data <- transdat[transdat$Strain==2,]
      exons <- exons[which(exons$EnsemblExonID %in% Exon),]
	exonorder <- as.character(exons$EnsemblExonID)
	#print(exonorder)
      #detach(transdat)

      attach(B6data)
      #tapply to get Exon means of Expression values
      B6exonmeans <- tapply(Expression, as.character(Exon), mean)
      detach(B6data)

	B6exonmeans <- B6exonmeans[exonorder]

	b6len <- length(B6exonmeans)

	if(b6len > 4){ 	

      deriv <- B6exonmeans[1:(b6len-1)] - B6exonmeans[2:b6len]

	derivlen <- b6len-1

	deriv2 <- deriv[1:(derivlen-1)] - deriv[2:derivlen]

	sign <- NULL
	oldsign <-NULL
	peaklist <- NULL

      for(i in 1:length(deriv)){
		sign <- deriv[i]
		#pick peak if sign changes from positive to negative
		if(!is.null(sign) & !is.null(oldsign)){
				
		if(sign <= 0 & oldsign > 0){

		  if(deriv2[i-1] > threshold){
		      ind <- i
		      peaklist <- c(peaklist, B6exonmeans[ind])
                
     		    }
		  }
            }
		oldsign <- sign

	    }
	}
	else {peaklist <- NULL
		}

peakpositions <- NULL

if(!is.null(peaklist)){
    for(i in 1:length(peaklist)){
       position <- grep(names(peaklist[i]), exons$EnsemblExonID)
       peakpositions <- c(peakpositions, position)
       }

    names(peakpositions) <- names(peaklist)
    peakout <- data.frame(peaklist, peakpositions)

    }

if(is.null(peaklist)){ peakout <- data.frame(c(NA),c(NA))
colnames(peakout) <- c("peaklist", "peakpositions")
rownames(peakout) <- c("NA")
}

peakout

}

FindExonShutoffs <- function(translist, ExpSet,threshold=5, uniqueTargets=T, noOverhang=F) {

     if(dbenv[["connected"]]==FALSE){
       stop("Please connect to a database using mapConnect()")
     }
  
     trLength <- length(translist)
     shutoffs <- NULL     
     exonNames <- NULL
     notfound <- NULL
     trout <- NULL

     for(i in 1:trLength) {
		outstring <- paste("Running %", i, ": ", translist[i],sep="")
		print(outstring)
		tdat <- fetchdata(translist[i], ExpSet)

		if(is.null(tdat)){notfound <- c(notfound, i)
		}

		if(!is.null(tdat)){
		   res <- getExonPeaks(tdat,threshold=threshold, uniqueTargets=uniqueTargets, noOverhang=noOverhang)
		   if(!is.null(res)){		
		      res <-  res[order(res$peaklist),]
		      }

		if(!is.na(res[1]$peaklist)){
 		   shutoffs <- rbind(shutoffs, res)
		   exonNames <- c(exonNames, rownames(res))
		   trout <- c(trout, rep(translist[i],nrow(res)))
               }
		}							
		
		}
   

   outframe <- data.frame(trout, exonNames, shutoffs)
   colnames(outframe)[1:2] <- c("EnsemblTranscriptID", "EnsemblExonID") 
   outframe
}



plotExonData <-  function(transdat,analysisUnit="exon",savePlot=TRUE,display=FALSE, 
	uniqueTargets=T, noOverhang=F, height=640, width=1024) {
  
      if(class(transdat)!="data.frame"){stop("Input is not data frame")}

      Id <- as.character(transdat$transid[1])

#grab all sorted exons for ID
      if(length(grep("T", Id))==1){ exons <-getExonLocationsTranscript(trId=Id, 
		uniqueTargets=uniqueTargets, noOverhang=noOverhang ) }
      else{exons <- getExonLocationsGene(geneId=Id, uniqueTargets=uniqueTargets, 
		noOverhang=noOverhang)}

#attach(transdat)


Exon <- transdat$Exon
      
#subset exon information to only those which exist in dataset
exons <- exons[which(exons$EnsemblExonID %in% Exon),]

Exonorder <- factor(Exon, levels=exons$Ensembl)

transdat <- merge(exons, transdat, by.x=1, by.y=5)


if("Start" %in% colnames(transdat)){
	if(exons$Strand[1] == 1){
		transdat <- transdat[order(transdat$Start),]

	}

	if(exons$Strand[1] == -1){
		transdat <- transdat[order(-transdat$Start),]
	}
}

else{
	if(exons$Strand[1] == 1){
		transdat <- transdat[order(transdat$ExonStart),]
	}

	if(exons$Strand[1] == -1){
		transdat <- transdat[order(-transdat$ExonStart),]
	}

}

Exon <- transdat$Exon
Probeset <- transdat$probesetid

#print(transdat)

ProbesetExon <- factor(paste(Exon, Probeset, sep="-"), levels = unique(paste(Exon, Probeset, sep="-")))

#print(ProbesetExon)

Probeset <- factor(Probeset, levels=unique(transdat$probesetid))


Expression <- transdat$Expression
Subject <- transdat$Subject
Strain <- transdat$Strain


sym <- getGeneSym(Id)

plotname <- paste(sym, "-", Id, sep = "")

filename <- paste(sym, "-", Id, ".exon.plot.png",sep="")
      
if(savePlot){
      png(filename, height=height, width=width)

      par(mar=c(12,3,1,1), las=2) 
	if(analysisUnit == "exon"){
		lineplot.CI(Exon, Expression, group=Strain, xlab="", ylab="Expression", main=as.character(plotname))
	}

	if(analysisUnit =="probeset"){
		lineplot.CI(ProbesetExon, Expression, group=Strain, xlab="", ylab="Expression", main =as.character(plotname))
		#axis(1, Exon)
	}

	#dev.print(file=filename, device=jpeg, height=640, width=900)

      dev.off()
}
 
if(display){

par(mar=c(10,3,1,1), las=2) 

if(analysisUnit=="exon"){
lineplot.CI(Exon, Expression, group=Strain, xlab="", ylab="Expression", main=as.character(plotname))
}

if(analysisUnit=="probesetid"){
lineplot.CI(ProbesetExon, Expression, group=Strain, xlab="", ylab="Expression", main =as.character(plotname))
}

}

      #clean up
 #     detach(transdat)
      rm(Exonorder)

if(savePlot)
         {workdir <- getwd()
          mess <- paste("plot is available as ", workdir, "/", filename, sep = "")
          print(mess)
   }
}


PlotExon <- function(id, ExpSet, display=TRUE, savePlot=FALSE, uniqueTargets=T, noOverhang=F, analysisUnit="exon") {
      if(class(id)!="character"){ stop("input is not id")}

#if(is.null(dBPackage)){ stop("Database Package needs to be specified") }
#else{library(dBPackage, character.only=TRUE)}
      if(dbenv[["connected"]]==FALSE){
        stop("Please connect to a database using mapConnect()")
      }

  #dBPackage <- dbenv[["dbPackage"]]

      
      if(class(ExpSet)!="ExpressionSet"){stop("input is not ExpressionSet")}

trdat <- fetchdata(id, ExpSet,uniqueTargets=uniqueTargets, noOverhang=noOverhang, grabProbesetStart=TRUE)
plotExonData(trdat, savePlot=savePlot, display=display, analysisUnit=analysisUnit)

}



PlotExonResults <- function(resultframe, ExpSet, savePlot=TRUE, uniqueTargets=T,noOverhang=F, analysisUnit="exon")
        {

          if(dbenv[["connected"]]==FALSE){
            stop("Please connect to a database using mapConnect()")
          }
          
	if(class(resultframe)!="data.frame"){
		stop("input is not a dataframe")
	}

	  dBPackage <- dbenv[["dbPackage"]]


        idlist <- as.character(resultframe[,1])
	for(id in idlist){
	     print(id)
	     trdat <- fetchdata(id,ExpSet, uniqueTargets=uniqueTargets)
#, noOverhang=noOverhang)
	     if(!is.null(trdat)){
	     plotExonData(trdat, savePlot=savePlot, uniqueTargets=uniqueTargets, analysisUnit=analysisUnit)
#, noOverhang=noOverhang)
           }
      }
}
