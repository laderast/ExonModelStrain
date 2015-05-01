\name{RunExonModelWorkflow}
\alias{RunExonModelWorkflow}
\title{Run Exon Model}
\description{Run Strain Specific Exon Modeling}

\usage{
  RunExonModelWorkflow(ExpSet,idlist=ExpSet,idlist = NULL,analysisType="transcript",dBPackage = "mouseexonensembl.db")
  RunExonModelWorkflow(ExpSet,idlist=NULL)
  RunExonModelWorkflow(ExpSet,analysisType="gene") 
}

\arguments{
  \item{ExpSet}{
    ExpressionSet containing probeset-level Exon Expression Data.  Currently, only
    Mouse Exon 1.0 data is supported.  The ExpressionSet phenoData should contain a
    column called 'Strain' where the two strains are coded 1 and 2.
  }
  \item{idlist}{
    Character Vector containing a list of valid Ensembl Transcript or Gene IDs.  
    Note a list of all Transcript or Gene IDs can be queried from the database 
    package by using getAllTranscripts() or getAllGenes().
 
    if idlist is NULL, then based on analysisType, the function will use 
    getAllTranscripts() or getAllGenes() to obtain a valid gene list.
  }
  \item{analysisType}{
    Character value, must be either "transcript" (for transcript-level analysis)
    or "gene" (gene-level analysis).
  }
  \item{dbPackage}{
    Name of database package containing Ensembl-Exon Array mapping (currently only 
    mouseexonensembl.db exists).
  }
}

\details{Given an ExpressionSet of core expression values for an Affymetrix Exon array, 
  RunExonWorkflow will attempt to model the expression data using one of two models.

  For a multiple-exon transcript/gene, the following model is used.

  \eqn{Expression ~ Strain + Exon + Subject in Strain + Exon:Strain }

  for a single-exon, the model reduces to:

  \eqn{Expression ~ Strain + Subject in Strain }

  For a given Ensembl transcript or gene ID, the function will attempt to gather the information
  for all probesets with associated Exons, and subset the ExpressionSet, producing a data frame
  appropriate for analysis.  The appropriate exon model is then run, and the location and
  identity of the exon with maximum strain difference is returned, along with the appropriate 
  raw p-values for that model as well as other flags and metrics (see section below).

  We strongly suggest the use of a FDR based method such as qvalue to adjust the 
  raw p-values for multiple comparisons before further analysis.}

\value{\code{RunExonModelWorkflow} returns a list with the following objects:
    \item{multi}{data frame that contains the following columns:
	  \item{ID}{Ensembl IDs - for gene level analysis, these are Ensembl Gene IDs.  For
	      transcript-level analysis, these are Ensembl Transcript IDs.}
        \item{pvals}{Three columns corresponding to the raw pvalues for
	  Exon, Strain, Exon:Strain.}
        \item{means}{Strain specific expression means (summarized at either transcript/gene
		level)}
	  \item{modelflag}{flag that indicates whether there are multiple probesets per exon.
		0 = 1 probeset per exon, 1 = at least 1 probeset has multiple probesets per exon.
	      Useful in further stratifying the analysis.}
        \item{maxexon}{Ensemble Exon ID for that corresponds to the maximum exon difference
		between the two strains.}
        \item{maxexondelta}{absolute value of the max exon expression difference between the two strains}
        \item{position}{position in transcript of maximum exon.  Note that for gene-level analysis,
            multiple exons for multiple transcripts may be mapped and thus this column may have no
		real meaning. }
        \item{positionflag}{flag indicating position of maximum exon difference - 1 = beginning of
		transcript (1st exon), 2 = end of transcript (final exon), 0 = middle of transcript. }
        \item{numexonsmapped}{Number of exons mapped in current data set for transcript/gene.  
            Note that this is different than the true number of exons for that transcript/gene.}
        \item{trueexonnum}{True number of exons for transcript.}
        \item{strand}{Strand (either -1 or 1) of transcript/gene.}
	}
    
    \item{singles}{data frame that contains the following columns for single-exon transcripts:
	  \item{ID}{Identical to above.}
        \item{pvals}{Single pvalue corresponding to strain-specific differences.  Note that
		for single-exon transcripts, a slightly different model is run.}
        \item{means}{Identical to above.}
        \item{res$numPrsetsSingleExon}{Number of probesets mapped for single-exon.  Useful for
		further stratifying the analysis.}
        \item{maxexon}{Identical to above.  Obviously, for a single exon transcript, this corresponds
		to the single exon.}
        \item{maxexondelta}{Identical to above.}
        \item{numexonsmapped}{Identical to above.  Obviously, this value is 1 for single-exon}
        \item{trueexonnum}{True number of exons per transcript.  For single-exon transcripts,
		this may not necessarily be 1, as the transcript may have multiple exons, but only
		1 exon is mapped to the data.}
        \item{strand}{identical to above.}
        
    }

    \item{notrun}{Character vector containing those ids that were not run. This
	usually is because the corresponding probeset values for that transcript
	do not exist in the ExpressionSet (due to masking or other reasons). }
  
}


\seealso{
  \code{\link{PlotExon}}
}

\examples{
  #mouseexonensembl.db package required
  library(mouseexonensembl.db)
  #load in sample dataset
  data(exontestdata)
  #show list of Transcript IDs
  testTrans
  results <- RunExonModelWorkflow(TestSetTrans, testTrans, dBPackage="mouseexonensembl.db")
  #9 out of the 20 transcripts are multiple-exon transcripts
  results$multi
  #5 single-exon results in test set
  results$singles
  #some transcripts are not run
  results$notrun
}

\keyword{models}
\keyword{regression}


