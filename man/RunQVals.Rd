\name{RunQVals}
\alias{RunQVals}
\title{Run qvalue on a resultframe}

\description{For a resultframe produced by \code{RunExonModelWorkflow}, produce a new
resultframe containing qvalues (FDR adjusted pvalues). }

\usage{
  RunQVals(resultframe)
}

\arguments{
  \item{resultframe}{
    Valid resultframe (either multi or singles) produced from \code{RunExonModelWorkflow}.  
  }
}

\details{
  \code{RunQVals} produces a dataframe that also contains the qvalues.  

  Note that if \code{qvalue} cannot produce qvalues for that set of pvalues, the function 
  will return a column of NAs for that set of pvalues, as well as the warning from qvalues.
}

\seealso{
  \code{\link{RunExonModelWorkflow}}
  \code{\link[qvalue]{qvalue}}
}

\examples{
  #mouseexonensembl.db package required
  #load in sample dataset
  data(exontestdata)
  #show output of RunExonModelWorkflow
  multi[1:5]
  output <- RunQVals(multi)
  output

  ##not run
  ##usually, would run with the output of RunExonModelWorkflow as below:
  ##results <- RunExonModelWorkflow(TestSetTrans, testTrans, dBPackage="mouseexonensembl.db")
  ##this example will not work - not enough pvalues to build a distribution 
  ##9 out of the 20 transcripts are multiple-exon transcripts
  ##run on trx with multiple exons
  ##mult <- RunQVals(results$multi)
  ##mult
  ##run on trx with single exons
  ##singles <- RunQVals(results$single)
  ##singles
}

\keyword{models}
