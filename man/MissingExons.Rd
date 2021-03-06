\name{MissingExons}
\alias{MissingExons}
\title{Find exons for a transcript/gene that are not mapped to the dataset}

\description{For a given transcript, or gene id, return all Ensembl Exon IDs that do not have
mapped probesets.}

\usage{
  MissingExons(id,ExpSet,dBPackage="mouseexonensembl.db")
  MissingExons(id,ExpSet)
}

\arguments{
  \item{id}{
    Valid Ensembl Transcript or Gene ID.  
    
    Note a list of all valid Transcript or Gene IDs can be queried from the database 
    package by using getAllTranscripts() or getAllGenes().
  }
  \item{ExpSet}{
    ExpressionSet containing probeset-level Exon Expression Data.  Currently, only
    Mouse Exon 1.0 data is supported.  The ExpressionSet slot phenoData should contain a
    column called 'Strain' where the two strains are coded 1 and 2.
  }
  \item{dbPackage}{
    Name of database package containing Ensembl-Exon Array mapping (currently only 
    mouseexonensembl.db exists).
  }
}

\details{
  \code{MissingExons} produces a list of Ensembl Exon IDs not mapped to the current
  Expression Set. This may be for two reasons: either the Exon ID does not have any
  probesets mapped to it, or these probesets do not exist in the current ExpressionSet
  (perhaps due to SNP masking, etc).

  If all existing Exons are mapped, \code{MissingExons} returns a NULL value.
}

\seealso{
  \code{\link{RunExonModelWorkflow}}
}

\examples{
  data(exontestdata)
  res <- MissingExons("ENSMUST00000062893",TestSetTrans)
  print(res)
}

\keyword{models}
\keyword{graphs}


