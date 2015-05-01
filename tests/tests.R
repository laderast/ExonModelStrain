testfetchdata <- function() {
   require(mouseexonensembl.db)
   m <- matrix(data=c(1,1,0,0),nrow=2)
   data("exontestdata.RData", package="ExonStrainModel")
   checkException(fetchdata("ENSMUST00000062893",ExpSet=m))
   checkException(fetchdata("ENSMUSG00000025651",ExpSet=m))
   checkEquals(nrow(fetchdata("ENSMUST00000081945",ExpSet=TestSetTrans)),88)
   checkEquals(nrow(fetchdata("ENSMUSG00000047641",ExpSet=TestSetGene)),88)
   checkEquals(fetchdata("blah", ExpSet=TestSetTrans), NULL)
}

testrunModel <- function() {
    data("exontestdata.RData", package="ExonStrainModel")
    m <- matrix(data=c(1,1,0,0),nrow=2)
    checkException(runModel(m))
}

testgetExonDeltas <- function() {
    m <- matrix(data=c(1,1,0,0),nrow=2)
    checkException(runModel(m))
    data("exontestdata.RData", package="ExonStrainModel")
    checkTrue(is.numeric(getExonDeltas(tdatt)$deltas))
    checkTrue(is.numeric(getExonDeltas(tdatg)$deltas))
    checkEquals(length(getExonDeltas(tdatt)),9)
    checkEquals(length(getExonDeltas(tdatg)),9)
}

testplotExonData <- function() {
    m <- matrix(data=c(1,1,0,0),nrow=2)
    checkException(plotExonData(m))
    data("exontestdata.RData", package="ExonStrainModel")
    
}

#integration tests
testRunExonModelWorkflow <- function(){

  
   m <- matrix(data=c(1,1,0,0),nrow=2)
   data("exontestdata.RData", package="ExonStrainModel")
   library(mouseexonensembl.db)
   checkException(RunExonModelWorkflow(ExpSet=m))
   checkException(RunExonModelWorkflow(ExpSet=TestSetTrans, analysisType="nottype"))
   checkException(RunExonModelWorkflow(ExpSet=TestSetGene, idlist = c(1,2,3)))
   checkEquals(length(RunExonModelWorkflow(ExpSet=TestSetTrans, idlist=testTrans[1:5],dBPackage="mouseexonensembl.db")$multi),9)
             
}
