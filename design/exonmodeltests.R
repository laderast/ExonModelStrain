setwd("c:/Users/laderast/Desktop/reference/exonarray/")
load("exondatatest.RData")
library(ExonModelStrainXmap)

mapConnect("exonmap", "mouse")
trans <- getAllTranscripts()

test <- RunExonModelWorkflow(ExSet, trans[1:100])

multi <- RunQvals(test$multi)

multi <- annotateResults(multi)

PlotExon(ExSet, multi[1,1])
