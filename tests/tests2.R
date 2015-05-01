library(ExonModelStrainXmap)
mapConnect(dbPackage="xmapcore")
trs <- getAllTranscripts()
res <- RunExonModelWorkflow(eset2, idlist=trs[1643:1650], analysisUnit="probeset")

res <- RunExonModelWorkflow(eset2, idlist=trs[500:510], analysisUnit="probeset")

resmulti <- NULL
ressingles <- NULL
i <- 1

while(i < length(trs)){

	minires <- RunExonModelWorkflow(eset2, idlist=trs[i:(i+500)], analysisUnit="probeset")
	resmulti <- rbind(resmulti, minires$multi)
	ressingles <- rbind(ressingles, minires$singles)
	i <- i + 501
	write.table(resmulti, file="probeset-model-multi.txt", quote=F, row.name=F, sep="\t")
	write.table(ressingles, file="probeset-model-singles.txt", quote=F, row.name=F, sep="\t")
}

write.table(out2, file="probeset-model-multi.txt", quote=F, row.name=F, sep="\t")


#out2 <- xmap_annotate(out)
#write.table(out2, file="probeset-model-multi.txt", quote=F, row.name=F, sep="\t")
