- 对转录组数据做层次聚类
- 1.获得转录组数据的FPKM,RPKM,TPM
- 2.对数据进行去批次效应处理，使用R语言sva包中的ComBat
- #2.1生成一个batch文件，每一组作为一个batch标记
```R
csif <- read.table("sif.txt", header = T, sep = "\t", row.names = 1)
modcombat = model.matrix(~1, data = csif)
batch = csif$batch
combat_edata = ComBat(dat=cdata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)​
tcombat <- combat_edata
out.dist=dist(tcombat,method="euclidean")
out.hclust=hclust(out.dist,method="complete")
plot(out.hclust,hang = -1,main="euclidean")
```