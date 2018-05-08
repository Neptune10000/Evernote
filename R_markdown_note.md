---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
plot(cars)
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).



-------------

#R
[TOC]

## 查看变量类型
mode(x)

## 从一个文件夹中逐个读取文件
```{r}
temp = list.files(pattern = "*.tsv")
for (i in 1:length(temp)) assign(temp[i],read.table(temp[i]))
```

## 矩阵合并
```{r}
cbind(a,b)

gene_FPKM<-cbind(`adrenal gland_ENCFF309JAC.tsv`[,0],`adrenal gland_ENCFF309JAC.tsv`[,1],`adrenal gland_ENCFF309JAC.tsv`[,7],`adrenal gland_ENCFF578HON.tsv`[,7],`bipolar neuron_ENCFF249CKU.tsv`[,7],`bipolar neuron_ENCFF688CEH.tsv`[,7],`cerebellum_ENCFF902DQW.tsv`[,7],`cerebellum_ENCFF947JHW.tsv`[,7],`diencephalon_ENCFF056KYM.tsv`[,7],`diencephalon_ENCFF224CON.tsv`[,7],`frontal cortex_ENCFF847LUK.tsv`[,7],`frontal cortex_ENCFF947SAO.tsv`[,7],`heart_ENCFF455HQT.tsv`[,7],`heart_ENCFF803YED.tsv`[,7],`lung_ENCFF573QQU.tsv`[,7],`lung_ENCFF649ECE.tsv`[,7],`mesenchymal stem cell_ENCFF290OQE.tsv`[,7],`mesenchymal stem cell_ENCFF693WRN.tsv`[,7],`metanephros_ENCFF397JTM.tsv`[,7],`metanephros_ENCFF904AZG.tsv`[,7],`occipital lobe_ENCFF450SDZ.tsv`[,7],`occipital lobe_ENCFF757CQE.tsv`[,7],`ovary_ENCFF132XQU.tsv`[,7],`ovary_ENCFF809AOV.tsv`[,7],`parietal lobe_ENCFF637ZPY.tsv`[,7],`parietal lobe_ENCFF924PET.tsv`[,7],`skeletal muscle tissue.tsv`[,7],`skeletal muscle tissue_2.tsv`[,7],`skin of body.tsv`[,7],`skin of body_2.tsv`[,7],`spinal cord_ENCFF126UNL.tsv`[,7],`spinal cord_ENCFF316BNE.tsv`[,7],`spleen_ENCFF474KYX.tsv`[,7],`spleen_ENCFF653LOC.tsv`[,7],`stomach.tsv`[,7],`stomach_2.tsv`[,7],`stomach_ENCFF046GFN.tsv`[,7],`stomach_ENCFF645CNL.tsv`[,7],`temporal lobe_ENCFF682JEP.tsv`[,7],`temporal lobe_ENCFF845GMQ.tsv`[,7],`testis_ENCFF857HXK.tsv`[,7],`testis_ENCFF863ERP.tsv`[,7],`thyroid gland_ENCFF734HNL.tsv`[,7],`thyroid gland_ENCFF860IOK.tsv`[,7],`tibial nerve_ENCFF046AFQ.tsv`[,7],`tibial nerve_ENCFF470ZWQ.tsv`[,7],`tibial nerve_ENCFF597HTT.tsv`[,7],`tongue.tsv`[,7],`tongue_2.tsv`[,7],`uterus.tsv`[,7],`uterus_2.tsv`[,7])

```
##矩阵大小

```{r}
dim(a)


```

##矩阵输出
1、导出文本文件
1）write.table函数语法：
write.table (x, file ="", sep ="", row.names =TRUE, col.names =TRUE, quote =TRUE)
x：需要导出的数据
file：导出的文件路径
sep：分隔符，默认为空格（" "），也就是以空格为分割列
row.names：是否导出行序号，默认为TRUE，也就是导出行序号
col.names：是否导出列名，默认为TRUE，也就是导出列名
quote：字符串是否使用引号表示，默认为TRUE，也就是使用引号表示
2）案列：
先生成一个数据框（生成数据！）
age <- c (22,23)
name <- c ("ken", "john")
f <- data.frame (age, name)
再导出数据框数据（导出数据！）
write.table (f, file ="f.csv")    ---以空格分隔数据列（默认），含行号（默认），含列名（默认），字符串带引号
write.table (f,file ="f.csv", sep =",")    ---以逗号分隔数据列，含行号（默认），含列名（默认），字符串带引号
write.table (f,file ="f.csv", sep =",", row.names = FALSE)    ---以逗号分隔数据列，不含行号，含列名（默认），字符串带引号
write.table (f,file ="f.csv", row.names = FALSE, col.names =FALSE)    ---以空格分隔数据列，不含行号，不含列名，字符串带引号
write.table (f,file ="f.csv", row.names = FALSE, col.names =FALSE, quote =FALSE)    ---以空格分隔数据列，不含行号，不含列名，字符串不带引号

##Heatmap
```{r}
library(RColorBrewer)
library(gplots)
table <- read.table("miRNA_expression.txt", header=TRUE)
display.brewer.all()
brewer.pal(11, "RdBu")
colorRampPalette(brewer.pal(11, "RdBu"))
colorRampPalette(brewer.pal(11, "RdBu"))(256)
################################
pdf("gene_heatmap.pdf")
heatmap.2(as.matrix(log(table[,2:29] + 1) / log(2)),labRow=table[,1], scale="row", distfun=function(x) as.dist((1-cor(t(x)))/2), trace="none", col=rev(colorRampPalette(brewer.pal(11, "BrBG"))(256)), margins=c(12,8))
dev.off()
################################
#pdf("minkowski.pdf")
#heatmap.2(as.matrix(log(table[,2:29] + 1) / log(2)),labRow=table[,1], scale="row", distfun=function(x) dist(x,method='minkowski'), trace="none", col=rev(colorRampPalette(brewer.pal(11, "BrBG"))(256)), margins=c(12,8))
#dev.off()
###############################
```
 
##读取csv文件
```{r}
library(readr)
dataset <- read_csv("gene_FPKM.csv")
```


##层次聚类

- 对转录组数据做层次聚类
- 1.获得转录组数据的FPKM,RPKM,TPM
- 2.对数据进行去批次效应处理，使用R语言sva包中的ComBat
    - 2.1生成一个batch文件，每一组作为一个batch标记
- 3.对TPM/FPKM/RPKM取对数，参考heatmap中取对数的方法
    - **注意**：输入数据的类型，excel直接导入为list，但是转置后可能会变成charactor，表格区域不能有字符型变量
- 4.进行层次聚类
    - **注意**：层次聚类的数据格式，行为不同样本，列为不用基因，如果反了内存会爆。
```R
#去除批次效应(可选，根据批次效应的程度而定)
csif <- read.table("sif.txt", header = T, sep = "\t", row.names = 1)
modcombat = model.matrix(~1, data = csif)
batch = csif$batch
combat_edata = ComBat(dat=cdata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)​
tcombat <- combat_edata
#对表达量数据取对数
tgene_expression_No0 <- t(gene_expression_No_0)
mode(tgene_expression_No0)
num<-t(gene_expression_No_0[,2:50])
mode(num)
log_tgene_expression_No0 <- as.matrix(log((num)+1)/log(2))
#层次聚类
## tcombat或log_tgene_expression_No0
out.dist=dist(tcombat,method="euclidean")
out.hclust=hclust(out.dist,method="complete")
plot(out.hclust,hang = -1,main="euclidean")
```
- 层次聚类距离选择
help(dist)
as.dist((1-cor(t(x)))/2)
as.dist((1-cor(t(x),method = pearson)))
1.欧氏距离
2.曼哈顿距离
3.切比雪夫距离
4.闵可夫斯基距离(Minkowski Distance)
5.标准化欧氏距离 (Standardized Euclidean distance )
6.马氏距离(Mahalanobis Distance)
7.巴氏距离（Bhattacharyya Distance）
8.汉明距离(Hamming distance)
9.夹角余弦(Cosine)
10.杰卡德相似系数(Jaccard similarity coefficient)
11.皮尔逊系数(Pearson Correlation Coefficient)

##R语言数据类型转换
函数一：as.character(x)

函数二：as.complex(x)

函数三：as.numeric(x)或者as.double(x)

函数四：as.integer(x)

函数五：as.logical(x)

#ggplot2

##PCA画图

```R
data.matrix <- matrix(nrow=100, ncol=10)
colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep=""))
rownames(data.matrix) <- paste("gene", 1:100, sep="")
for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}
head(data.matrix)
dim(data.matrix)h

pca <- prcomp(t(data.matrix), scale=TRUE) 

## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")


# PCA plot
PCAplot1 <- ggplot(PCA, aes(PC1, PC2)) + aes(colour=Group)+geom_point(size=3) + geom_text(aes(label=rownames(PCA)), size=4, vjust =0.5, hjust=-0.1 ,nudge_x=1,show.legend=F) + scale_colour_discrete(l = 55) + coord_cartesian(xlim=c(-150, 150)) + geom_hline(yintercept=0, linetype="dashed", col="#D3D3D3", size=0.8) + geom_vline(xintercept=0, linetype="dashed", col="#D3D3D3", size=0.8) + labs(title="PCA plot"  ,x="PC1 (38.3%)", y="PC2 (20.9%)") + theme_minimal() + theme(plot.title=element_text(family="Calibri", face="bold", size=20, hjust=0.5), axis.title=element_text(face="bold", size=12), legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold"), legend.key.height=unit(1,"cm"), legend.key.width=unit(0.5, "cm"))

###################################

reads <- read.table("D:/reads.txt")
reads <- read.table("D:/reads.txt",header = TRUE,row.names = 1)
pca <- prcomp(t(reads), scale=TRUE)
plot(pca$x[,1], pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

library(ggplot2)
value <- pca$x
value <- as.data.frame(value)
value2 <- cbind(value[,1:7],Group)
value$PC1
PCAplot1 <- ggplot(value, aes(PC1, PC2)) +geom_point(size=3)


PCAplot1 <- ggplot(value2, aes(PC1, PC2)) + aes(colour=Group)+geom_point(size=3) + geom_text(aes(label=rownames(value2)), size=4, vjust =0.5, hjust=-0.1 ,nudge_x=1,show.legend=F) + scale_colour_discrete(l = 55) + coord_cartesian(xlim=c(-500, 500)) + geom_hline(yintercept=0, linetype="dashed", col="#D3D3D3", size=0.8) + geom_vline(xintercept=0, linetype="dashed", col="#D3D3D3", size=0.8) + labs(title="PCA plot"  ,x="PC1 (38.3%)", y="PC2 (20.9%)") + theme_minimal()

PCAplot2 <- ggplot(value2, aes(PC1, PC2)) + 
aes(colour=Group)+geom_point(size=3) + 
geom_text(aes(label=rownames(value2)), size=4, vjust =0.5, hjust=-0.1 ,nudge_x=1,show.legend=F) + 
scale_colour_discrete(l = 55) + 
coord_cartesian(xlim=c(-500, 500)) + 
geom_hline(yintercept=0, linetype="dashed", col="#D3D3D3", size=0.8) + 
geom_vline(xintercept=0, linetype="dashed", col="#D3D3D3", size=0.8) + 
labs(title="PCA plot"  ,x="PC1 (38.3%)", y="PC2 (20.9%)") + 
theme_minimal() + 
theme(plot.title=element_text(family="Calibri", face="bold", size=20, hjust=0.5), 
      axis.title=element_text(face="bold", size=12), 
      legend.title=element_blank(), 
      legend.text = element_text(size = 12, face = "bold"), 
      legend.key.height=unit(1,"cm"), legend.key.width=unit(0.5, "cm"))

Group<-c('Aortic smooth muscle cell','Aortic smooth muscle cell','Fibroblast of arm','Fibroblast of arm','Fibroblast of lung','Fibroblast of lung','Skeletal muscle myoblast','Skeletal muscle myoblast','Occipital lobe','Occipital lobe','Parietal lobe','Parietal lobe','Temporal lobe','Temporal lobe','Cerebellum ','Cerebellum ','Frontal cortex','Frontal cortex','Diencephalon','Diencephalon','Spinal cord','Spinal cord','Bipolar neuron','Bipolar neuron','Tibial nerve','Tibial nerve','Tibial nerve','Tibial nerve')

```

- **去除背景色和网格**
    - 方法一
        ```
        > PCAplot + 
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        + panel.background = element_blank(),axis.line = element_line(colour = "black"))
        ```
    - 方法二
        ```
        PCAplot + theme_bw() +
        theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        + panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
        ```


# VSC
- 多行编辑：alt + shift + 鼠标左键
- 