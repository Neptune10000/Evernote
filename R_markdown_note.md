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


#R


* 从一个文件夹中逐个读取文件
```R
temp = list.files(pattern = "*.tsv")
for (i in 1:length(temp)) assign(temp[i],read.table(temp[i]))
```

* 矩阵合并
```{r}
cbind(a,b)

gene_FPKM<-cbind(`adrenal gland_ENCFF309JAC.tsv`[,0],`adrenal gland_ENCFF309JAC.tsv`[,1],`adrenal gland_ENCFF309JAC.tsv`[,7],`adrenal gland_ENCFF578HON.tsv`[,7],`bipolar neuron_ENCFF249CKU.tsv`[,7],`bipolar neuron_ENCFF688CEH.tsv`[,7],`cerebellum_ENCFF902DQW.tsv`[,7],`cerebellum_ENCFF947JHW.tsv`[,7],`diencephalon_ENCFF056KYM.tsv`[,7],`diencephalon_ENCFF224CON.tsv`[,7],`frontal cortex_ENCFF847LUK.tsv`[,7],`frontal cortex_ENCFF947SAO.tsv`[,7],`heart_ENCFF455HQT.tsv`[,7],`heart_ENCFF803YED.tsv`[,7],`lung_ENCFF573QQU.tsv`[,7],`lung_ENCFF649ECE.tsv`[,7],`mesenchymal stem cell_ENCFF290OQE.tsv`[,7],`mesenchymal stem cell_ENCFF693WRN.tsv`[,7],`metanephros_ENCFF397JTM.tsv`[,7],`metanephros_ENCFF904AZG.tsv`[,7],`occipital lobe_ENCFF450SDZ.tsv`[,7],`occipital lobe_ENCFF757CQE.tsv`[,7],`ovary_ENCFF132XQU.tsv`[,7],`ovary_ENCFF809AOV.tsv`[,7],`parietal lobe_ENCFF637ZPY.tsv`[,7],`parietal lobe_ENCFF924PET.tsv`[,7],`skeletal muscle tissue.tsv`[,7],`skeletal muscle tissue_2.tsv`[,7],`skin of body.tsv`[,7],`skin of body_2.tsv`[,7],`spinal cord_ENCFF126UNL.tsv`[,7],`spinal cord_ENCFF316BNE.tsv`[,7],`spleen_ENCFF474KYX.tsv`[,7],`spleen_ENCFF653LOC.tsv`[,7],`stomach.tsv`[,7],`stomach_2.tsv`[,7],`stomach_ENCFF046GFN.tsv`[,7],`stomach_ENCFF645CNL.tsv`[,7],`temporal lobe_ENCFF682JEP.tsv`[,7],`temporal lobe_ENCFF845GMQ.tsv`[,7],`testis_ENCFF857HXK.tsv`[,7],`testis_ENCFF863ERP.tsv`[,7],`thyroid gland_ENCFF734HNL.tsv`[,7],`thyroid gland_ENCFF860IOK.tsv`[,7],`tibial nerve_ENCFF046AFQ.tsv`[,7],`tibial nerve_ENCFF470ZWQ.tsv`[,7],`tibial nerve_ENCFF597HTT.tsv`[,7],`tongue.tsv`[,7],`tongue_2.tsv`[,7],`uterus.tsv`[,7],`uterus_2.tsv`[,7])

```
矩阵大小

```R
dim(a)


```

* 矩阵输出
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

* Heatmap
```R
library(RColorBrewer)
library(gplots)
table <- read.table("miRNA_expression.txt", header=TRUE)
display.brewer.all()
brewer.pal(11, "RdBu")
colorRampPalette(brewer.pal(11, "RdBu"))
colorRampPalette(brewer.pal(11, "RdBu"))(256)
################################
pdf("gene_pickout_heatmap.pdf")
heatmap.2(as.matrix(log(table[,2:29] + 1) / log(2)),labRow=table[,1], scale="row", distfun=function(x) as.dist((1-cor(t(x)))/2), trace="none", col=rev(colorRampPalette(brewer.pal(11, "BrBG"))(256)), margins=c(12,8))
dev.off()
################################
#pdf("minkowski.pdf")
#heatmap.2(as.matrix(log(table[,2:29] + 1) / log(2)),labRow=table[,1], scale="row", distfun=function(x) dist(x,method='minkowski'), trace="none", col=rev(colorRampPalette(brewer.pal(11, "BrBG"))(256)), margins=c(12,8))
#dev.off()
###############################
```


