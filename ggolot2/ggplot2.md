#ggplot2

[TOC]

#PCA画图

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
