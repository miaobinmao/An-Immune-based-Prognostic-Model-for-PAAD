#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
rm(list = ls())
library(limma)
#setwd(" ")         #设置工作目录

#读取GTEx文件，并整理数据
rt1=read.table("GTExNormalExp.txt",sep="\t",header=T,check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=avereps(data1)
data1=data1[rowMeans(data1)>0,]

#读取TCGA文件，并整理数据
rt2=read.table("tcgaSymbol.txt",sep="\t",header=T,check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)
data2=data2[rowMeans(data2)>0,]
#对基因取交集
sameGene=intersect( row.names(data1),row.names(data2) )
data=cbind(data1[sameGene,],data2[sameGene,])

#数据矫正
outTab=normalizeBetweenArrays(data)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="merge.txt",sep="\t",quote=F,col.names=F)

