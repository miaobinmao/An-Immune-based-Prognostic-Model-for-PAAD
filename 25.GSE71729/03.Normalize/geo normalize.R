

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")

library("limma")
normalCount=46               #正常样品数目
tumorCount=145                #肿瘤样品数目

#setwd("   ")          #设置工作目录
rt=read.table("sampleExp.txt",sep="\t",header=T,check.names=F)           #读取文件
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

group=c(rep("normal",normalCount),rep("tumor",tumorCount))
design <- model.matrix(~factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(data)
out=normalizeBetweenArrays(data)
out=rbind(ID=colnames(out),out)
write.table(out,file="normalize.txt",sep="\t",quote=F,col.names=F)        #输出文件

