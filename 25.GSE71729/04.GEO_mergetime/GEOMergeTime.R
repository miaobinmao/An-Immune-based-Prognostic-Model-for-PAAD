rm(list = ls())
library(limma)

expFile="normalize.txt"       
geneFile="intersectGenes.txt"    
cliFile="time.txt"             
#setwd("   )   

#读取geo标准化后表达矩阵
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#读取交集基因（此处EMC2改为别名TTC35）
rownames(data)=gsub("CARS","CARS1",rownames(data))
geneRT=read.table(geneFile,header=F,sep="\t",check.names=F)
data=data[as.vector(geneRT[,1]),]
data <- t(data)
#读取生存数据
cli=read.table(cliFile,header=T,sep="\t",check.names=F,row.names=1)    

#数据合并取交集输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="GSE71729.expTime.txt",sep="\t",row.names=F,quote=F)



