### 设置镜像，存在第二种方式
#options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
   #install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("scatterplot3d")
rm(list = ls())
library(limma)
library(scatterplot3d)

#setwd("C:CA")                  

myPCA=function(input=null,output=null)
{
		#读取基因表达文件
		rt=read.table(input,sep="\t",header=T,check.names=F)
		rt=as.matrix(rt)
		rownames(rt)=rt[,1]
		exp=rt[,2:ncol(rt)]
		dimnames=list(rownames(exp),colnames(exp))
		data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
		data=avereps(data)
		data=data[rowMeans(data)>0.5,]
		type=sapply(strsplit(colnames(data),"\\-"),"[",4)
		type=sapply(strsplit(type,""),"[",1)
		type=gsub("2","1",type)
		data=t(data[,type==0])
		rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
		
		#读取风险文件
		risk=read.table("tcgaRisk.txt",sep="\t",header=T,row.names=1)
		sameSample=intersect(rownames(data),rownames(risk))
		data=data[sameSample,]
		risk=risk[sameSample,]
		group=as.vector(risk[,"risk"])
		
		#PCA分析
		data.class <- rownames(data)
		data.pca <- prcomp(data, scale. = TRUE)
		
		#可视化
		color=ifelse(group=="low",3,2)
		pcaPredict=predict(data.pca)
		pdf(file=output,width=5.5,height=5)
		s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
		legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
		dev.off()
}

#所有基因
myPCA(input="merge.txt",output="allGene.PCA.pdf")
#免疫基因
myPCA(input="immuneGeneExp.txt",output="immuneGene.PCA.pdf")


#模型相关基因绘制PCA
risk=read.table("tcgaRisk.txt",sep="\t",header=T,row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])
		
#PCA分析
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)
		
#可视化
color=ifelse(group=="low",3,2)
pcaPredict=predict(data.pca)
pdf(file="riskGene.PCA.pdf",width=5.5,height=5)
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
dev.off()
