#install.packages("glmnet")
#install.packages("survival")

rm(list=ls())
library("glmnet")
library("survival")

coxSigFile="uniSigExp.txt"       #输入文件
#setwd("   ")            #设置工作目录
rt=read.table(coxSigFile,header=T,sep="\t",row.names=1)              #读取文件
rt$futime=rt$futime/365

#构建模型
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 10000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 10000)

#输出相关基因系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

#输出训练组风险值
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table("GSE71729.expTime.txt",header=T,sep="\t",row.names=1)
#rt$futime=rt$futime/365
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="GSE71729Risk.txt",sep="\t",quote=F,row.names=F)


