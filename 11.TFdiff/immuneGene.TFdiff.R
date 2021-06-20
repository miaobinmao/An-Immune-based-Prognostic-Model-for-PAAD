

#install.packages("pheatmap")

#setwd("  ")                 #设置工作目录
fdrFilter=0.01                                                    #fdr临界值
logFCfilter=2                                                     #logFC临界值
conNum=171                                                         #normal组样品数目
treatNum=172                                                      #tumor组样品数目

rt=read.table("all.txt",sep="\t",header=T,check.names=F,row.names=1)
diffExp=read.table("diffGeneExp.txt",sep="\t",header=T,check.names=F,row.names=1)
gene=read.table("TF.txt",sep="\t",header=F)
immuneDiffAll=rt[intersect(gene[,1],row.names(rt)),]
immuneDiffGene=intersect(gene[,1],row.names(diffExp))
hmExp=diffExp[immuneDiffGene,]
immuneDiffResult=immuneDiffAll[immuneDiffGene,]

#输出TF的差异结果
immuneDiffResult=rbind(ID=colnames(immuneDiffResult),immuneDiffResult)
write.table(immuneDiffResult,file="TFdiff.xls",sep="\t",col.names=F,quote=F)

#输出差异TF的表达量
immuneGeneExp=rbind(ID=colnames(hmExp),hmExp)
write.table(immuneGeneExp,file="TFgeneExp.txt",sep="\t",col.names=F,quote=F)

#绘制火山图
pdf(file="vol.pdf",height=5,width=5)
xMax=max(abs(as.numeric(as.vector(immuneDiffAll$logFC))))
yMax=max(-log10(immuneDiffAll$fdr))+1
plot(as.numeric(as.vector(immuneDiffAll$logFC)), -log10(immuneDiffAll$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(immuneDiffAll, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=0.8)
diffSub=subset(immuneDiffAll, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

#绘制差异基因热图
library(pheatmap)
hmExp=log2(hmExp+0.001)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(diffExp)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=12,width=15)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()
