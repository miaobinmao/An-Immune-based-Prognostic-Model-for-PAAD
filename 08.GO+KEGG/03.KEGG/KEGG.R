
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("enrichplot")

rm(list = ls())
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(enrichplot)
library(GOplot)
library(DOSE)
library(stringr)
library(export)#输出图片
setwd("    ")             #设置工作目录
rt=read.table("immuneDiff.txt",sep="\t",check.names=F,header=T)    

genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]      #去除基因id为NA

pvalueFilter=0.05        #p值过滤条件
qvalueFilter=0.05           #矫正后的p值过滤条件，结果太少的话改为1（原0.05）
#KEGG富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$ID[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)


#气泡图
pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 20)
dev.off()


####通路-基因网络图 
library(dplyr)
#首先转为数据框
KEGG <- KEGG %>% 
  arrange(p.adjust)
KEGG <- KEGG[1:10,]
KEGGOplotIn<-KEGG[,c(1,2,6,8)] #我们先提取KEGG提取ID,Description,p.adjust,GeneID四列。
KEGGOplotIn$geneID <-str_replace_all(KEGGOplotIn$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(KEGGOplotIn)<-c('ID','Term','adj_pval','Genes')
KEGGOplotIn$Category="KEGG"
#读取数据
fr = data.table::fread(file = "immuneDiff.txt",data.table = F)
genedata<-data.frame(ID=fr$ID,logFC=fr$logFC)
result<-GOplot::circle_dat(KEGGOplotIn,genedata) 

#导出
write.table(result,file="KEGGnetwork.txt",sep="\t",quote=F,row.names = F)

#输出调控网络节点属性
immuneGeneUp=result[result$logFC>0,]
immuneGeneDown=result[result$logFC<0,]
KEGGLabel=cbind(rownames(KEGG),"KEGG")
immuneGeneupLabel=cbind(immuneGeneUp$genes,"UP")
immuneGenedownLabel=cbind(immuneGeneDown$genes,"DOWN")
nodeLabel=rbind(c("ID","type"),KEGGLabel,immuneGeneupLabel,immuneGenedownLabel)
write.table(nodeLabel,file="nodeType.txt",sep="\t",quote=F,col.names=F,row.names=F)



