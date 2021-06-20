

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
setwd("D:\\immuneGene\\PAAD\\IRGS_GO+KEGG\\GO")                     #设置工作目录
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           #读取id.txt文件
rt=rt[is.na(rt[,"entrezID"])==F,]                                 #去除基因id为NA的基因
gene=rt$entrezID

#GO富集分析
GO<- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果


###制作和弦图  
library(dplyr)
#首先转为数据框
GO=as.data.frame(GO)
GO <- GO %>% 
  arrange(p.adjust)
GOplotIn<-GO[1:5,] 
go=data.frame(Category = GOplotIn$ONTOLOGY,ID = GOplotIn$ID,Term = GOplotIn$Description, Genes = gsub("/", ", ", GOplotIn$geneID), adj_pval = GOplotIn$p.adjust)
#构建基因表达矩阵
#读取数据
fr = data.table::fread(file = "immuneDiff.txt",data.table = F)
genedata<-data.frame(ID=fr$ID,logFC=fr$logFC)
row.names(genedata)=genedata[,1]
circ <- circle_dat(go, genedata)
write.table(circ,file="circ.txt",sep="\t",quote=F,row.names = F)
chord<-chord_dat(data = circ,genes = genedata) #生成含有选定基因的数据框

###好了，到此数据都整理好了，开始画和弦图。
GOChord( #GO富集和弦图，fig9
  data = chord,
  title = 'GOchord plot',
  space = 0,#GO Term间距
  limit = c(1,1),
  gene.order = 'logFC',
  gene.space = 0.25,
  gene.size = 5,
  lfc.col = c('red','white','blue'), #上下调基因颜色
  #ribbon.col = brewer.pal(length(go$Term)),#GO Term colors
  process.label = 10) #GO Term字体大小
print(GOChord)
graph2ppt(file="GOChord.pptx", width=20, height=20)

