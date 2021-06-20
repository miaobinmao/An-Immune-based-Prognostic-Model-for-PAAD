rm(list = ls())
#install.packages("digest")
#install.packages("GOplot")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(enrichplot)
library(GOplot)
library(DOSE)
library(stringr)
library(export)#输出图片

setwd("D:\\immuneGene\\PAAD\\IRGS_GO+KEGG\\GO\\GOBubble")             

ego=read.table("GO.txt", header = T,sep="\t",check.names=F)          
go=data.frame(Category = ego$ONTOLOGY,ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)


id.fc <- read.table("id.txt", header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$ID, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
circ <- circle_dat(go, genelist)

#绘制GO圈图
pdf(file="GOCircle.pdf",width = 11,height = 6)
GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=5)            #nsub=10中10代表显示GO的数据，可修改
dev.off()

