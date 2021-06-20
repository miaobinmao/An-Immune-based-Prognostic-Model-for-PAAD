
corFilter=0.4              #相关系数过滤标准
pvalueFilter=0.001         #p值过滤标准

#setwd("  ")                     #设置工作目录

#读取文件，并对样品取交集
TF = read.table("TFgeneExp.txt", row.names=1 ,header=T,sep="\t",check.names=F)   #读取TF表达文件
immuneGene = read.table("uniSigExp.txt", row.names=1 ,header=T,sep="\t",check.names=F)   #读取生存相关的immuneGene文件
immuneGene=t(immuneGene[,3:ncol(immuneGene)])
rownames(immuneGene)=gsub("\\|","\\-",rownames(immuneGene))
group=sapply(strsplit(colnames(TF),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
TF=TF[,group==0]
colnames(TF)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(TF))
sameSample=intersect(colnames(TF),colnames(immuneGene))
TF1=TF[,sameSample]
immuneGene1=immuneGene[,sameSample]

#相关性检验
outTab=data.frame()
for(i in row.names(TF1)){
  if(sd(TF1[i,])>1){
    for(j in row.names(immuneGene1)){
      x=as.numeric(TF1[i,])
      y=as.numeric(immuneGene1[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(TF=i,immuneGene=j,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(TF=i,immuneGene=j,cor,pvalue,Regulation="negative"))
      }
    }
  }
}
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)        #输出相关性结果

#输出调控网络节点属性
immuneGeneSig = read.table("uniCox.txt", row.names=1 ,header=T,sep="\t",check.names=F)
rownames(immuneGeneSig)=gsub("\\|","\\-",rownames(immuneGeneSig))
immuneGeneUp=immuneGeneSig[immuneGeneSig$HR>1,]
immuneGeneDown=immuneGeneSig[immuneGeneSig$HR<1,]
TFLabel=cbind(rownames(TF),"TF")
immuneGeneupLabel=cbind(rownames(immuneGeneUp),"highRiskIRG")
immuneGenedownLabel=cbind(rownames(immuneGeneDown),"lowRiskIRG")
nodeLabel=rbind(c("ID","type"),TFLabel,immuneGeneupLabel,immuneGenedownLabel)
write.table(nodeLabel,file="nodeType.txt",sep="\t",quote=F,col.names=F,row.names=F)
