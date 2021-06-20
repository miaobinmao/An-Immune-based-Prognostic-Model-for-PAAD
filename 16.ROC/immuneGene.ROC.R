

#install.packages("survivalROC")
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if(!requireNamespace("timeROC",quietly = TRUE)) install.packages("timeROC",update = F,ask = F)
if(!requireNamespace("timereg",quietly = TRUE)) install.packages("timereg",update = F,ask = F)
library(survival)
library(survminer)
library(timeROC)
library(survivalROC)
#setwd("   ")      #设置工作目录
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取cox回归风险文件
bioROC=function(inputFile=null,rocFile=null){
  
  rt=read.table("risk.txt",header=T,sep="\t")
  
  ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
                 marker=rt$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green",'blue','red'),lwd=2,bty = 'n')
  dev.off()
}

bioROC(inputFile="risk.txt",rocFile="ROC.pdf")
