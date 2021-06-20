#GEO芯片数据的探针ID转换
rm(list = ls())
#setwd("   ")
## 读取GEO数据
exprSet <- read.table("GSE71729_series_matrix.txt",comment.char="!",stringsAsFactors=F,header=T)
#本矩阵无需id转换

#合并删去重复值
library(dplyr)
library(tidyr)
library(tibble)


exprSet <- exprSet %>%
    #求出平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% #mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(真的是画龙点睛)
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>%
  # 去重，symbol留下第一个
  distinct(ID_REF,.keep_all = T) %>%
  #反向选择去除rowMean这一列
  select(-rowMean)%>%
#列名转行名
column_to_rownames("ID_REF")
  
write.table(exprSet,file = "geneMtrix.txt",sep="\t",quote=F,col.names=T)


