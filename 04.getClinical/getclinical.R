options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
rm(list = ls())


library(TCGAbiolinks)
library(dplyr)
library(DT)


request_cancer=c("PAAD")
i = request_cancer[1] 
cancer_type=paste("TCGA",i,sep="-")
print(cancer_type)
query <- GDCquery(project = cancer_type, data.category = 'Clinical', file.type = 'xml')
GDCdownload(query)
clinical <- GDCprepare_clinic(query,clinical.info = 'patient')
GDCdownload(query)
dir.create(cancer_type)
write.csv(clinical,file = paste(cancer_type,"/",cancer_type,"-clinical.csv",sep = ""))
