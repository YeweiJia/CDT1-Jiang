setwd("/Users/jiayewei/Documents/投稿/GSE46517 Necrosis/8. Estimate")
rm(list=ls())

library(ggsci)
library(tidyr)
library(ggpubr)
library(utils)
library(estimate)
library(tidyverse)

expr<-read.table("SKCM根据LILRB3分组.txt",sep="\t",row.names = 1,check.names = F,
                 stringsAsFactors = F,header = T)

##计算免疫评分
filterCommonGenes(input.f = "SKCM根据LILRB3分组.txt", ##输入文件名
                  output.f = "SKCM根据LILRB3分组.gct", ###输出文件名
                  id="GeneSymbol")
estimateScore("SKCM根据LILRB3分组.gct",###刚才的输出文件名
              "SKCM-LILRB3 estimate score.txt",
              platform = "affymetrix")##评分文件
##3输出每个样本打分
result<-read.table("SKCM-LILRB3 estimate score.txt",sep="\t",row.names = 1,check.names = F,
                   stringsAsFactors = F,header = T)
result<-result[,-1]
colnames(result)<-result[1,]              
result<-as.data.frame(t(result[-1,]))

write.table(result,file="SKCM-LILRB3 estimate score.txt",sep="\t",row.names = T,col.names = NA,
            quote=F) ##保存并覆盖文件
