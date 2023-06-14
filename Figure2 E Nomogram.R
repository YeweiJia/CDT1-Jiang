setwd("/Users/jiayewei/Desktop/UCSC LUAD/COX")
rm(list=ls())

## 载入R包
library(rms)
library(survival)
#3.读入数据
a=read.table('',header = T,sep = '\t')
head(a)
dd=datadist(a)
options(datadist="dd") 

## 构建logist模型,绘制诺莫图
## 构建logist模型,绘制诺莫图
f <- lrm(OS ~ age + pathologic_M + pathologic_N+pathologic_T+ CDT1 , data = a)
nom <- nomogram(f, fun=plogis, lp=F, funlabel="Risk")
plot(nom)

###COX回归中位生存时间的Nomogram
## 构建COX比例风险模型
f2 <- psm(Surv(OS.time,OS) ~ age + pathologic_M + pathologic_N+pathologic_T+ CDT1,data =  a, dist='lognormal') 
med <- Quantile(f2) # 计算中位生存时间
surv <- Survival(f2) # 构建生存概率函数
## 绘制COX回归中位生存时间的Nomogram图
nom <- nomogram(f2, fun=function(x) med(lp=x),funlabel="Median Survival Time")
plot(nom)

###绘制COX回归生存概率的Nomogram图
## LIHC数据的time是以”天“为单位,此处绘制1年，5年的生存概率
nom <- nomogram(f2, fun=list(function(x) surv(365, x),
                             function(x) surv(1825, x),
                             function(x) med(lp=x)),
                funlabel=c("1-year Survival Probability", "5-year Survival Probability","Median Survival Time"))
plot(nom, xfrac=.2)
