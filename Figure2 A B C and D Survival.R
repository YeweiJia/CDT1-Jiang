
setwd("/Users/jiayewei/Desktop/GSE46517 Necrosis/Necrosis verus Normal Survival")
rm(list=ls())
options(stringsAsFactors = F)
# install.packages("survminer")
library(survminer)
library(survival)
a=read.table('',header = T,sep = '\t')
head(a)
fit<-survfit(Surv(OS.time,OS)~gender,data=a);fit
ggsurvplot(fit,data=a,surv.median.line = "hv",conf.int = TRUE)

fit2<-survfit(Surv(OS.time,OS)~group,data=a);fit
ggsurvplot(fit2,data=a,surv.median.line = "hv",conf.int = TRUE)


##绘制累积风险曲线
ggsurvplot(fit2, data=a, conf.int = TRUE,#增加置信区间
           fun="cumhaz")##绘制累计风险曲线

####添加风险表

ggsurvplot(fit2,data=a, conf.int = TRUE, risk.table=TRUE)

ggsurvplot(fit2,
           dat=a, 
           conf.int = TRUE,
           pval=TRUE,##添加值
           surv.median.line = "hv",##添加中位生存时间线
           risk.table=TRUE,##添加风险表
           ) ##设置X轴刻度间距
