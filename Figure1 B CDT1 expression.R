setwd("/Users/jiayewei/Desktop/For Yuanpang/DLGAP5 expression")
rm(list=ls())

df=read.table('gene',header = T,sep = '\t')
head(df)
library(ggplot2)
library(ggsignif)

#### t.test() 比较两组（参数) wilcox.test()比较两组（非参数）
#### aov()或anova()比较多组（参数）kruskal.test()比较多组

###them设置
###all line elements (element_line())； all rectangular elements (element_rect())；
###all text elements (element_text())；aspect.ratio：aspect ratio of the panel

compaired<-list(c("Primary Tumor","Solid Tissue Normal"))
p<-ggplot(data=df,aes(x=group, y=gene))+
  theme_bw()+ geom_boxplot(aes(fill=group))+ geom_signif(comparisons=compaired,step_increase = 0.1,map_signif_level = F,test=t.test)
p+theme(text=element_text(colour = "black",size=12),
        axis.text.x = element_text(face = "bold",color = "black",size=12),
        axis.title.y = element_text(color="black",hjust=0.5,size=12),
        aspect.ratio = 1.5)
###小提琴图
e <- ggplot(df, aes(x = group, y = Neutrophils))
#基础
e + geom_violin()
#旋转
e + geom_violin() + coord_flip()
#不修剪小提琴的尾部
e + geom_violin(trim = FALSE, fill = "steelblue")

###stat_summary() 可以加均值/中位值等

###用fun= mean/median加均值/中位值

e + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
###加均值及标准差
e + geom_violin(trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")


###改变颜色
e + geom_violin(aes(color = group), trim = FALSE)
e + geom_violin(aes(fill = group), trim = FALSE)
