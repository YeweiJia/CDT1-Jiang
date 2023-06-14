setwd("/Users/jiayewei/Documents/Bioinformation/Lung cancer/LUAD/1. KRT80 in LUAD/6. Cibersort/Heatmap")
rm(list=ls())

library(pheatmap)
library(vegan)

data<-read.table('LUAD heatmap.txt',head=T,sep="\t",row.names = 1)
head(data)
data<-as.matrix(data)

#更改热图颜色
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
#color.key 是构建的颜色变化尺，只需要改变color.key


pheatmap(data,scale = "row") #参数归一化
pheatmap(data,treeheight_row = 50,treeheight_col = 20) #设置col,row方向的树高
pheatmap(data,cellwidth = 15,cellheight = 12,main="Example heatmap") #设定每一个热图格子的宽度和高度，main参数添加主标题
pheatmap(data,cluster_cols = FALSE, color=colorRampPalette(c("navy","white","firebrick3"))(50)) #取消聚类，更改颜色
pheatmap(data,show_rownames = F,show_colnames = F)#设置不显示行名和列名
pheatmap(data,display_numbers = TRUE,number_color = "dark")#在每个热图格子里显示相应数值

#注释添加行名信息
annotation_col=data.frame(Group=factor(rep(c('Low',"High"),c(263,264)))) #注释，按照分组，4Normal组，4Necrosis组
rownames(annotation_col)=colnames(data)
pheatmap(data,scale = "row",show_colnames = F,color=colorRampPalette(color.key)(50),border_color = NA,
         annotation_col = annotation_col)
#注释2 不同分组 添加行名信息
annotation_col=data.frame(Group=factor(rep(c("Normal","Necrosis"),4))) #注释，按照分组，交叉分组



#注释添加列名信息
annotation_row=data.frame(Genes=rep(c("up","down"),c(19,20)))
rownames(annotation_row)=rownames(data)
pheatmap(data,color=colorRampPalette(color.key)(50),border_color = NA,
         annotation_row = annotation_row)
