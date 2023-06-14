setwd('/Users/jiayewei/Desktop')

rm(list = ls())
library(ggplot2)

dataset <- read.table(file = "", 
                      header = TRUE, sep = "",row.names = 1)
# 设置p_value和logFC的阈值
cut_off_pvalue = 0.05 #统计显著性
cut_off_logFC = 1.5       #差异倍数值

# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
dataset$change = ifelse(dataset$adj.P.Val< cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, 
                        ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(dataset)

geneList0 <- c('CDT1')
geneList <- dataset[geneList0,]
p <- ggplot(# 数据、映射、颜色
  dataset, aes(x = logFC, y = -log10(adj.P.Val), colour=change)) +
  geom_point(alpha=0.5, size=3.5) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  #突出表示差异基因
  geom_point(data=geneList,aes(x = logFC, y = -log10(adj.P.Val)),colour="yellow",size=3.5)+
  #辅助线
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)")+   # 坐标轴# 坐标轴和图标题title="Volcano plot",
  theme_bw()+    #去除背景色
  theme(panel.grid = element_blank())+  #去除网格线
  #xlim(-2, 2)+   #设置坐标轴范围
  #图例
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="bottom", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
        axis.title.x =element_text(size=18), 
        axis.title.y=element_text(size=18),
        axis.text=element_text(size=14,face = "bold"))
p

dataset$gene <- rownames(dataset)
rownames(dataset) <- 1:nrow(dataset)

#标记出5个基因的label
geneList <- as.data.frame(geneList0)
geneList[,2] <- geneList
colnames(geneList) <- c('gene','label')


c <-merge(dataset,geneList,by='gene',all.x=TRUE)  #增加label列，以突出显示指定基因


library(ggrepel)
p + geom_label_repel(data = c, 
                     aes(x = logFC, y = -log10(adj.P.Val), label = label),
                     size = 4,color="black",
                     #box.padding = unit(0.5, "lines"),
                     #point.padding = unit(0.8, "lines"), 
                     #segment.color = "black",   #连线的颜色
                     #segment.size = 0.4,  #连线的粗细
                     #arrow = arrow(length=unit(0.01, "npc")), #标签、点之间连线的箭头
                     show.legend = FALSE)
