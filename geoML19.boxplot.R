#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("PerformanceAnalytics")


#���ð�
library(limma)
library(reshape2)
library(ggpubr)
library(PerformanceAnalytics)

expFile="merge.normalize.txt"      #���������ļ�
geneFile="modelGene.list.txt"      #�����б��ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\19boxplot")     #���ù���Ŀ¼

#��ȡ���������ļ�
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#��ȡ�����б��ļ�, ��ȡģ�ͻ���ı�����
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

#��ȡ��Ʒ�ķ�����Ϣ(�������ʵ����)
Type=gsub("(.*)\\_(.*?)", "\\2", row.names(data))
treatData=data[Type=="Treat",]
rt=cbind(as.data.frame(data), Type)

#������ת��Ϊ����ͼ�������ļ�
data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#��������ͼ
p=ggboxplot(data, x="Gene", y="Expression", fill = "Type",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Type",
	     palette = c("#0088FF", "#FF5555"), width=0.75)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#���ͼ��
pdf(file="boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()

#���������ͼ��(PerformanceAnalytics��)
pdf(file="cor.pdf", width=7, height=6.5)
chart.Correlation(treatData, histogram=TRUE, pch=19, method="pearson")
dev.off()


