#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#���ð�
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      #pֵ��������
adjPvalFilter=0.05     #�������pֵ��������

#����ͼ�ε���ɫ
colorSel="p.adjust"
if(adjPvalFilter>0.05){
	colorSel="pvalue"
}

drugFile="DSigDB_All_detailed.txt"      #ҩ�����ݿ��ļ�
hubFile="hubGenes.txt"                  #���Ļ�����ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\32drugEnrich")      #���ù���Ŀ¼

#��ȡ���Ļ�����б��ļ�
rt=read.table(hubFile, header=T, sep="\t", check.names=F)
#��ȡ���Ļ��������
genes=unique(as.vector(rt[,1]))

#��ȡҩ�����ݿ��ļ�
drugRT=read.table(drugFile, header=T, sep="\t", check.names=F, quote="", comment.char="")
drugRT=drugRT[,1:2]

#ҩ�︻������
kk=enricher(genes,
	pvalueCutoff=1, qvalueCutoff=1,
	minGSSize = 10, maxGSSize = 500,
	TERM2GENE=drugRT)

#������������Ľ��
DRUG=as.data.frame(kk)
DRUG=DRUG[(DRUG$pvalue<pvalueFilter & DRUG$p.adjust<adjPvalFilter),]
write.table(DRUG[,-c(3,4)], file="DRUG.enrich.xls", sep="\t", quote=F, row.names = F)

#����չʾҩ�����Ŀ
showNum=30
if(nrow(DRUG)<showNum){
	showNum=nrow(DRUG)
}

#��״ͼ
pdf(file="barplot.pdf", width=7, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#����ͼ
pdf(file="bubble.pdf", width=7, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

#�����ҩ��Ĺ�ϵͼ
pdf(file="cnetplot.pdf", width=9, height=7)
cnet=cnetplot(kk, circular=TRUE, showCategory=10, colorEdge=TRUE)
print(cnet)
dev.off()


