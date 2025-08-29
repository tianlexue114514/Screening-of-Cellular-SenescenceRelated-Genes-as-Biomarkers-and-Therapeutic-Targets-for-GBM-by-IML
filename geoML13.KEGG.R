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
adjPvalFilter=1        #�������pֵ��������

#����ͼ�ε���ɫ
colorSel="p.adjust"
if(adjPvalFilter>0.05){
	colorSel="pvalue"
}

setwd("G:\\Users\\x1525\\Desktop\\new\\12KEGG")      #���ù���Ŀ¼
rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)     #��ȡ����������б��ļ�

#��ȡ���������, ����������ת��Ϊ����id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt, entrezIDs)
colnames(rt)[1]="gene"
rt=rt[rt[,"entrezIDs"]!="NA",]      #ȥ������idΪNA�Ļ���
gene=rt$entrezID
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg��������
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<adjPvalFilter),]
#������������Ľ��
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#����չʾͨ·����Ŀ
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#��״ͼ
pdf(file="barplot.pdf", width=8.5, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#����ͼ
pdf(file="bubble.pdf", width=8.5, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()


