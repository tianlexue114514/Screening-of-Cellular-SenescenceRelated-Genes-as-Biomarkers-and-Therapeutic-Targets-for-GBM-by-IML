#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#���ð�
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

gene="TGFBI"      #���������(�����о���Ŀ���������޸�)
expFile="merge.normalize.txt"          #���������ļ�
gmtFile="c2.cp.kegg.Hs.symbols.gmt"    #�����ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\37GSEA\\T\\KEGG")      #���ù���Ŀ¼

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=meanH-meanL
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
#����pvalue<0.05�����ݽ��й���,�õ����������Ľ��
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#���Ƹ߱����鸻����ͼ��
termNum=5     #����չʾͨ·����Ŀ
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in high expression group")
	pdf(file="GSEA.highExp.pdf", width=6.5, height=5.5)
	print(gseaplot)
	dev.off()
}

#���Ƶͱ����鸻����ͼ��
termNum=5     #����չʾͨ·����Ŀ
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low expression group")
	pdf(file="GSEA.lowExp.pdf", width=6.5, height=5.5)
	print(gseaplot)
	dev.off()
}


