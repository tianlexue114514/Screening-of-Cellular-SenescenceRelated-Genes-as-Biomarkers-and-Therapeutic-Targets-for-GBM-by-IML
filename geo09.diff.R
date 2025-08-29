#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("ggplot2")


#���ð�
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)

logFCfilter=0.585          #logFC�Ĺ�������(logFC=0.585,���챶��1.5��;logFC=1,����2��;logFC=2,����4��)
adj.P.Val.Filter=0.05      #������pֵ�Ĺ�������
inputFile="merge.normalize.txt"      #���������ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\09diff")      #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#��ȡ��Ʒ�ķ�����Ϣ
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,order(Type)]       #������Ʒ�ķ�����Ϣ����Ʒ��������
Project=gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))    #���GEO�����о���id
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data)=gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

#�������
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("Control","Treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(Treat-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#������л���Ĳ������
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#��������Ĳ������
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#���������������
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#���Ʋ��������ͼ
geneNum=50     #����չʾ�������Ŀ
diffUp=diffSig[diffSig$logFC>0,]
diffDown=diffSig[diffSig$logFC<0,]
geneUp=row.names(diffUp)
geneDown=row.names(diffDown)
if(nrow(diffUp)>geneNum){geneUp=row.names(diffUp)[1:geneNum]}
if(nrow(diffDown)>geneNum){geneDown=row.names(diffDown)[1:geneNum]}
hmExp=data[c(geneUp,geneDown),]
#׼��ע���ļ�
names(Type)=colnames(data)
Type=as.data.frame(Type)
Type=cbind(Project, Type)
#���ͼ��
pdf(file="heatmap.pdf", width=10, height=7)
pheatmap(hmExp, 
         annotation_col=Type, 
         color = colorRampPalette(c("blue2", "white", "red2"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5.5,
         fontsize_col=8)
dev.off()


#����������
rt=read.table("all.txt", header=T, sep="\t", check.names=F)
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
#���ƻ�ɽͼ
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green2", "grey","red2"))+
    labs(title = " ")+ theme_bw()+ 
    theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))
#���ͼ��
pdf(file="vol.pdf", width=5.5, height=4.5)
print(p)
dev.off()



