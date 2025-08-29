#install.packages("ggplot2")
#install.packages("ggrepel")


#���ð�
library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter=0.585              #logFC��������
adj.P.Val.Filter=0.05          #�������pֵ��������
diffFile="all.txt"             #���л����������Ľ���ļ�
geneFile="model.genes.txt"     #ģ�ͻ�����ļ�
method="Stepglm[backward]+Lasso"           #ѡ�����ѧϰ�ķ���(��Ҫ������ͼ�����޸�)
setwd("G:\\Users\\x1525\\Desktop\\new\\18vol")       #���ù���Ŀ¼

#��ȡ��������Ľ���ļ�
rt=read.table(diffFile, header=T, sep="\t", check.names=F)
row.names(rt)=rt[,1]
#����������
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

#���ƻ�ɽͼ
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "grey","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))

#��ͼ���б�ע�������������
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
geneRT=geneRT[geneRT$algorithm==method,]
sameGene=intersect(as.vector(geneRT[,1]), row.names(rt))
showData=rt[sameGene,]
p1=p+geom_label_repel(data=showData,
                    box.padding=0.2, point.padding=0.2, min.segment.length=0.1,
                    size=3, aes(label=id)) + theme_bw()
#���ͼ��
pdf(file="vol.pdf", width=5.25, height=4.5)
print(p1)
dev.off()

#�������ģ�ͻ�����б�
write.table(sameGene, file="modelGene.list.txt", sep="\t", quote=F, row.names=F, col.names=F)
#�������ģ�ͻ���Ĳ������
write.table(showData, file="modelGene.diff.txt", sep="\t", quote=F, row.names=F)



