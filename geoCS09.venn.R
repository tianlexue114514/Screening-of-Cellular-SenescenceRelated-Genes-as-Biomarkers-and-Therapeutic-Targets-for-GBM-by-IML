#install.packages("VennDiagram")


library(VennDiagram)      #���ð�
setwd("G:\\Users\\x1525\\Desktop\\new\\10diffCS")     #���ù���Ŀ¼
geneList=list()

#��ȡ��������Ľ���ļ�
diffRT=read.table("diff.txt", header=T, sep="\t", check.names=F, quote="", comment.char="")
geneNames=as.vector(diffRT[,1])          #��ȡ������������
geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
uniqGene=unique(geneNames)               #�Բ������ȡunique
geneList[["DEG"]]=uniqGene

#��ȡϸ��˥�ϻ�����б��ļ�
rt=read.table("cellage3.tsv", header=T, sep="\t", check.names=F, quote="", comment.char="")
geneNames=as.vector(rt[,2])              #��ȡϸ��˥�ϻ��������
geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
uniqGene=unique(geneNames)               #��ϸ��˥�ϻ���ȡunique
geneList[["CS"]]=uniqGene

#����vennͼ
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#�������������б�
interGenes=Reduce(intersect, geneList)
outTab=diffRT[diffRT[,1] %in% interGenes,]
write.table(outTab, file="CS.diff.txt", sep="\t", quote=F, row.names=F)




