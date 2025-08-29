#install.packages("glmnet")
#install.packages("pROC")


#���ð�
library(glmnet)
library(pROC)

expFile="merge.normalize.txt"      #���������ļ�
geneFile="modelGene.list.txt"      #�����б��ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\20geneROC")    #���ù���Ŀ¼

#��ȡ���������ļ�
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#��ȡ��Ʒ�ķ�����Ϣ(�������ʵ����)
y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(rt))
y=ifelse(y=="Control", 0, 1)

#��ȡ�����б��ļ�
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#����ͼ�ε���ɫ
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)

#�Ի������ѭ��������ROC����
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
	k=k+1
	#����ROC����
	roc1=roc(y, as.numeric(rt[x,]))     #�õ�ROC���ߵĲ���
	if(k==1){
		pdf(file="ROC.genes.pdf", width=5, height=4.5)
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}else{
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3, add=TRUE)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}
}
#����ͼ�����õ�ROC�����µ����
legend("bottomright", aucText, lwd=3, bty="n", cex=0.8, col=bioCol[1:(ncol(rt)-1)])
dev.off()




