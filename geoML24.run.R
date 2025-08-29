#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="merge.normalize.txt"      #�����ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\28CIBERSORT")      #���ù���Ŀ¼
source("geoML24.CIBERSORT.R")       #���ð�

#����ϸ���������
outTab=CIBERSORT("ref.txt", inputFile, perm=1000)

#�����߽��������ˣ����ұ�������ϸ��������
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)




