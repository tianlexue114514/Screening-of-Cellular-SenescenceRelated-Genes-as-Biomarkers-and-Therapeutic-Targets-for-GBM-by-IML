#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#���ð�
library(limma)
library(sva)

geneFile="interGenes.txt"      #�����б��ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\14MLdata")      #���ù���Ŀ¼

#��ȡĿ¼������"normalize.txt"��β���ļ�
files=dir()
files=grep("normalize.txt$", files, value=T)
geneList=list()

#��ȡ���б��������ļ��еĻ�����Ϣ�����浽geneList
for(file in files){
    rt=read.table(file, header=T, sep="\t", check.names=F)      #��ȡ�����ļ�
    geneNames=as.vector(rt[,1])      #��ȡ��������
    uniqGene=unique(geneNames)       #����ȡunique
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

#��ȡ�������ݵĽ�������
interGenes=Reduce(intersect, geneList)

#���ݺϲ�
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    #��ȡ�����ļ������������ļ���������
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))

    #���ݺϲ�
    if(i==1){
    	allTab=rt[interGenes,]
    }else{
    	allTab=cbind(allTab, rt[interGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#��ȡ��������ı�����
svaTab=ComBat(allTab, batchType, par.prior=TRUE)
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
geneTab=svaTab[intersect(row.names(svaTab), as.vector(geneRT[,1])),]
geneTab=t(geneTab)

#��ȡѵ�������֤�������
train=grepl("^merge", rownames(geneTab), ignore.case=T)
trainExp=geneTab[train,,drop=F]
testExp=geneTab[!train,,drop=F]
rownames(trainExp)=gsub("merge_", "Train.", rownames(trainExp))
trainType=gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(trainExp))
testType=gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(testExp))
trainType=ifelse(trainType=="Control", 0, 1)
testType=ifelse(testType=="Control", 0, 1)
trainExp=cbind(trainExp, Type=trainType)
testExp=cbind(testExp, Type=testType)

#���ѵ�������֤�������
trainOut=rbind(id=colnames(trainExp), trainExp)
write.table(trainOut, file="data.train.txt", sep="\t", quote=F, col.names=F)
testOut=rbind(id=colnames(testExp), testExp)
write.table(testOut, file="data.test.txt", sep="\t", quote=F, col.names=F)



