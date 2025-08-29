#install.packages("reshape2")
#install.packages("ggplot2")


#���ð�
library(reshape2)
library(ggplot2)

setwd("G:\\Users\\x1525\\Desktop\\new\\07boxplot")     #���ù���Ŀ¼

#��������ͼ�ĺ���
bioBoxplot=function(inputFile=null, outFile=null, titleName=null){
	#��ȡ�����ļ�,��ȡ����
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))            #��ȡ�����о���id
	Sample=gsub("(.+)\\_(.+)\\_(.+)", "\\2", rownames(data))     #��ȡ��Ʒ������
	data=cbind(as.data.frame(data), Sample, Project)
	
	#������ת����ggplot2�����ļ�
	rt1=melt(data, id.vars=c("Project", "Sample"))
	colnames(rt1)=c("Project","Sample","Gene","Expression")

	#��������ͼ
	pdf(file=outFile, width=10, height=5)
	p=ggplot(rt1, mapping=aes(x=Sample, y=Expression))+
  		geom_boxplot(aes(fill=Project), notch=T, outlier.shape=NA)+
  		ggtitle(titleName)+ theme_bw()+ theme(panel.grid=element_blank())+ 
  		theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5,size=2), plot.title=element_text(hjust = 0.5))
	print(p)
	dev.off()
}

#���ú���, �������ν���ǰ������ͼ
bioBoxplot(inputFile="merge.preNorm.txt", outFile="boxplot.preNorm.pdf", titleName="Before batch correction")
#���ú���, �������ν����������ͼ
bioBoxplot(inputFile="merge.normalize.txt", outFile="boxplot.normalzie.pdf", titleName="After batch correction")




