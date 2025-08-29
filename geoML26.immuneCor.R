#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#���ð�
library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)

gene="TGFBI"                         #���������(�����о���Ŀ���������޸�)
expFile="merge.normalize.txt"       #���������ļ�
immFile="CIBERSORT-Results.txt"     #����ϸ���������ļ�
setwd("G:\\Users\\x1525\\Desktop\\new\\IMMU\\26immuneCor\\TGFBI")     #���ù���Ŀ¼

#��ȡ���������ļ�,�������ݽ��д���
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#ȥ�����������Ʒ
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

#��ȡĿ�����ı�����
data=t(data[gene,,drop=F])
data=as.data.frame(data)

#��ȡ����ϸ������Ľ���ļ�
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])

#������ϸ������ѭ�������������ɢ��ͼ
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-1)]){
	x=as.numeric(rt[,gene])
	y=as.numeric(rt[,i])
	if(sd(y)==0){y[1]=0.00001}
	cor=cor.test(x, y, method="spearman")
	
	outVector=cbind(Gene=gene, Cell=i, cor=cor$estimate, pvalue=cor$p.value)
	outTab=rbind(outTab,outVector)
	
	#����������������ϸ�����п��ӻ�, ���������ͼ��
	if(cor$p.value<0.05){
		outFile=paste0("cor.", i, ".pdf")
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab(paste0(gene, " expression")) + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="density", xparams=list(fill = "#FD7446FF"), yparams=list(fill = "#3B4992FF"))
		#�����ͼ��
		pdf(file=outFile, width=5, height=4.3)
		print(p2)
		dev.off()
	}
}
#�������ϸ������Եı����ļ�
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)


#####################��������԰�����ͼ#####################
#��ȡ����ԵĽ��
data = read.table("cor.result.txt", header=T, sep="\t", check.names=F)

#����ԲȦ��ɫ�ĺ���
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                ifelse(x>0.2,p.col[4], p.col[5])
                )))
  return(color)
}

#��������ԲȦ��С�ĺ���
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
              ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

#��������Լ����pvalue����ԲȦ��ɫ
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color

#�������ϵ������ԲȦ��С
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
data=data[order(data$cor),]

########����ͼ��########
pdf(file="Lollipop.pdf", width=9.5, height=7)      #���ͼ��
xlim = ceiling(max(abs(data$cor))*10)/10           #����x�᷶Χ
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2)
#��ͼ���������߶�
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
#��ͼ��������ԲȦ
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
#��ͼ����չʾ����ϸ��������
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5)
#չʾ����Լ����pvalue
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

#����ԲȦ��С��ͼ��
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#����ԲȦ��ɫ��ͼ��
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()


