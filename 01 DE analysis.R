#前期准备
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("affyPLM")
biocLite("RColorBrewer")
biocLite("impute")
biocLite("limma")
biocLite("pheatmap")
biocLite("ggplot2")


setwd("I:\\Pre_data\\\GEO\\...GSE_RAW\\Control")  #设置工作目录


library(affyPLM)
Data<-ReadAffy()
Pset<-fitPLM (Data)
image(Data[,1])###  image(Data[,2])
image(Pset,type="weights",which=1,main="Weights")
image(Pset,type="resids",which=1,main="Residuals")
image(Pset,type="sign.resids",which=1,main="Residuals.sign")
# 芯片质量  箱式图
library(affyPLM)
library(RColorBrewer)
Pset<-fitPLM (Data)
colors<-brewer.pal(12,"Set3")
Mbox(Pset,col=colors,main="RLE",las=3) 
Mbox(Pset,ylim=c(-1,1),col=colors,main="RLE",las=3)

library(affyPLM)
library(RColorBrewer)
Pset<-fitPLM (Data)
colors<-brewer.pal(12,"Set3")
boxplot(Pset,col=colors,main="NUSE",las=3)


#RNA降解
library(affy)
Data<-ReadAffy()
data.deg<-AffyRNAdeg(Data)
plotAffyRNAdeg(data.deg,col=colors)
legend("topleft",sampleNames(Data),col=colors,lwd=1,inset=0.5,cex=1)

#RMA背景化处理
setwd("") #设置工作目录
library(affyPLM)
library(affy)
Data<-ReadAffy()
Pset<-fitPLM (Data)
sampleNames(Data)
N=length(Data)
eset.rma<-rma(Data)
#读取探针值
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="CAD_exprs probeid.txt",sep='\t',quote=F,row.names=F)

#RMA
setwd("....")
library(affyPLM)
library(affy)
Data<-ReadAffy()
sampleNames(Data)
N=length(Data)
eset.rma<-rma(Data)
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="C_exprs.txt",sep='\t',quote=F,row.names=F)


#读取两个文件
#setwd(" ")  #设置工作目录
normal_exprs<-read.table("GPL85.txt",header=T,sep="\t")
tumor_exprs<-read.table("GSE738_Trans.txt",header=T,sep="\t")

#合并
probe_exprs<-merge(normal_exprs,tumor_exprs,by="probeid")
write.table(probe_exprs,file="NCvsET1_probeid.txt",sep='\t',quote=F,row.names=F)
###########################
###############################lincRNA###################
#####################################################
#Probe ID转Gene symbol
probe_exp<-read.csv("GSE738_Trans.txt",header=T,sep="\t",row.names=1)
probe_exp<- as.matrix(probe_exp)
#读取注释文件
probeid_geneid<-read.csv("GPL85.txt",header=T,sep="\t")
probe_name<-rownames(probe_exp)
#probe匹配
loc<-match(probeid_geneid[,1],probe_name)
#提取表达值
probe_exp<-probe_exp[loc,]
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]

exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean)) 

rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="exprs_exprs.txt",sep='\t',quote=F,row.names=F)
#gene id 转为gene symbol
loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<-rownames(gene_exp_matrix)
gene_exp_matrix3<-cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file="AB_all_genesybmol.txt",sep='\t',quote=F,row.names=F)

#处理缺失值 
#k-Nearest Neighbor
#impute
library(impute)
#取值
gene_exp_matrix<-read.table("Control vs14D.txt",header=T,sep="\t")
gene_exp_matrix <-gene_exp_matrix[!duplicated(gene_exp_matrix[,1]),] 
rownames(gene_exp_matrix) <- gene_exp_matrix[,1]
gene_exp_matrix <- gene_exp_matrix[,-1]
gene_exp_matrix<-as.matrix(gene_exp_matrix)
#############
#gene_exp_matrix <- read.table('GSE33970_WB_genesyb.exprs.txt',header=T,sep='\t') 
#gene_exp_matrix <-gene_exp_matrix[!duplicated(gene_exp_matrix[,1]),]
#rownames(gene_exp_matrix) <- gene_exp_matrix[,1] 
#gene_exp_matrix <- gene_exp_matrix[,-1]

#KNN缺失值
imputed_gene_exp<-impute.knn(gene_exp_matrix,k=10,rowmax=0.5,colmax=0.8,maxp=3000,rng.seed=362436069)
#
GeneExp<-imputed_gene_exp$data
#写
genesymbol<-rownames(GeneExp)
GeneExp<-cbind(genesymbol,GeneExp)
write.table(GeneExp,file="Control vs14D_KNN.txt",sep='\t',quote=F,row.names=F)
#######################差异##########################################
library(limma)
rt<-read.table("Control vs14D_KNN.txt",header=T,sep="\t",row.names=1)
#differential
class<-c(rep("control",3),rep("Is",3))
design<-model.matrix(~factor(class))
colnames(design)<-c("control","Is")
fit<-lmFit(rt,design)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffLab<-allDiff[with(allDiff, ((logFC>1.5 |logFC<(-1.5)) & adj.P.Val<0.05)),]
write.table(diffLab,file="diffEXp_1.5.xls",sep="\t",quote=F)
#读取每个样本的DEGs表达水平
diffExpLevel<-rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)
write.table(diffExpQvalue,file="diffExpLeve_1.5.xls",sep="\t",quote=F)

#上下调基因提取
Upgene = allDiff[(allDiff$adj.P.Val < 0.05 & (allDiff$logFC>1.5)),]
write.table(Upgene, "Upgene1.5.xls",sep="\t",quote=F)
#
Downgene = allDiff[(allDiff$adj.P.Val < 0.05 & (allDiff$logFC<(-1.5))),]
write.table(Downgene, "Downgene1.5.xls",sep="\t",quote=F)

