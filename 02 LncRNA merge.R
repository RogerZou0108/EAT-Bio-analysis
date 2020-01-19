#############################################################
############设置目录#########################################
setwd("")

deg=read.csv('Enriched.txt',sep = '\t',header=T)
exp=read.csv('After Filterring Exp_F.txt',sep = '\t',header=T)
##提取
#exp = read.table("Post_Combat.genesymbol.txt",header=T,sep="\t")
deg.exp = exp[match(deg$ID,exp$ID),]
write.table(deg.exp,"After Filterring Exp_F_hubDEGs.txt",col.names = T,row.names = F,quote=F,sep="\t")

###################合并######################################
normal_exprs<-read.table("GSE32062_series_matrix.txt",header=T,sep="\t")
tumor_exprs<-read.table("GPL6480Anno.txt",header=T,sep="\t")

probe_exprs<-merge(normal_exprs,tumor_exprs,by="ID")
write.table(probe_exprs,file="Co-Gene.txt",sep='\t',quote=F,row.names=F)

