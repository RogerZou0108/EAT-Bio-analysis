#���þ���װ��
#options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
#install.packages("ClassDiscovery")
#�ⲿ����
library(ClassDiscovery) 
library(pheatmap) #������ͼ
#��ͼ��ɫ
library(gplots) 
Sys.setenv(LANGUAGE = "en") #��ʾӢ�ı�����Ϣ
options(stringsAsFactors = FALSE) #��ֹchrת��factor
mat <- read.csv("lncRNA_exprs.csv",row.names = 1,check.names = F,stringsAsFactors = F)
mat[1:3,1:3]
annCol <- read.csv("annotation.csv",row.names = 1,check.names = F,stringsAsFactors = F)
head(annCol)
# �ȱʧֵ
annCol[is.na(annCol) | annCol == ""] <- "N/A"
# Ϊannotation.csv��ÿһ����ɫ
annColors <- list()
annColors[["gender"]] <- c("MALE"="blue","FEMALE"="red","N/A"="white")
annColors[["vital_status"]] <- c("Alive"="yellow","Dead"="black","N/A"="white")
annColors[["colon_polyps_present"]] <- c("YES"="red","NO"="black","N/A"="white")
annColors

#ֱ����pheatmap���ಢ��ͼ�����νṹ�۲����������ķ���
pheatmap(mat,
         scale = "row",
         color = bluered(64),
         annotation_col = annCol,
         annotation_colors = annColors,
         show_rownames = T, show_colnames = F)

pheatmap(mat,
         color=colorRampPalette(c("white","blue"))(100),
         annotation_col = annCol,
         annotation_colors = annColors,
         show_rownames = T, show_colnames = F)

#�ⲿ���������
# distanceMatrix()�Լ�hclust()
#�˽�distance��Ⱥ�linkage����
hcs <- hclust(distanceMatrix(as.matrix(mat), "pearson"), "ward.D") 
hcg <- hclust(distanceMatrix(t(as.matrix(mat)), "pearson"), "ward.D") # ע����뺯��������еģ����Զ��о���Ҫת��
group <- cutree(hcs,k=3)

# ����annotation����ɫ
annCol$Clust2 <- paste0("C",group[rownames(annCol)])
annColors[["Clust2"]] <- c("C1"="#fedcbd","C2"="#7f7522","C3"="#2b4490")
head(annCol)
pheatmap(mat,
         scale = "row",
         color = bluered(64),
         cluster_rows = hcg,
         cluster_cols = hcs,
         annotation_col = annCol,
         annotation_colors = annColors,
         show_rownames = T,show_colnames = F,
         filename = "heatmap_with_outside_Cluster1.pdf")

#��������ļ�
sample_order <- data.frame(row.names = seq(1:length(hcs$labels)), sample = hcs$labels, group = group)
sample_order <- sample_order[hcs$order,] #����������
sample_order$ori.order <- row.names(sample_order) #���������˳��
write.csv(sample_order, "sample_order.csv", quote = F, row.names = F)

#�ⲿ�������������ͼ�в���ʾ���ṹ��
pheatmap(mat,
         scale = "row",
         color = bluered(64),
         cluster_rows = hcg,
         cluster_cols = hcs,
         treeheight_col = 0, # �е�����Ϊ0����ͼ�н�����ʾ���ṹ������ʱ�������Ǳ������Ƶ�
         annotation_col = annCol,
         annotation_colors = annColors,
         show_rownames = T,show_colnames = F,
         filename = "heatmap_with_outside_Cluster_noTree.pdf")
#�������ṹ, ������noTree��ͼ��һ�µ�
index <- order.dendrogram(as.dendrogram(hcs))
sam_order <- colnames(mat)[index] 

pheatmap(mat[,sam_order], #ע���������Ҫ����˳��
         scale = "row",
         color = bluered(64),
         cluster_cols = F,
         cluster_rows = hcg,
         #treeheight_col = 0, #������佫û���κ�Ч��
         annotation_col = annCol[sam_order,], #ע���ʱע���ļ�ҲҪ�޸�˳��ҲҪ�޸�˳��ҲҪ�޸�˳�򣡣���
         annotation_colors = annColors,
         show_rownames = T, show_colnames = F,
         filename = "heatmap_with_outside_Cluster_discardTree.pdf")
