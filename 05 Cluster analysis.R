#设置镜像安装包
#options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
#install.packages("ClassDiscovery")
#外部聚类
library(ClassDiscovery) 
library(pheatmap) #绘制热图
#热图颜色
library(gplots) 
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
mat <- read.csv("lncRNA_exprs.csv",row.names = 1,check.names = F,stringsAsFactors = F)
mat[1:3,1:3]
annCol <- read.csv("annotation.csv",row.names = 1,check.names = F,stringsAsFactors = F)
head(annCol)
# 填补缺失值
annCol[is.na(annCol) | annCol == ""] <- "N/A"
# 为annotation.csv的每一列配色
annColors <- list()
annColors[["gender"]] <- c("MALE"="blue","FEMALE"="red","N/A"="white")
annColors[["vital_status"]] <- c("Alive"="yellow","Dead"="black","N/A"="white")
annColors[["colon_polyps_present"]] <- c("YES"="red","NO"="black","N/A"="white")
annColors

#直接用pheatmap聚类并画图：树形结构观察样本基本的分类
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

#外部聚类产生树
# distanceMatrix()以及hclust()
#了解distance测度和linkage方法
hcs <- hclust(distanceMatrix(as.matrix(mat), "pearson"), "ward.D") 
hcg <- hclust(distanceMatrix(t(as.matrix(mat)), "pearson"), "ward.D") # 注意距离函数是针对列的，所以对行聚类要转置
group <- cutree(hcs,k=3)

# 增加annotation及配色
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

#保存聚类文件
sample_order <- data.frame(row.names = seq(1:length(hcs$labels)), sample = hcs$labels, group = group)
sample_order <- sample_order[hcs$order,] #按聚类排序
sample_order$ori.order <- row.names(sample_order) #保存最初的顺序
write.csv(sample_order, "sample_order.csv", quote = F, row.names = F)

#外部聚类产生树（热图中不显示树结构）
pheatmap(mat,
         scale = "row",
         color = bluered(64),
         cluster_rows = hcg,
         cluster_cols = hcs,
         treeheight_col = 0, # 列的树高为0，热图中将不显示树结构，但此时样本还是被树限制的
         annotation_col = annCol,
         annotation_colors = annColors,
         show_rownames = T,show_colnames = F,
         filename = "heatmap_with_outside_Cluster_noTree.pdf")
#丢弃树结构, 和上面noTree的图是一致的
index <- order.dendrogram(as.dendrogram(hcs))
sam_order <- colnames(mat)[index] 

pheatmap(mat[,sam_order], #注意输入矩阵要更改顺序
         scale = "row",
         color = bluered(64),
         cluster_cols = F,
         cluster_rows = hcg,
         #treeheight_col = 0, #这行语句将没有任何效果
         annotation_col = annCol[sam_order,], #注意此时注释文件也要修改顺序，也要修改顺序，也要修改顺序！！！
         annotation_colors = annColors,
         show_rownames = T, show_colnames = F,
         filename = "heatmap_with_outside_Cluster_discardTree.pdf")

