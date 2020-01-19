
#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.8")


setwd("F:/新生儿心外膜脂肪/GSE82155_RAW/mRNA_DEGs and Immune analysis/Immune analsysis/CIBERSORT")
source("GEOimmune.CIBERSORT.R")
results=CIBERSORT("ref.txt", "normalize.txt", perm=1000, QN=TRUE)
