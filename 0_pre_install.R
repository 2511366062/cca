options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
#options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))#清华镜像
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F,version = "3.21")
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

cran_packages <- c('tidyverse',
                   'msigdbr',
                   'patchwork',
                   'SeuratObject',
                   'Seurat'
                   ) 
Biocductor_packages <- c('sva',
                         'monocle',
                         'GOplot',
                         'GSVA',
                         'plotmo',
                         'regplot',
                         'scRNAseq',
                         'BiocStyle',
                         'celldex',
                         'SingleR',
                         'BiocParallel'
)

for (pkg in cran_packages){
  if (! require(pkg,character.only=T,quietly = T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T,quietly = T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}

.libPaths(c("/sharelib/x86_64-pc-linux-gnu-library/4.3","/usr/local/lib/R/site-library",
 "/usr/lib/R/site-library","/usr/lib/R/library" ,.libPaths()))
setwd("/home/liniannian/lxk/workspace/scrna/cca")

packageVersion("Seurat")

install.packages("wesanderson")
source("MyFunctions.R")
library(scRNAtoolVis)###umap图、热图等工具
library(paletteer)##
library(RColorBrewer)#扩展颜色

#富集分析
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
#library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(harmony)
# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)  # 人类数据，其他物种需更换（如小鼠：org.Mm.eg.db）
library(DOSE)
library(enrichplot)

library(RColorBrewer)
#rm(list = ls()[! ls() %in% c("sce.all")])
source("MyFunctions.R")



library(ggprism)##
library(ggrepel)
library(tidyverse)


source("function/bar.R")
source("function/fuji.R")
source("function/genes_compare.R")
source("function/heatmap.R")




###

