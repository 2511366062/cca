




# 寻找好看配色
getClusterColor= function(n){
  
  colors <-c(
    # ---- 基础科学色（色盲安全）----
    "#4477AA",  # 深蓝 (Paul Tol) 1 
    "#EE7733",  # 橙 (Paul Tol) 2
    "#CC3311",  # 红 (Paul Tol) 3
    "#009988",  # 青 (Paul Tol) 4
    "#E69F00",  # 橙黄 (Okabe-Ito) 5
    "#56B4E9",  # 天蓝 (Okabe-Ito) 6
    "#009E73",  # 深绿 (Okabe-Ito) 7
    
    # ---- 自然/生物学色 ----
    "#8C510A",  # 棕 (Earth Tones) 8
    "#01665E",  # 深青绿 (Earth) 9
    "#023858",  # 深海蓝 (Ocean) 10
    "#3690C0",  # 中蓝 (Ocean) 11
    
    # ---- 现代简约色 ----
    "#90A4AE",  # 灰蓝 (Material Design) 12
    "#FFAB91",  # 浅橙 (Material) 13
    "#C5E1A5",  # 淡绿 (Material) 14 
    "#B39DDB",  # 淡紫 (Material) 15
    
    # ---- 高区分度扩展色 ----
    "#FBB4AE",  # 粉红 (Pastel) 16
    "#B3CDE3",  # 浅蓝 (Pastel) 17
    "#CCEBC5",  # 薄荷绿 (Pastel) 18
    "#FB8072",  # 橙红 (ColorBrewer Set3) 19
    "#80B1D3",  # 浅钢蓝 (Set3) 20
    "#BC80BD",  # 紫灰 (Set3) 21
    
    # ---- 冷/暖对比色 ----
    "#2B8CBE",  # 蓝 (冷色系) 22
    "#E34A33"   # 橙红 (暖色系) 23
  )
  
  colors[n]
  
  
  
}

#'@describeIn 修改差异分析p_val=0、p_val_adj=0为浮动的随机数 将0的行替换为10^-100. -  10^-300
#'@param p_values 从差异分析列表中直接获取：demo$p_val,demo$p_adj
#
#
replace_zero_with_random = function(p_values) {
  # 找到所有p=0的位置
  zero_indices <- which(p_values == 0)
  # 生成与零值数量相同的随机极小值（对数均匀分布）
  random_values <- 10^runif(length(zero_indices), min = -300, max = -100)
  p_values[zero_indices] <- random_values
  return(p_values)
}



# 修改p_val 、p_val_adj为0的项，使其为浮动的数
# 增加neg_log10_pval、neg_log10_padj、tag(down\up\no)列

diff_we = function(object,category,ident_1,ident_2,log2FC = 0.5){
  
  
  Idents(object) = category#可以直接写列名
  object = FindMarkers(object,ident.1 = ident_1,ident.2 = ident_2,min.pct = 0.25)
  
  object$p_val = replace_zero_with_random(object$p_val)
  object$p_val_adj = replace_zero_with_random(object$p_val_adj)
  
  object$neg_log10_pval = -log10(object$p_val)
  object$neg_log10_padj = -log10(object$p_val_adj)
  object$tag = "no"
  
  object$tag = ifelse(object$avg_log2FC > log2FC & object$p_val_adj < 0.05,"up",ifelse(object$avg_log2FC < -1 & object$p_val_adj <0.05,"down","no"))
  object$symbol = rownames(object)
  object
}

#获取上调的基因集，问AI哪些可以用
get_upgene = function(object){
  
  tagGene = rownames(object)[which(object$tag == "up" & object$avg_log2FC > 1 & object$p_val_adj <0.01)]
  tagGene
}

#通过传入的字符串向量，通过差异分析数据获取过滤掉的数据框
get_upgene_data = function(object,tagGene){
  
  tagGene_data = filter(object,rownames(object) %in% tagGene)
  tagGene_data
}


#画图火山图：

getMyHSplot = function(object,genes){
  data = filter(object,rownames(object) %in% genes)#数据框
  
  p1 = ggplot(object, aes(x = avg_log2FC, y = neg_log10_padj)) +
    geom_point(aes(color = tag), alpha = 0.5,size = 3) +  # 按显著性着色
    scale_color_manual(values = c("#0F542FFF","grey", "#950404FF")) +      # 自定义颜色
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +  # 添加log2FC阈值线
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p值阈值线
    labs(
      x = "log2(Fold Change)",
      y = "-log10(p_val_adj)",
      #title = "Volcano Plot"
    ) +
    theme_bw() + theme(panel.grid = element_blank())
  library(ggrepel)
  p1 = p1 + geom_label_repel(
    data = data,
    aes(label = rownames(data)),
    box.padding = 0.5,
    segment.color = "red",
    max.overlaps = 40) #显示的基因数限制
  
  
  p2 = ggplot(object, aes(x = avg_log2FC, y = neg_log10_padj)) +
    geom_point(aes(color = neg_log10_padj,
                   size = neg_log10_padj), 
               alpha = 0.8
    ) +  
    scale_color_gradientn(values = seq(0,1,0.2),
                          colors = c("#0F542FFF","#388F30FF","#9F5630FF","#950404FF","#E04B28FF"))+ 
    scale_size_continuous(range = c(2,4))+#点的大小
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +  # 添加log2FC阈值线
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p值阈值线
    labs(
      x = "log2(Fold Change)",
      y = "-log10(p_val_adj)",
      
    ) +
    theme_bw()+theme(panel.grid = element_blank())
  
  p2 = p2 + geom_label_repel(
    data = data,
    aes(label = rownames(data)),
    box.padding = 0.5,
    segment.color = "red",
    max.overlaps = 40 #显示的基因数限制
  )
  
  plots = list("p1" = p1,"p2" = p2)
  plots
  
  
}




#画KEGG、GO富集图
do_GO_KEGG = function(object,filego,filekegg,tag){
  
  up_gene = rownames(object)[which(object$tag == tag)]
  up_ENTREZID = bitr(up_gene, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
  up_ENTREZID = up_ENTREZID$ENTREZID
  
  ego <- enrichGO(gene          = up_ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "ALL",       # 可选 "BP", "MF", "CC"
                  pAdjustMethod = "BH",        # 多重检验校正方法
                  pvalueCutoff  = 0.05,        # p值阈值
                  qvalueCutoff  = 0.1,         # q值阈值
                  readable      = TRUE) 
  write.csv(ego, filego)
  kk <- enrichKEGG(gene         = up_ENTREZID,
                   organism     = "hsa",       # 人类：hsa，小鼠：mmu
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   use_internal_data = FALSE #FALSE为自动下载
  )
  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(kk, filekegg)  
  
  p1 = barplot(ego, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") #柱状图
  p2 = dotplot(ego,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
  p3 = barplot(kk,showCategory = 40,title = 'KEGG Pathway')
  p4 = dotplot(kk)
  
  plots = list("go_p1" = p1,
               "go_p2" = p2,
               "kegg_p1" = p3,
               "kegg_p2" = p4)
  
  
  
  
  
  
}




#获取细胞亚群总的normal-turn。以及多个样本的UMAP图

get_yaqun_umap = function(object){
  
  p1 = clusterCornerAxes(object = object,
                        reduction = 'umap',
                        noSplit = T,
                        pSize = 0.1,
                        arrowType='open',
                        clusterCol = 'category',
                        base_size = 25,
                        aspect.ratio = 1,#绘图的宽高比(正方形)，默认NULL
                        relLength = 0.5 #x轴y轴的相对某个的长度，0-1
  )+scale_color_manual(values = c("#6c77d1","#6ca9d1"))
  
  p2 = clusterCornerAxes(object = object,
                    reduction = 'umap',
                    noSplit = F,
                    pSize = 0.05,
                    arrowType='open',
                    groupFacet = 'sample',
                    clusterCol = 'seurat_annotation',
                    axes = 'one',
                    base_size = 15,
                    aspect.ratio = 1,#绘图的宽高比(正方形)，默认NULL
                    relLength = 0.5, #x轴y轴的相对某个的长度，0-1
  )+scale_color_paletteer_d("ggthemes::Red_Blue_Brown")
  
  plots = list("category" = p1,"sample" = p2)
  
  
}


# 多基因多分组的小提琴图
#'@param object serurat对象
#'@param genes 基因组，字符串向量
#'@param group 分组，用于对照
duosampes_duogenes_vlnPlot = function(object,genes,group){
  m = as.matrix(object@assays$RNA$data)#layers去掉不用加列明行名
  
  vln.df = m %>%
    t() %>%#转制
    as.data.frame()%>%
    dplyr::select(genes) %>% 
    rownames_to_column("CB") %>% 
    mutate(cluster = object@meta.data[,colnames(object@meta.data) == group])%>%#以什么分的页
    pivot_longer(cols = 2:(ncol(.)-1),
                 names_to = "gene",
                 values_to = "exp") %>% 
    mutate(gene = factor(gene,levels = genes))
  head(vln.df)
  
  library(paletteer)
  my_color = paletteer_d(`"ggsci::nrc_npg"`)
  my_color = colorRampPalette(my_color)(length(unique(vln.df$cluster)))
  # 画图
  p1 <- ggplot(vln.df,aes(exp,cluster),color=factor(cluster))+
    geom_violin(aes(fill=cluster),scale = "width")+
    scale_fill_manual(values = my_color)+
    facet_grid(.~gene,scales = "free_y", switch = "x")+
    scale_x_continuous(expand = c(0,0),position = "top")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title.x.top = element_blank(),
      #axis.ticks.x.bottom = element_blank(),
      axis.text.x.top= element_text(hjust = 1,vjust = NULL,color = "black",size = 7),
      #axis.title.y.left = element_blank(),
      #axis.ticks.y.left = element_blank(),
      #axis.text.y.left = element_blank(),
      legend.position = "none",
      panel.spacing.y = unit(0, "cm"),
      strip.text.y = element_text(angle=0,size = 14,hjust = 0),
      strip.background.y = element_blank()
    )
  p1
}


##根据meta.data的idents统计细胞个数
#'@param category_ 统计分组细胞个数
#'@param sample_ 统计细胞个数
numberCellsss = function(object,category_ = "",sample_ = ""){
  category_
  if(!category_ == ""){
    "das222"
    object = subset(object,category == category_)
  }
  
  if(!sample_ == ""){
    object = subset(object,sample == sample_)
  }
  
  #统计细胞个数
#Idents(object) = ident #修改统计方式
dat = as.data.frame(table(Idents(object)))
#lable自动类型
dat$label = paste(dat$Var1,dat$Freq,sep = ":")
head(dat)
library(ggplot2)
library(paletteer)
#View(palettes_d_names)

#ggthemes::Red_Blue_Brown
p = ggplot(dat,aes(x = Freq,fill = Var1,y = Var1))+
  scale_fill_paletteer_d("ggsci::nrc_npg")+
  geom_bar(stat = "identity")+
  theme_bw()+
  geom_text(aes(x = 0,label = label),hjust = 0)+
  theme(axis.text.y = element_blank(),   # 隐藏纵坐标刻度文字
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),# 隐藏纵坐标刻度线
        legend.title = element_blank(),#隐藏图例标题
    panel.grid.major = element_blank(),#隐藏背景网格
    panel.grid.minor = element_blank(),#隐藏背景网格
    legend.text = element_text(size = 18),
    #panel.border = element_blank(),
        )  
p
}



#' Title
#'
#' @param object seurat对象
#' @param gene 单个字符向量
#' @param category meta.data用于区分normal vs turned的列
#' @param cellname 细胞名，在哪个组内提供
#'
#' @returns
#' @export 声明该函数需要导出到包的命名空间
#'
#' @examples 提供可运行的示例代码
compare_OneGene = function(object,gene,category,cellname){
  
  library(ggpubr)            #引用包
  
  #设置比较租
  rt = FetchData(object, vars = c(gene, category))
  colnames(rt) = c("Expression","Type")
  rt = rt[rt$Expression > 0,] #要么是数字代表那个位置，要么是T/F
  x = cellname
  y = gene
  group=levels(factor(rt$Type))
  rt$Type=factor(rt$Type, levels=group)#确定分组顺序
  comp=combn(group,2)#生成组合
  my_comparisons = as.vector(group)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #绘制boxplot
  boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                    xlab=x,
                    ylab=y,
                    legend.title=x,#图例名称
                    palette = c("#e06c8b","#6ce0c0"),
                    #add = "jitter"
                    )+ 
    stat_compare_means(comparisons = my_comparisons)
  boxplot
}




fuji_pdf = function(){
  
  c(width = 12,height = 14)
  
}

#亚群聚类
#有celltype这一列
#resolution 分辨率,默认0.5
one_to_all = function(object,celltype,resolution = 0.5){
  
  obj = subset(object,celltype == celltype)
  obj = obj %>% 
    NormalizeData() %>%  
    FindVariableFeatures() %>%  
    ScaleData(features = rownames(.)) %>%  
    RunPCA(pc.genes = VariableFeatures(.))  %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(dims = 1:15, reduction = "harmony") %>% 
    FindClusters(resolution = resolution) %>% 
    RunUMAP(dims = 1:15,reduction = "harmony") #%>% 
  
  obj
}



#注释基因
zhushi =  function(object,file) {
  
  celltype = read.csv(file = file,sep=",",check.names =FALSE,header = F) #自己照着DotPlot图填的
  celltype
  new.cluster.ids <- celltype$V2
  names(new.cluster.ids) <- levels(object)
  object <- RenameIdents(object, new.cluster.ids)
  object$celltype = Idents(object)
  object
  
}

#绘制：
#两个组的UMAP
#多个样本的UMAP 
#细胞统计图：组间、样本间 
#比例图

yaqun_paint_umap = function(obj) {
  
  p1 = clusterCornerAxes(object = obj,
                        reduction = 'umap',
                        noSplit = F,
                        pSize = 0.1,
                        arrowType='open',
                        groupFacet = 'category',
                        clusterCol = 'celltype',
                        base_size = 25,
                        aspect.ratio = 1,#绘图的宽高比(正方形)，默认NULL
                        relLength = 0.5 #x轴y轴的相对某个的长度，0-1
  )+scale_color_manual(values = palette)
  
  p2 = clusterCornerAxes(object = obj,
                        reduction = 'umap',
                        noSplit = F,
                        pSize = 0.05,
                        arrowType='open',
                        groupFacet = 'sample',
                        clusterCol = 'celltype',
                        
                        themebg = "bwCorner",
                        base_size = 15,
                        aspect.ratio = 1,#绘图的宽高比(正方形)，默认NULL
                        relLength = 0.5, #x轴y轴的相对某个的长度，0-1
  )+scale_color_manual(values = palette)
  
  list("umap_comp" = p1,"umap_sample" = p2)

}


#==========统计细胞数量=========
yaqun_paint_tiaoxingtu = function(obj,category_yes,category_no) {
  
  #统计细胞个数
  dat = as.data.frame(table(Idents(obj)))
  #lable自动类型
  dat$label = paste(dat$Var1,dat$Freq,sep = ":")
  head(dat)
  library(ggplot2)
  library(paletteer)
  #View(palettes_d_names)
  
  p1 = ggplot(dat,aes(x = Freq,fill = Var1,y = Var1))+
    scale_fill_manual(values = palette)+#颜色
    geom_bar(stat = "identity")+
    theme_bw()+
    geom_text(aes(x = 0,label = label),hjust = 0)+
    theme(axis.text.y = element_blank(),   # 隐藏纵坐标刻度文字
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),# 隐藏纵坐标刻度线
          
          panel.grid.major = element_blank(),#隐藏背景网格
          panel.grid.minor = element_blank(),#隐藏背景网格
          legend.title = element_text(size = 0),
          legend.text = element_text(size = 18),
          #panel.border = element_blank(),
    )  

  
  #统计细胞个数
  dat = obj@meta.data[obj$category == category_yes,]
  dat =  as.data.frame(table(dat$celltype))
  
  
  #lable自动类型
  dat$label = paste(dat$Var1,dat$Freq,sep = ":")
  head(dat)
  library(ggplot2)
  library(paletteer)
  #View(palettes_d_names)
  
  p2 = ggplot(dat,aes(x = Freq,fill = Var1,y = Var1))+
    scale_fill_manual(values = palette)+
    geom_bar(stat = "identity")+
    theme_bw()+
    geom_text(aes(x = 0,label = label),hjust = 0)+
    theme(axis.text.y = element_blank(),   # 隐藏纵坐标刻度文字
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),# 隐藏纵坐标刻度线
          
          panel.grid.major = element_blank(),#隐藏背景网格
          panel.grid.minor = element_blank(),#隐藏背景网格
          legend.title = element_text(size = 0),
          legend.text = element_text(size = 18),
          #panel.border = element_blank(),
    )  
  
  #统计细胞个数
  
  dat = obj@meta.data[obj$category == category_no,]
  dat =  as.data.frame(table(dat$celltype))
  #lable自动类型
  dat$label = paste(dat$Var1,dat$Freq,sep = ":")
  head(dat)
  library(ggplot2)
  library(paletteer)
  #View(palettes_d_names)
  
  p3 = ggplot(dat,aes(x = Freq,fill = Var1,y = Var1))+
    scale_fill_manual(values = palette)+
    geom_bar(stat = "identity")+
    theme_bw()+
    geom_text(aes(x = 0,label = label),hjust = 0)+
    theme(axis.text.y = element_blank(),   # 隐藏纵坐标刻度文字
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),# 隐藏纵坐标刻度线
          
          panel.grid.major = element_blank(),#隐藏背景网格
          panel.grid.minor = element_blank(),#隐藏背景网格
          legend.title = element_text(size = 0),#图例标题字体大小
          legend.text = element_text(size = 18),#图例字体大小
          #panel.border = element_blank(),
    )  
  
  list("all" = p1,"yes" = p2, "no" = p3)
  
}


## 细胞比例
yaqun_paint_bili = function(obj){
  
  p1 = ggplot(obj@meta.data,aes(category,fill = celltype))+
    geom_bar(position = "fill", alpha = 0.9,width = 0.5)+
    scale_fill_manual(values = palette)+
    theme_classic()+
    coord_flip()+
    coord_fixed(ratio = 4)+ #纵轴长度是横轴的4倍
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(size = 16)
    )
  
  
  
  p2 = ggplot(obj@meta.data,aes(sample,fill = celltype))+
    geom_bar(position = "fill", alpha = 0.9,width = 0.5)+
    scale_fill_manual(values = palette)+
    theme_classic()+
    coord_flip()+
    coord_fixed(ratio = 4)+ #纵轴长度是横轴的4倍
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(size = 16)
    )
  
  list("two" = p1,"all" = p2)

}


##marker基因热图和气泡图
heatmap_paopaotu = function(obj,file){
  
  # 事例数据，以他为例top3pbmc.markers
  library(tidyverse)
  #构建这样的数据框
  genes = read.table(file, header=F, sep=",",comment.char = "", check.names =FALSE)
  genes = genes$V2
  genes = str_remove_all(genes," ")#删除空格
  p1 = jjDotPlot(object = obj,
                gene  = genes,
                anno = T,
                plot.margin = c(3,1,1,1),
                point.geom = T,
                tile.geom = F,
                id = "celltype",#纵坐标
                dot.col = c("#edf8b1","#2c7fb8")
  )
  
  p2 = jjDotPlot(object = obj,
                 gene  = genes,
                 anno = T,
                 plot.margin = c(3,1,1,1),
                 point.geom = F,
                 tile.geom = T,
                 id = "celltype",#纵坐标
                 dot.col = c("#edf8b1","#2c7fb8")
  )
  
  list("paopao" = p1,"heatmap" = p2)
  
}


## ==============火山图==========
getGoodHSplot = function(diff,genes,title){
  
  #定义颜色 
  pal <- c('#2128cf','#5d6c93','#c927a9')
  pal.grad <- colorRampPalette(pal[c(1,2,3)])(20)
  # 上下、左右边距
  width <- c(.5,1)
  # 长宽最大值
  ymax <- max(diff$neg_log10_padj)
  xmax <- max(abs(diff$avg_log2FC))
  # 背景渐变色插值，值越大渐变越平滑
  number <-1000
  # 渐变色数据
  rect.data <- expand.grid(
    x = seq(-xmax - width[2], xmax + width[2],
            length.out = number),
    y = seq(-width[1], ymax + width[1],
            length.out = number)
  )
  # 自定义 y 轴
  y_breaks <- scales::extended_breaks()(diff$neg_log10_padj)[-1]
  
  diff$symbol = rownames(diff)
  top.data = diff[diff$symbol %in% genes,]
  diff %>%
    ggplot(aes(avg_log2FC, neg_log10_padj, colour = tag)) +
    geom_raster(
      aes(x = x, y = y, fill = x),
      data = rect.data,
      alpha =   0.3,
      inherit.aes =   FALSE,
      show.legend =   FALSE
    ) +
    geom_point(
      aes(size = abs(avg_log2FC)), alpha =   0.6,   
      stroke =   0,
      show.legend =   FALSE
    ) +
    geom_text_repel(
      aes(avg_log2FC, neg_log10_padj, label = symbol),
      data = top.data,
      inherit.aes =   FALSE
    ) +
    geom_vline(xintercept =   0, linewidth =   0.5) +
    annotate(
      geom =   'segment',
      x = -0.2,
      y = y_breaks,
      xend =   0,
      yend = y_breaks,
      linewidth =   0.5
    ) +
    annotate(
      "text",
      x = -0.3,   
      y = y_breaks,
      label = y_breaks,
      hjust =   1, # 控制水平位置（负值向左，正值向右）
      color =   "black"
    ) +
    scale_color_manual(values = pal) +
    scale_fill_gradientn(colours = pal.grad) +
    scale_size_continuous(range = c(1,   4)) +
    scale_x_continuous(expand = expansion(c(0,   0))) +
    scale_y_continuous(expand = expansion(c(0,   0))) +
    labs(y =   '-log10(p.adjust)', title =  title) +
    theme_prism() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_line(size = 0.5, colour = "black"),#坐标轴
      axis.ticks = element_line(size = 0.5, colour = "black"),
      axis.text.x = element_text(size = 10 )
    )
  
}


##富集分析图================展示最富集的5条通路，可以自定义自己想展示的
goodFuji = function(diff,tag,colors){
  #devtools::install_github("dxsbiocc/gground")
  #install.packages("ggprism")
  library(gground)
  library(ggprism)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  pal = colors
  up_gene = rownames(diff)[which(diff$tag == tag)]
  up_ENTREZID = bitr(up_gene, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
  up_ENTREZID = up_ENTREZID$ENTREZID
  
  GO <- enrichGO(gene          = up_ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = "ENTREZID",
                 ont           = "ALL",       # 可选 "BP", "MF", "CC"
                 pAdjustMethod = "BH",        # 多重检验校正方法
                 pvalueCutoff  = 0.05,        # p值阈值
                 qvalueCutoff  = 0.1,         # q值阈值
                 readable      = TRUE) 
  
  KEGG <- enrichKEGG(gene         = up_ENTREZID,
                     organism     = "hsa",       # 人类：hsa，小鼠：mmu
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1,
                     use_internal_data = FALSE #FALSE为自动下载
  )
  KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  GO = as.data.frame(GO)
  KEGG = as.data.frame(KEGG)
  
  use_pathway <- group_by(GO, ONTOLOGY) %>%
    top_n(5, wt = -p.adjust) %>%
    group_by(p.adjust) %>%
    top_n(1, wt = Count) %>%
    rbind(
      top_n(KEGG,5, -p.adjust) %>%
        group_by(p.adjust) %>%
        top_n(1, wt = Count) %>%
        mutate(ONTOLOGY ='KEGG')
    ) %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY,
                             levels = rev(c('BP','CC','MF','KEGG')))) %>%
    dplyr::arrange(ONTOLOGY, p.adjust) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    tibble::rowid_to_column('index')
  
  
  # 左侧分类标签和基因数量点图的宽度
  width <-0.5
  # x 轴长度
  xaxis_max <- max(-log10(use_pathway$p.adjust)) +1
  # 左侧分类标签数据
  rect.data <- group_by(use_pathway, ONTOLOGY) %>%
    reframe(n = n()) %>%
    ungroup() %>%
    mutate(
      xmin = -3* width,
      xmax = -2* width,
      ymax = cumsum(n),
      ymin = lag(ymax, default =0) +0.6,
      ymax = ymax +0.4
    )
  
  
  
  p = use_pathway %>%
    ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Description), width =0.6, alpha =0.8
    ) +
    geom_text(
      aes(x =0.05, label = Description),
      hjust =0, size =5
    ) +
    
    #geom_text(
      #aes(x =0.1, label = geneID,y = index-0.1, colour = ONTOLOGY),
      #hjust =0, vjust =2.6, size =3.5, fontface ='italic',
      #show.legend =FALSE
    #) +
    # 基因数量
    geom_point(
      aes(x = -width, size = Count),
      shape =21
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name ='Count', range = c(5,16)) +
    # 分类标签
    geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      radius = unit(2,'mm'),
      inherit.aes =FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes =FALSE
    ) +
    geom_segment(
      aes(x =0, y =0, xend = xaxis_max, yend =0),
      linewidth =1.5,
      inherit.aes =FALSE
    ) +
    labs(y =NULL) +
    scale_fill_manual(name ='Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max,2),
      expand = expansion(c(0,0))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  p
  
}


##==========基因对比图==========
gene_to_gene = function(obj,celltypes,genes,category){
  
 plots = list()
 
 for (gene in genes){
   
   genelist = list()
   for(celltype in celltypes ) {
     celltype
     p = compare_OneGene(object = subset(obj,celltype == celltype),
                     gene = gene,
                     category = category,
                     cellname = celltype)
     
     genelist["ddaaa"] = p
   }
   
   plots[gene] = genelist
   
 }
 plots
}


#UMAP图展示
marker_umap = function(obj,genes,pt.size = 0.05){
  
  lists = list()
  for(gene in genes){
    p = plot_density(obj, features = c(gene),size = pt.size) +
      #scale_color_gradientn(colors = c("#A6ADCC", "#EDCAE0", "#F494BE"))+
      theme_void()+theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")  # 上、右、下、左
      )
    lists[[gene]] = p
    
  }
  
  lists
  
  
}
#########################################################################################################
#########################################################################################################
#小提琴图

marker_vlnplot = function(obj,genes,colors){
  
  
  library(reshape2)
  exp=GetAssayData(obj,slot = "data")
  vln.df=as.data.frame(exp[genes,])
  vln.df$gene=rownames(vln.df)
  vln.df=melt(vln.df,id="gene")
  colnames(vln.df)[c(2,3)]=c("barcode","exp")
  
  anno=sce.all@meta.data
  anno$barcode=rownames(anno)
  vln.df=inner_join(vln.df,anno,by="barcode")
  
  p = vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
    facet_grid(vln.df$gene~.,scales = "free_y")+
    scale_fill_manual(values = colors)+
    scale_x_discrete("")+scale_y_continuous("")+
    theme_classic()+
    theme(
      axis.text.x.bottom = element_text(angle = 90,hjust = 1,vjust = 1),
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  
  p
  
  
}

#########################################################################################################
#########################################################################################################
#气泡图
markerDot = function(obj,genes,colors,id = "celltype" ){
  
  
  p = jjDotPlot(object = obj,
                gene  = genes,
                anno = T,
                plot.margin = c(3,1,1,1),
                point.geom = T,
                tile.geom = F,
                id = id,#纵坐标
                dot.col = colors, 
                ytree = F, #去掉右边显示
                gene.order = genes
                
  )+theme(panel.grid = element_blank())
  p
  
}

#########################################################################################################
#########################################################################################################
#热图

marker_heatmap = function(obj,genes,colors,id = "celltype"){
  
  
  p = jjDotPlot(object = obj,
                gene  = genes,
                anno = T,
                plot.margin = c(3,1,1,1),
                point.geom = F,
                tile.geom = T,
                id = id,#纵坐标
                dot.col = colors, 
                ytree = F, #去掉右边显示
                gene.order = genes
                
  )+theme(panel.grid = element_blank())
  p
  
  
}

#########################################################################################################
#########################################################################################################


myggsave = function(file,plot,device = cairo_pdf,dpi = 600,width,height){
  
  
  dirname <- sub("/[^/]+$", "", file)
  print(dirname)
  if (!dir.exists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  }
  
  ggsave(file = file, plot = plot, device = device, dpi = dpi, width = width, height = height)
}