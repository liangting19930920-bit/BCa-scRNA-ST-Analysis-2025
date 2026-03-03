# ==============================================================================
# Project: Bladder Cancer Single-Cell & Spatial Analysis (ScBLCA)
# Description: Spatial Transcriptomics Pipeline (Harmony Integration, GSVA, CellChat)
# ==============================================================================

# ==============================================================================
# 0. 环境准备与全局参数配置 (Dependencies & Config)
# ==============================================================================
# 基础分析与可视化包
library(Seurat)
library(harmony)
library(hdf5r)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(tidydr)
library(pheatmap)
library(cowplot)

# 富集分析与高级可视化包
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(msigdbr)
library(ggforce)     # 用于流星图
library(ggraph)
library(tidygraph)

# 空间细胞通讯包
library(CellChat)

# 开启多线程并突破内存限制
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 2000 * 1024^2) 

# 设置通用相对路径 (指向标准的 10X SpaceRanger 输出目录)
DIR_H1 <- "data/High_Grade/outs/"
DIR_L1 <- "data/Low_Grade/outs/"
OUT_DIR <- "results/"

# 全局调色板
my_cols <- c("#EDB931","#eb6841","#cc2a36","#00a0b0","#7A989A", "#849271", 
             "#CF9546", "#C67052", "#C1AE8D", "#3F6F76")

# ==============================================================================
# 1. 数据读取与基础合并 (Data Loading & Merge)
# ==============================================================================
message(">>> Loading H5 Spatial Data...")
ST_H1 <- Load10X_Spatial(data.dir = DIR_H1, filename = "filtered_feature_bc_matrix.h5", 
                         assay = "Spatial", slice = "High_Grade")
ST_H1$sample <- "High_Grade"
ST_H1$status <- "High_Grade"

ST_L1 <- Load10X_Spatial(data.dir = DIR_L1, filename = "filtered_feature_bc_matrix.h5", 
                         assay = "Spatial", slice = "Low_Grade")
ST_L1$sample <- "Low_Grade"
ST_L1$status <- "Low_Grade"

# 基础合并 (不直接改变表达量，保护图像信息)
spatial_list <- list(High_Grade = ST_H1, Low_Grade = ST_L1)
Spatial_integrated <- merge(x = spatial_list[[1]], y = spatial_list[-1], add.cell.ids = names(spatial_list))

# 修复 Seurat 合并后易丢失多图像的 Bug
Spatial_integrated@images <- list(
  High_Grade = spatial_list[["High_Grade"]]@images$High_Grade,
  Low_Grade  = spatial_list[["Low_Grade"]]@images$Low_Grade
)
Spatial_integrated$CellID <- rownames(Spatial_integrated@meta.data)

# ==============================================================================
# 2. SCTransform 与 Harmony 批次校正 (Harmony Integration)
# ==============================================================================
message(">>> Running SCTransform and Harmony...")
Spatial_integrated <- SCTransform(Spatial_integrated, assay = "Spatial", verbose = FALSE)
Spatial_integrated <- RunPCA(Spatial_integrated, assay = "SCT", verbose = FALSE)

# 核心：使用 Harmony 消除多样本技术批次效应
Spatial_integrated <- RunHarmony(Spatial_integrated, group.by.vars = "sample", 
                                 assay.use = "SCT", project.dim = FALSE)

Spatial_integrated <- RunUMAP(Spatial_integrated, reduction = "harmony", dims = 1:20)
Spatial_integrated <- FindNeighbors(Spatial_integrated, reduction = "harmony", dims = 1:20)
Spatial_integrated <- FindClusters(Spatial_integrated, resolution = 0.3)

# 区域 (Region) 规范化命名
region_map <- setNames(paste0("region ", 0:max(as.numeric(as.character(Spatial_integrated$seurat_clusters)))), 
                       0:max(as.numeric(as.character(Spatial_integrated$seurat_clusters))))
Spatial_integrated$region <- unname(region_map[as.character(Spatial_integrated$seurat_clusters)])
Idents(Spatial_integrated) <- "region"



# 降维与原位可视化
p_umap <- DimPlot(Spatial_integrated, reduction = "umap", group.by = "region", cols = my_cols, pt.size = 0.1) + 
  theme_dr() + theme(panel.grid = element_blank())

p_spatial <- SpatialDimPlot(Spatial_integrated, group.by = "region", stroke = 0.1, ncol = 2) & 
  scale_fill_manual(values = my_cols) & theme_bw() & 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

# ==============================================================================
# 3. 肿瘤微环境功能基因集打分 (Functional Gene Set Scoring)
# ==============================================================================
message(">>> Scoring TME Signatures...")
signatures <- list(
  Hypoxia = c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "CA9", "EPO", "BNIP3", "ENO1", "HK2"),
  Stemness = c("CD44", "ALDH1A1", "SOX2", "NANOG", "OCT4", "BMI1", "ABCG2", "NOTCH1", "SDC1"),
  Angiogenesis = c("VEGFA", "ANGPT2", "PDGFB", "FGF2", "HIF1A", "TGFB1", "FLT1", "KDR", "PECAM1")
)

for (sig_name in names(signatures)) {
  Spatial_integrated <- AddModuleScore(Spatial_integrated, features = list(signatures[[sig_name]]), name = sig_name)
}

# 空间评分映射图
SpatialFeaturePlot(Spatial_integrated, features = "Hypoxia1", alpha = c(0.1,1), min.cutoff = 0) & theme_bw()

# 评分小提琴图
ggviolin(Spatial_integrated@meta.data, x = "region", y = "Stemness1", width = 0.8, color = "black", 
         fill = "region", add = 'mean_sd', bxp.errorbar = TRUE, palette = "npg", legend = "right")

# ==============================================================================
# 4. 空间差异表达与定制化富集流星图 (Spatial DEGs & Meteor Plot)
# ==============================================================================
message(">>> Calculating DEGs and running Enrichment Analysis...")
Spatial_integrated <- PrepSCTFindMarkers(Spatial_integrated)
Markers <- FindAllMarkers(Spatial_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  filter(!str_starts(gene, "^(MT-|RPL|RPS)"))

# 获取 Region 3 的差异基因进行富集
H_de3 <- Markers %>% filter(cluster == "region 3" & p_val_adj < 0.05)
H_de_ego3 <- enrichGO(gene = H_de3$gene, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', 
                      ont = "ALL", qvalueCutoff = 0.05, pvalueCutoff = 0.05)



# --- 4.1 绘制高水平 GO 富集流星图 ---
plot_go_meteor <- function(go_result) {
  df <- as.data.frame(go_result) %>% group_by(ONTOLOGY) %>% slice_head(n = 10) %>% arrange(desc(pvalue))
  df$ratio <- sapply(df$GeneRatio, function(x) eval(parse(text = x)))
  df$Description <- factor(df$Description, levels = df$Description)
  
  ggplot(df) +
    geom_link(aes(x = 0, y = Description, xend = -log10(pvalue), yend = Description,
                  alpha = stat(index), color = ONTOLOGY, size = after_stat(index)), n = 500, show.legend = FALSE) +
    geom_point(aes(x = -log10(pvalue), y = Description), color = "black", fill = "white", size = 6, shape = 21) +
    geom_line(aes(x = ratio*100, y = Description, group = 1), orientation = "y", linewidth = 1, color = "#FFCC00") +
    scale_x_continuous(sec.axis = sec_axis(~./100, labels = scales::label_percent(), name = "Percent of geneRatio")) +
    theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(color = "black")) +
    labs(y = "", x = "-log10 Pvalue") +
    facet_wrap(~ONTOLOGY, scales = "free", ncol = 1) + scale_color_brewer(palette = "Set1")
}
p_meteor <- plot_go_meteor(H_de_ego3)

# ==============================================================================
# 5. 特定区域高级别 vs 低级别差异分析 (Region-Specific Group Comparison)
# ==============================================================================
# 自定义富集柱状图函数
GO_a <- function(interest_gene, global_gene) {
  gene_entrez <- bitr(interest_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID
  universe_entrez <- bitr(global_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID
  
  enr_res <- enrichGO(gene = gene_entrez, universe = universe_entrez, ont = "BP", OrgDb = 'org.Hs.eg.db')
  enr_res2 <- simplify(enr_res)
  list(go_res = enr_res, go_res_simplify = enr_res2)
}

Spatial_sub <- subset(Spatial_integrated, region %in% c("region 3", "region 5"))
Idents(Spatial_sub) <- "status"
Spatial_sub <- PrepSCTFindMarkers(Spatial_sub)

de_h_v_l <- FindMarkers(Spatial_sub, ident.1 = "High_Grade", ident.2 = "Low_Grade", 
                        test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25) %>%
  rownames_to_column("gene") %>% filter(!str_starts(gene, "^(MT-)"))

# ==============================================================================
# 6. Spatial CellChat 空间组细胞通讯推断 (Spatial Cell Communication)
# ==============================================================================
message(">>> Building Spatial CellChat Network for High Grade...")
ST_H1_chat <- subset(Spatial_integrated, sample == "High_Grade")
ST_H1_chat@images[-which(names(ST_H1_chat@images) %in% "High_Grade")] <- NULL
Idents(ST_H1_chat) <- "region"

data.input <- GetAssayData(ST_H1_chat, slot = "data", assay = "SCT")
meta <- data.frame(labels = Idents(ST_H1_chat), row.names = names(Idents(ST_H1_chat)))
spatial.locs <- GetTissueCoordinates(ST_H1_chat, scale = NULL, cols = c("imagerow", "imagecol"))

# 根据 10X Visium 官方定义转换物理距离
spot.size <- 65 
spatial.factors <- data.frame(ratio = spot.size/100, tol = spot.size/2) # 简化比例

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, 
                           spatial.factors = spatial.factors)

cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat) %>%
  identifyOverExpressedGenes() %>%
  identifyOverExpressedInteractions(variable.both = FALSE) %>%
  computeCommunProb(type = "truncatedMean", trim = 0.1, distance.use = TRUE, 
                    interaction.range = 250, scale.distance = 0.01, 
                    contact.dependent = TRUE, contact.range = 100) %>%
  filterCommunication(min.cells = 10) %>%
  computeCommunProbPathway() %>%
  aggregateNet()



# 空间通讯可视化展示 (例如 VEGF 通路)
pathways.show <- c("VEGF")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", 
                    edge.width.max = 2, alpha.image = 0.2, vertex.size.max = 4, vertex.label.cex = 3.5)

# saveRDS(cellchat, file = "ST_Cellchat_HighGrade.rds")
message(">>> Pipeline Successfully Completed.")