# ==============================================================================
# Project: Bladder Cancer Single-Cell Analysis (ScBLCA)
# Description: Fibroblast Sub-clustering, Annotation, and Functional Enrichment
# ==============================================================================

# ==============================================================================
# 0. 环境准备与参数配置 (Environment Setup)
# ==============================================================================
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(msigdbr)
library(irGSEA)
library(UCell)
library(cols4all)
library(tidydr)
library(plyr)

# 全局路径配置
WORK_DIR <- "/path/to/your/working/directory"
setwd(WORK_DIR)

# 设置疾病分级标准顺序
status_levels <- c("NAT", "GC", "Low_Grade", "High_Grade")

# ==============================================================================
# 1. 成纤维细胞子集提取与重聚类 (Sub-clustering)
# ==============================================================================
# 假设基础对象 seob 已加载
F_seob <- subset(seob, subset = cell_type == "Fibroblast cells")
F_seob$status <- factor(F_seob$status, levels = status_levels)

# 数据标准化与降维
F_seob <- F_seob %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(object = .), reduction.name = "pca") %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.05, random.seed = 1) %>%
  RunTSNE(dims = 1:20) %>%
  RunUMAP(dims = 1:20)

# 剔除潜在的噪声或混合群 (例如 Cluster 9, 10)
F_seob2 <- subset(F_seob, idents = c("9", "10"), invert = TRUE)

# ==============================================================================
# 2. 亚群注释与特征可视化 (Annotation & Visualization)
# ==============================================================================
# 基于特征基因进行亚群精确注释
F_seob2$cell_type <- plyr::mapvalues(
  F_seob2$seurat_clusters,
  from = 0:8,
  to = c("mCAF_NME2", "iCAF_CFD", "iCAF_CCL11", "mCAF_ECRG4", 
         "mCAF_CCL5", "mCAF_POSTN", "iCAF_HLA-DRB1", "iCAF_STC1", "mCAF_SEPT11")
)
Idents(F_seob2) <- "cell_type"

# --- 2.1 高级 UMAP 可视化 ---

UMAP_df <- cbind(
  as.data.frame(Embeddings(F_seob2, "umap")),
  F_seob2@meta.data[, c("seurat_clusters", "cell_type", "status")]
)

mytheme <- theme_void() + theme(plot.margin = margin(5.5, 15, 5.5, 5.5)) 
mycol <- c4a('wright25', length(unique(F_seob2$cell_type)))

p_umap <- ggplot(data = UMAP_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cell_type), size = 0.5, alpha = 0.8) + 
  mytheme + 
  theme_dr(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) + 
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  scale_color_manual(values = mycol)
# ggsave("UMAP_Fibroblast_Subtypes.pdf", plot = p_umap, width = 7, height = 6)

# --- 2.2 标记基因点图 (DotPlot) ---
marker_genes <- c("NME2", "CFD", "CCL11", "ECRG4", "CCL5", "POSTN", "HLA-DRB1", "STC1", "SEPT11")
p_dot <- DotPlot(F_seob2, features = marker_genes, group.by = "cell_type") +
  coord_flip() + theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 10)) +
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

# ==============================================================================
# 3. 细胞组成比例分析 (Cell Fraction Analysis)
# ==============================================================================
celltype_ratio <- F_seob2@meta.data %>%
  group_by(status, cell_type) %>%
  dplyr::summarise(n = n(), .groups = 'drop') %>%
  mutate(relative_freq = n / sum(n))

celltype_ratio$cell_type <- factor(celltype_ratio$cell_type)
celltype_ratio$status <- factor(celltype_ratio$status, levels = status_levels)

df_percent <- ddply(celltype_ratio, 'cell_type', transform, percent = n / sum(n))

# 绝对数量分布图
p_count <- ggplot(celltype_ratio, aes(x = cell_type, y = n, fill = status)) + 
  geom_bar(stat = 'identity', position = "stack") + 
  scale_fill_manual(values = c("#6EB467","#80689C","#CEB481","#D86F6F")) + 
  labs(x = "", y = "Cell Number") + coord_flip() + theme_classic() +
  theme(axis.title.x = element_text(size = 12, face = 'bold'), 
        axis.text.y = element_text(size = 10, face = 'bold'), 
        legend.position = "bottom")

# 相对比例分布图
p_ratio <- ggplot(df_percent, aes(x = cell_type, y = percent, fill = status)) + 
  geom_bar(stat = 'identity', position = "stack") + 
  scale_fill_manual(values = c("#6EB467","#80689C","#CEB481","#D86F6F")) + 
  labs(x = "", y = "Sample Ratio") + coord_flip() + theme_classic() +
  theme(axis.title.x = element_text(size = 12, face = 'bold'),
        axis.text.y = element_blank(), 
        legend.position = "none")

# p_combined_ratio <- grid.arrange(p_count, p_ratio, ncol = 7, layout_matrix = rbind(c(1,1,1,1,1,2,2)))

# ==============================================================================
# 4. 功能富集分析 (Functional Enrichment)
# ==============================================================================

# --- 4.1 GSVA 通路评分 ---
genesets <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  subset(select = c("gs_name", "gene_symbol")) %>% 
  as.data.frame() %>% 
  split(.$gene_symbol, .$gs_name)

expr_mat <- AverageExpression(F_seob2, assays = "RNA", slot = "data")[[1]]
expr_mat <- as.matrix(expr_mat[rowSums(expr_mat) > 0, ])

gsva_res <- gsva(expr_mat, genesets, method = "gsva")
# pheatmap::pheatmap(gsva_res, show_colnames = TRUE, scale = "row", angle_col = "45", filename = "GSVA_heatmap.pdf")

# --- 4.2 特异性亚群 GO 分析 (例如: mCAF_POSTN) ---
F_markers <- FindAllMarkers(F_seob2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75) %>%
  filter(!str_starts(gene, "MT-"))

postn_degs <- F_markers %>% filter(cluster == "mCAF_POSTN" & abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)

ego_postn <- enrichGO(gene = pull(postn_degs, gene), OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL', ont = "ALL", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# --- 4.3 irGSEA 单细胞水平通路富集 ---
F_score <- irGSEA.score(object = F_seob2, assay = "RNA", slot = "data", seeds = 123, ncores = 1,
                        msigdb = TRUE, species = "Homo sapiens", category = "H",
                        method = c("AUCell", "UCell", "singscore"), kcdf = 'Gaussian')

result_dge <- irGSEA.integrate(object = F_score, group.by = "cell_type", 
                               method = c("AUCell", "UCell", "singscore"))

# irGSEA.heatmap(object = result_dge, method = "RRA", top = 50)
# irGSEA.ridgeplot(object = F_score, method = "UCell", show.geneset = "HALLMARK-ANGIOGENESIS", group.by = "cell_type")