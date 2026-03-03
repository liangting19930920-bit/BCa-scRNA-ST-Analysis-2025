# ==============================================================================
# Project: Bladder Cancer Single-Cell & Spatial Analysis (ScBLCA)
# Description: Spatial Mapping, Colocalization, and Distance Analysis using CellTrek
# ==============================================================================

# ==============================================================================
# 0. 环境准备与依赖加载 (Dependencies)
# ==============================================================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(stringr)
library(viridis)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

# CellTrek 及其依赖 (如果没有安装，请使用 devtools::install_github("navinlabcode/CellTrek"))
library(CellTrek)
library(ConsensusClusterPlus)

# 全局路径配置
WORK_DIR <- "data/" # 建议用户将工作目录设置在项目根目录的 data 文件夹下
setwd(WORK_DIR)

# 全局调色板
my_cols <- c("#EDB931","#eb6841","#cc2a36","#00a0b0","#7A989A", "#849271", 
             "#CF9546", "#C67052", "#C1AE8D", "#3F6F76", "#C65840", "#62496F")

# ==============================================================================
# 1. 空间亚群特异性基因签名打分 (Spatial Signature Scoring)
# ==============================================================================
message(">>> Scoring Specific Subpopulation Signatures...")
# 假设 E_seob2 是已经注释好的上皮细胞单细胞对象，ST_H1 是高级别空间对象
Idents(E_seob2) <- "seurat_clusters"
all_markers <- FindAllMarkers(E_seob2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75) %>%
  filter(!str_starts(gene, "^(MT-|RPL|RPS)"))

# 提取特定 Cluster (如 Cluster 1) 的核心 Signature
cluster1_degs <- all_markers %>% filter(cluster == "1" & abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)
signature_cluster1 <- list(pull(cluster1_degs, gene))

# 在空间转录组中进行打分
ST_H1 <- AddModuleScore(object = ST_H1, features = signature_cluster1, name = 'Epi_Cluster1_Score')

[Image of spatial transcriptomics gene signature score map]

# 高级可视化：使用梯度颜色展示打分
SpatialFeaturePlot(ST_H1, features = "Epi_Cluster1_Score1", alpha = c(0.1, 1), 
                   min.cutoff = 0, max.cutoff = 1) &
  scale_fill_gradientn(colors = c("#000080", "#0000FF", "#00FF00", "#FFFF00", "#FFD700", "#FF0000", "#800000"),
                       limits = c(0, 1)) &
  theme_bw() &
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), legend.position = "top", 
        legend.title = element_blank(), legend.text = element_text(size = 10)) &
  guides(fill = guide_colorbar(barwidth = 7, barheight = 0.5))

# ==============================================================================
# 2. 单细胞与空间数据准备 (scRNA-seq & ST Prep for CellTrek)
# ==============================================================================
message(">>> Preparing Data for CellTrek Mapping...")
# 降采样单细胞数据以控制 CellTrek 内存消耗
set.seed(123)
Sc_H <- subset(seob_2, subset = status == "High_Grade")
Sc_H_sub <- subset(Sc_H, cells = sample(colnames(Sc_H), min(20000, ncol(Sc_H))))

# 确保条形码规范化 (非常重要，防止 Seurat 合并时引发命名冲突)
ST_H1 <- RenameCells(ST_H1, new.names = make.names(Cells(ST_H1)))
Sc_H_sub <- RenameCells(Sc_H_sub, new.names = make.names(Cells(Sc_H_sub)))

# ==============================================================================
# 3. CellTrek 空间映射 (CellTrek Co-embedding & Mapping)
# ==============================================================================
message(">>> Running CellTrek traint and celltrek...")
# 联合降维
H1_traint <- CellTrek::traint(st_data = ST_H1, sc_data = Sc_H_sub, 
                              sc_assay = 'RNA', cell_names = 'cell_type')

# 细胞映射
H1_celltrek <- CellTrek::celltrek(st_sc_int = H1_traint, int_assay = 'traint', 
                                  sc_data = Sc_H_sub, sc_assay = 'RNA', 
                                  reduction = 'pca', intp = TRUE, intp_pnt = 5000, 
                                  intp_lin = FALSE, nPCs = 30, ntree = 1000, 
                                  dist_thresh = 0.55, top_spot = 5, spot_n = 5, 
                                  repel_r = 20, repel_iter = 20, keep_model = TRUE)$celltrek

[Image of single cells mapped to spatial tissue spots using CellTrek]

# 可视化特定细胞亚群在空间的分布 (例如: Mono_TYMP, mCAF_POSTN)
Idents(H1_celltrek) <- "cell_type"
H1_subset_celltrek <- subset(H1_celltrek, idents = c("Mono_TYMP", "mCAF_POSTN"))

SpatialDimPlot(H1_subset_celltrek, pt.size = 1.5, stroke = 0.1) & 
  theme_bw() &
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank()) & ggtitle("High_Grade Mapping")

# saveRDS(H1_celltrek, file = "results/H1_CellTrek_Mapped.rds")

# ==============================================================================
# 4. 空间共定位分析 (Spatial Colocalization Network)
# ==============================================================================
message(">>> Running Spatial Colocalization (KL divergence)...")
# 计算共定位网络
H1_KL <- CellTrek::scoloc(H1_celltrek, col_cell = 'cell_type', use_method = 'KL', eps = 1e-50)

# 提取最小生成树 (MST) 矩阵
H1_mst_cons <- H1_KL$mst_cons
rownames(H1_mst_cons) <- colnames(H1_mst_cons)

# 构建 metadata (包含细胞类型及频率)
H1_cell_class <- H1_celltrek@meta.data %>% dplyr::select(id = cell_type) %>% unique()
H1_counts <- data.frame(freq = table(H1_celltrek$cell_type))
H1_class_meta <- merge(H1_cell_class, H1_counts, by.x = "id", by.y = "freq.Var1")

[Image of spatial cell colocalization minimum spanning tree network]

# 绘制共定位网络图
CellTrek::scoloc_vis(H1_mst_cons, meta_data = H1_class_meta)

# ==============================================================================
# 5. 跨细胞亚群物理距离分析 (Spatial Distance - kdist)
# ==============================================================================
message(">>> Calculating Spatial Distances (Endothelial vs. Myeloid/CAF)...")
# 提取坐标与细胞注释
inp_df <- H1_celltrek@meta.data %>% dplyr::select(cell_names = cell_type, coord_x, coord_y)

# 提示: 在 10X Visium 中，有时需要翻转 x 坐标以匹配图像方向 (按需开启)
# inp_df$coord_x <- 270 - inp_df$coord_x 

# 计算内皮细胞与特定髓系/CAF亚群的空间距离
query_cells <- c("Mono_TYMP", "Mono_SNCG", "Mac_SPP1", "Mac_APOE", "Mac_IGHG4", 
                 "Mac_NME2", "cDC_FCER1A", "cDC_LAMP3", "mCAF_POSTN")

dist_output <- kdist(inp_df = inp_df, ref = "Endothelial cells", ref_type = 'all', 
                     que = query_cells, k = 10, 
                     new_name = "Distance_to_Endothelial", keep_nn = FALSE)

# 合并距离结果与元数据
res_df <- dist_output$kdist_df
res_df$barcode <- rownames(res_df)
inp_df$barcode <- rownames(inp_df)
res_df <- left_join(res_df, inp_df, by = "barcode")

[Image of spatial distance boxplot between cell types]

# 绘制距离分布箱线图并进行克鲁斯卡尔-沃利斯检验 (Kruskal-Wallis)
p_dist <- ggboxplot(data = res_df, x = "cell_names", y = "Distance_to_Endothelial", 
                    fill = "cell_names", title = "K-distance to Endothelial cells") + 
  stat_compare_means(method = "kruskal.test", label.y.npc = "top") +
  theme(plot.title = element_text(color = "black", hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = "none") + 
  labs(y = "Spatial Distance", x = "Cell Type")

print(p_