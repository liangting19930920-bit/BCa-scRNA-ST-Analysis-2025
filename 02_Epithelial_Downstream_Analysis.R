# ==============================================================================
# Project: Bladder Cancer Single-Cell Analysis (ScBLCA)
# Description: Epithelial Cells Downstream Analysis (Sub-clustering, GSVA, InferCNV, Trajectory)
# ==============================================================================

# ==============================================================================
# 0. 环境准备与全局配置 (Dependencies & Config)
# ==============================================================================
# 数据与可视化
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

# 下游分析专用包
library(UCell)           # 干性评分
library(clusterProfiler) # GO/KEGG 富集
library(org.Hs.eg.db)    # 人类基因注释
library(GSVA)            # 通路评分
library(msigdbr)         # 基因集
library(infercnv)        # 拷贝数变异
library(monocle)         # 拟时序分析
library(CytoTRACE2)      # 细胞分化潜能

# 全局路径配置 (请其他用户运行前修改此处)
WORK_DIR <- "/path/to/your/working/directory"
INFERCNV_DIR <- paste0(WORK_DIR, "/InferCNV_output")
TRAJECTORY_DIR <- paste0(WORK_DIR, "/Trajectory_output")
dir.create(INFERCNV_DIR, showWarnings = FALSE)
dir.create(TRAJECTORY_DIR, showWarnings = FALSE)

setwd(WORK_DIR)

# ==============================================================================
# 1. 上皮细胞亚群提取与重聚类 (Sub-clustering)
# ==============================================================================
message(">>> Starting Epithelial Sub-clustering...")
# 提取上皮细胞并重新降维聚类
E_seob <- subset(seob, subset = cell_type == "Epithelial cells")

E_seob <- E_seob %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(object = .)) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.1, random.seed = 1) %>%
  RunUMAP(dims = 1:20)

# 可视化重聚类结果
my_cols <- c("#EDB931","#eb6841","#cc2a36","#00a0b0","#7A989A", "#849271", "#CF9546", "#C67052")
p_umap <- DimPlot(E_seob, reduction = "umap", group.by = "seurat_clusters", cols = my_cols, label = TRUE) +
  theme_classic() + theme(panel.grid = element_blank())
print(p_umap)

# Checkpoint 1: 保存上皮细胞对象
# saveRDS(E_seob, file = "E_seob_annotated.rds")

# ==============================================================================
# 2. 差异基因与 GO 富集分析 (DEGs & GO Enrichment)
# ==============================================================================
message(">>> Identifying DEGs and running GO Enrichment...")
Idents(E_seob) <- "seurat_clusters"
markers <- FindAllMarkers(E_seob, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75)

# 过滤线粒体与核糖体基因
sig_markers <- markers %>% 
  filter(p_val_adj < 0.05) %>%
  filter(!str_starts(gene, "^MT-|^RPL|^RPS"))

top5_markers <- sig_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# 选择 Cluster 10 作为示例进行 GO 分析
cluster10_degs <- sig_markers %>% filter(cluster == "10" & abs(avg_log2FC) > 0.5)
gene_list <- pull(cluster10_degs, gene)

de_ego <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                   ont = "ALL", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# ==============================================================================
# 3. 肿瘤干性评分 (Stemness Scoring by UCell)
# ==============================================================================
message(">>> Calculating Stemness Scores...")
stem_genes <- list(Stemness = c("CD44", "ALDH1A1", "SOX2", "NANOG", "OCT4", "BMI1", "ABCG2", "NOTCH1"))
E_seob <- AddModuleScore_UCell(E_seob, features = stem_genes, name = "_score")

stem_df <- FetchData(E_seob, vars = c("statues", "Stemness_score"))
stem_df$statues <- factor(stem_df$statues, levels = c("NAT", "Low_Grade", "High_Grade", "GC"))

p_stem <- ggplot(stem_df, aes(x = statues, y = Stemness_score, fill = statues)) +
  geom_violin(color = 'black', size = 1) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Stemness Score across Phenotypes", y = "Score", x = "")
print(p_stem)

# ==============================================================================
# 4. InferCNV 拷贝数变异分析 (InferCNV Analysis)
# ==============================================================================
message(">>> Running InferCNV...")
# 下采样以控制计算资源
set.seed(123)
epi_cells <- sample(rownames(E_seob@meta.data), min(5000, ncol(E_seob)))
ref_T_cells <- sample(rownames(seob@meta.data)[seob$cell_type == 'T/NK cells'], 500)
ref_Fib_cells <- sample(rownames(seob@meta.data)[seob$cell_type == 'Fibroblast cells'], 500)

subset_cells <- c(epi_cells, ref_T_cells, ref_Fib_cells)
sc_cnv <- subset(seob, cells = subset_cells)

# 准备表达矩阵和注释文件
counts_matrix <- GetAssayData(sc_cnv, layer = 'counts')
anno_df <- data.frame(cell = colnames(sc_cnv), cell_type = sc_cnv$cell_type, row.names = 1)

# 注意: hg38_gencode_v27.txt 需要提前下载并在目录中
# infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
#                                      annotations_file = anno_df,
#                                      gene_order_file = "hg38_gencode_v27.txt",
#                                      ref_group_names = c("T/NK cells", "Fibroblast cells"))
# infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=INFERCNV_DIR, 
#                               cluster_by_groups=TRUE, denoise=FALSE, HMM=FALSE)

# ==============================================================================
# 5. 单细胞轨迹推断 (Monocle 2 & CytoTRACE 2)
# ==============================================================================
message(">>> Running Monocle 2 & CytoTRACE 2 Trajectory Analysis...")
E_seob_sub <- subset(E_seob, cells = sample(colnames(E_seob), 6000))

# 构建 Monocle 对象
pd <- new("AnnotatedDataFrame", data = E_seob_sub@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(E_seob_sub), row.names = rownames(E_seob_sub)))
cds <- newCellDataSet(as.matrix(GetAssayData(E_seob_sub, layer = "counts")), 
                      phenoData = pd, featureData = fd, 
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds) %>% estimateDispersions()
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)$gene_id
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

# 拟时序密度图封装函数
create_density_plot <- function(data, x, colour, fill) {
  ggplot(data, aes_string(x = x, colour = colour, fill = fill)) +
    geom_density(bw = 0.5, size = 1, alpha = 0.5) + 
    theme_classic() + theme(legend.position = "right")
}
p_traj_density <- create_density_plot(pData(cds), "Pseudotime", "statues", "statues") + scale_fill_npg() + scale_color_npg()
print(p_traj_density)

# 结合 CytoTRACE 2 评估分化潜能
cytotrace2_sce <- cytotrace2(E_seob_sub, is_seurat = TRUE, slot_type = "counts", species = 'human', seed = 1234)
p_cyto <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  scale_colour_gradientn(colours = c("#000004FF", "#8C2981FF", "#FE9F6DFF", "#FCFDBFFF"),
                         name = "Relative\nOrder") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Checkpoint 2: 保存最终工作空间
# saveRDS(cds, file = paste0(TRAJECTORY_DIR, "/Monocle_cds.rds"))
message(">>> Analysis Complete.")
