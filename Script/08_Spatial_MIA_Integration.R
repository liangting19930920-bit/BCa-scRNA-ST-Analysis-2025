# ==============================================================================
# Project: Bladder Cancer Single-Cell & Spatial Analysis (ScBLCA)
# Description: Multimodal Intersection Analysis (MIA) & Functional Region Mapping
# ==============================================================================

# ==============================================================================
# 0. 环境准备与全局配置 (Dependencies & Config)
# ==============================================================================
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(patchwork)

# 全局路径配置 (请根据实际情况修改)
WORK_DIR <- "data/"
MIA_SCRIPT_PATH <- "scripts/syMIA.R" # 将自定义 MIA 脚本放在 scripts 目录下
setwd(WORK_DIR)

# ==============================================================================
# 1. 单细胞数据准备与特征基因提取 (scRNA-seq Markers)
# ==============================================================================
message(">>> Extracting scRNA-seq Markers for MIA...")

# 提取高级别单细胞数据并合并精细注释的亚群 (例如髓系、成纤维)
# 假设 seob_2 为基础注释，Mye_seob2 和 F_seob2 为精细注释
Idents(seob_2) <- "cell_type"
seob_other <- subset(seob_2, idents = c("Myeloid cells", "Fibroblast cells"), invert = TRUE)

# 合并精细化注释的对象
Sc_merged <- merge(x = seob_other, y = list(Mye_seob2, F_seob2))

# 提取高级别样本 (High Grade) 
Sc_H <- subset(Sc_merged, subset = status == "High_Grade")

# (可选) 下采样以加速计算
set.seed(123)
Sc_H_sub <- subset(Sc_H, cells = sample(colnames(Sc_H), min(10000, ncol(Sc_H))))
Sc_H_sub <- SCTransform(Sc_H_sub, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

# 寻找单细胞各亚群的特征基因
Idents(Sc_H_sub) <- "cell_type"
sc_markers <- FindAllMarkers(Sc_H_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)

# 严格过滤单细胞 Marker (p_val_adj < 0.05 且表达比例差异显著)
sc_markers_sig <- sc_markers %>% 
  filter(p_val_adj < 0.05) %>%
  mutate(diff_pct = pct.1 - pct.2) %>%
  filter(diff_pct > 0.2) %>%
  arrange(cluster, desc(avg_log2FC))

# ==============================================================================
# 2. 空间转录组特征基因提取 (Spatial Transcriptomics Markers)
# ==============================================================================
message(">>> Extracting Spatial Region Markers for MIA...")

# 假设 ST_H1 已经完成了区域划分 (region 0-6)
Idents(ST_H1) <- "region"

sp_markers <- FindAllMarkers(ST_H1, only.pos = TRUE, min.pct = 0.25, test.use = 'wilcox', logfc.threshold = 0.2)

# 严格过滤空间 Marker (去除线粒体，要求显著)
sp_markers_sig <- sp_markers %>% 
  filter(!str_starts(gene, "^MT-")) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(diff_pct = pct.1 - pct.2) %>%
  filter(diff_pct > 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

# ==============================================================================
# 3. 多模态交叉分析 (Multimodal Intersection Analysis - MIA)
# ==============================================================================
message(">>> Running MIA Integration...")

# 准备 MIA 格式输入
celltype_specific <- sc_markers_sig[, c("cluster", "gene")]
colnames(celltype_specific) <- c("celltype", "gene")

region_specific <- sp_markers_sig[, c("cluster", "gene")]
colnames(region_specific) <- c("region", "gene")

# 计算背景基因总数 (Universe N)
N_universe <- length(union(rownames(Sc_H_sub), rownames(ST_H1)))

# 设置区域颜色
color_region <- c("#666666", "#766fb1", "#e42a88", "#189d77", "#ff9933", "#33ccff", "#ffcc00")
names(color_region) <- paste0("region ", 0:6)



# 运行 MIA (调用外部脚本)
source(MIA_SCRIPT_PATH)
mia_results <- syMIA(region_specific, celltype_specific, N_universe, color_region)

# ==============================================================================
# 4. 空间功能区域定义与可视化 (Functional Region Mapping)
# ==============================================================================
message(">>> Mapping Tumor vs Stroma Regions...")

# 基于 MIA 的映射结果，将特定 cluster 定义为 Stroma (间质) 或 Cancer (肿瘤实质)
# 注意：原代码中 region 是基于 seurat_clusters 映射的，这里直接对 cluster 进行分类
ST_H1@meta.data$Functional_Region <- NA
ST_H1@meta.data$Functional_Region[ST_H1@meta.data$seurat_clusters %in% c('0', '1', '5')] <- "Stroma"
ST_H1@meta.data$Functional_Region[ST_H1@meta.data$seurat_clusters %in% c('2', '3', '4', '6')] <- "Cancer"



# 空间可视化 (使用更新后的 SpatialDimPlot 替代废弃的 SpatialPlot)
p_functional_map <- SpatialDimPlot(ST_H1, group.by = 'Functional_Region', 
                                   cols = c("Cancer" = '#e6353f', "Stroma" = '#00a9e0'), 
                                   stroke = 0, pt.size.factor = 1.2) +
  theme_bw() + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = "right", legend.title = element_text(face = "bold")) +
  ggtitle("High Grade: Cancer vs Stroma Architecture")

print(p_functional_map)

# ggsave("results/Functional_Region_Map.pdf", plot = p_functional_map, width = 6, height = 5)
message(">>> MIA Analysis and Mapping Complete.")