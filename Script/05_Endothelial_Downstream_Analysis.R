# ==============================================================================
# Project: Bladder Cancer Single-Cell Analysis (ScBLCA)
# Description: Endothelial Cells Downstream Analysis (Sub-clustering, Metabolism, PROGENy)

# ==============================================================================

# ==============================================================================
# 0. 环境准备与全局配置
# ==============================================================================
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(grid)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(stringr)

# 专有分析工具包
library(clusterProfiler)
library(org.Hs.eg.db)
library(scRNAtoolVis)
library(GSVA)
library(msigdbr)
library(irGSEA)
library(UCell)
library(scMetabolism)
library(progeny)
library(tidydr)
library(cols4all)

# 全局路径与参数设定
WORK_DIR <- "/path/to/your/working/directory"
setwd(WORK_DIR)
status_levels <- c("NAT", "GC", "Low_Grade", "High_Grade")

# ==============================================================================
# 1. 内皮细胞亚群提取与重聚类
# ==============================================================================
# 假设基础对象 seob_2 已加载
En_seob <- subset(seob_2, subset = cell_type == "Endothelial cells")

# 基础降维与聚类
En_seob <- En_seob %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(reduction.name = "pca") %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.1, random.seed = 1) %>%
  RunTSNE(dims = 1:20) %>%
  RunUMAP(dims = 1:20)

# 剔除异常群 (Cluster 9)
En_seob2 <- subset(En_seob, idents = c("9"), invert = TRUE)

# ==============================================================================
# 2. 亚群精确注释与 Marker 可视化
# ==============================================================================
# 亚群重命名
En_seob2$cell_type <- plyr::mapvalues(
  En_seob2$seurat_clusters,
  from = 0:8,
  to = paste0("En", 0:8)
)
En_seob2$status <- factor(En_seob2$status, levels = status_levels)
Idents(En_seob2) <- "cell_type"



# 经典 Marker DotPlot (附带定制化 Grid 分界线)
genes_to_check <- list(
  Arteries = c('GJA5', 'FBLN5', "GJA4"),
  Veins = c("ACKR1", 'SELP', "CLU"),
  Capillaries = c("CA4", "CD36", "RGCC"),
  Tip_vEC = c("COL4A1", 'KDR', "ESM1")
)

p_dot <- DotPlot(En_seob2, features = genes_to_check, cols = c("grey", "red"), cluster.idents = TRUE) +
  RotatedAxis() +
  theme(panel.border = element_rect(color="black"),
        panel.spacing = unit(1, "mm"),
        strip.text = element_text(margin=margin(b=3, unit="mm")),
        strip.placement = 'outside',
        axis.line = element_blank()) + labs(x="", y="")

p_grob <- ggplotGrob(p_dot)
lg <- linesGrob(x=unit(c(0,1),"npc"), y=unit(c(0,0)+0.2,"npc"), gp=gpar(col="black", lwd=4))
grid.newpage()
for (k in grep("strip-t", p_grob$layout$name)) {
  p_grob$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}
grid.draw(p_grob)

# ==============================================================================
# 3. 差异基因与 GO 富集分析
# ==============================================================================
En_markers <- FindAllMarkers(En_seob2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75) %>%
  filter(!str_starts(gene, "MT-"))

# 特定亚群 (En2, En3, En4, En7) GO 富集
target_clusters <- c("En2", "En3", "En4", "En7")
H_de_target <- En_markers %>% 
  filter(cluster %in% target_clusters & abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)

H_de_ego_target <- enrichGO(gene = H_de_target$gene, OrgDb = org.Hs.eg.db, 
                            keyType = 'SYMBOL', ont = "ALL", 
                            qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# 转换并筛选绘图数据
dt_ego <- as.data.frame(H_de_ego_target)
dt_ego$logp <- -log10(dt_ego$pvalue)
dt_ego$Description <- factor(dt_ego$Description, levels = rev(dt_ego$Description))



ggdotchart(dt_ego[1:20, ], x = "Description", y = "Count", color = "logp",
           sorting = "descending", add = "segments", rotate = TRUE, dot.size = 5,
           label = round(dt_ego$Count[1:20]), font.label = list(color = "white", size = 9, vjust = 0.5),
           ggtheme = theme_pubr()) +
  scale_color_viridis_c(option = "D")


# ==============================================================================
# 4. PROGENy 通路活性打分
# ==============================================================================
En_seob2 <- progeny(En_seob2, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
En_seob2 <- Seurat::ScaleData(En_seob2, assay = "progeny")

progeny_scores_df <- as.data.frame(t(GetAssayData(En_seob2, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

CellsClusters <- data.frame(Cell = names(Idents(En_seob2)), 
                            CellType = as.character(Idents(En_seob2)), stringsAsFactors = FALSE)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters, by = "Cell")

summarized_progeny_scores_df <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), .groups = 'drop') %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE)



pheatmap(t(summarized_progeny_scores_df), fontsize = 12, fontsize_row = 10, 
         color = turbo(90), main = "PROGENy Pathway Activity", 
         angle_col = 90, treeheight_col = 0, border_color = NA)