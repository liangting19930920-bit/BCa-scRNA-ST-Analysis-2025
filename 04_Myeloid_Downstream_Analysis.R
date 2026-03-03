# ==============================================================================
# Project: Bladder Cancer Single-Cell Analysis (ScBLCA)
# Description: Myeloid Cells Downstream Analysis (Clustering, Trajectory, SCENIC)

# ==============================================================================

# ==============================================================================
# 0. 环境准备与依赖加载 (Dependencies)
# ==============================================================================
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(grid)
library(gridExtra)
library(cols4all)
library(tidydr)
library(plyr)
library(ggsci)
library(viridis)
library(ggrastr)
library(ggforce)
library(clusterProfiler)
library(org.Hs.eg.db)
library(monocle)
library(irGSEA)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(ggheatmap) # devtools::install_github("XiaoLuo-boy/ggheatmap")
library(ggraph)
library(tidygraph)
library(stringr)
library(cowplot)

# 全局路径配置
WORK_DIR <- "/path/to/your/working/directory"
setwd(WORK_DIR)
# 注：在此将原代码中的 statues 统一更正为了规范拼写 status
status_levels <- c("NAT", "GC", "Low_Grade", "High_Grade")

# ==============================================================================
# 1. 髓系细胞提取与重聚类 (Sub-clustering)
# ==============================================================================
Mye_seob <- subset(seob_2, subset = cell_type == "Myeloid cells")
Mye_seob$status <- factor(Mye_seob$status, levels = status_levels)

Mye_seob <- Mye_seob %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(reduction.name = "pca") %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.1, random.seed = 1) %>%
  RunTSNE(dims = 1:20) %>%
  RunUMAP(dims = 1:20)

# 剔除异常群 (Cluster 12)
Mye_seob2 <- subset(Mye_seob, idents = c("12"), invert = TRUE)

# ==============================================================================
# 2. 亚群精确注释与 Marker 可视化 (Annotation & Marker Viz)
# ==============================================================================
Mye_seob2$cell_type <- plyr::mapvalues(
  Mye_seob2$seurat_clusters,
  from = 0:11,
  to = c("Mac_NME2", "Mac_SPP1", "Mac_IGHG4", "Neu_IL1B", "cDC_FCER1A", "Mac_APOE",
         "Mono_TYMP", "cDC_LAMP3", "Mono_SNCG", "Neu_CD79A", "Mast_KIT", "Mast_KIT")
)
Idents(Mye_seob2) <- "cell_type"

# --- 2.1 经典 Marker DotPlot ---
list_genes <- list(
  Mac = c('C1QA', 'C1QB', "CD68"),
  Mono = c("S100A9", 'S100A8', "CD14", "CST3"),
  Neu = c("CSF3R", "IFITM2", "G0S2"),
  Mast = c("TPSAB1", 'TPSB2', "KIT"),
  cDC_FCER1A = c("FCER1A", 'CD1C'),
  cDC_LAMP3 = c('LAMP3', 'IDO1')
)

p_dot <- DotPlot(Mye_seob2, features = list_genes, cols = c("grey", "red"), cluster.idents = TRUE) +
  RotatedAxis() +
  theme(panel.border = element_rect(color="black"),
        panel.spacing = unit(1, "mm"),
        strip.text = element_text(margin=margin(b=3, unit="mm")),
        strip.placement = 'outside',
        axis.line = element_blank()) + 
  labs(x="", y="")

# 使用 grid 增加分界线
p_grob <- ggplotGrob(p_dot)
lg <- linesGrob(x=unit(c(0,1),"npc"), y=unit(c(0,0)+0.2,"npc"), gp=gpar(col="black", lwd=4))
grid.newpage()
for (k in grep("strip-t", p_grob$layout$name)) {
  p_grob$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}
grid.draw(p_grob)

# ==============================================================================
# 3. 细胞比例计算与可视化 (Cell Fractions)
# ==============================================================================
celltype_ratio <- Mye_seob2@meta.data %>%
  group_by(status, cell_type) %>%
  dplyr::summarise(n = n(), .groups = 'drop') %>%
  mutate(relative_freq = n / sum(n))

celltype_ratio$status <- factor(celltype_ratio$status, levels = status_levels)
df_percent <- ddply(celltype_ratio, 'cell_type', transform, percent = n / sum(n))

# 细胞数目累积柱状图
p_count <- ggplot(celltype_ratio, aes(x = cell_type, y = n, fill = status)) + 
  geom_bar(stat = 'identity', position = "stack") + 
  scale_fill_manual(values = c("#6EB467","#80689C","#CEB481","#D86F6F")) + 
  labs(x = "", y = "Cell Number") + coord_flip() + theme_classic() +
  theme(axis.title.x = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=10, face='bold'),
        legend.position = c(0.70, 0.72))

# 细胞百分比柱状图
p_percent <- ggplot(df_percent, aes(x = cell_type, y = percent, fill = status)) + 
  geom_bar(stat = 'identity', position = "stack") + 
  scale_fill_manual(values = c("#6EB467","#80689C","#CEB481","#D86F6F")) + 
  labs(x = "", y = "Sample Ratio") + coord_flip() + theme_classic() +
  theme(axis.text.y = element_blank(), legend.position = "none")

grid.arrange(p_count, p_percent, ncol = 7, layout_matrix = rbind(c(1,1,1,1,1,2,2)))

# ==============================================================================
# 4. Monocle 2 拟时序分析 (Trajectory Analysis)
# ==============================================================================
Mye_traj_seob <- subset(Mye_seob2, cell_type %in% c("Mac_NME2", "Mac_SPP1", "Mac_IGHG4", 
                                                    "Mac_APOE", "Mono_TYMP", "Mono_SNCG"))

pd <- new("AnnotatedDataFrame", data = Mye_traj_seob@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(Mye_traj_seob), row.names = rownames(Mye_traj_seob)))
cds <- newCellDataSet(as.matrix(GetAssayData(Mye_traj_seob, layer = "counts")), 
                      phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds) %>% estimateDispersions()
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)$gene_id

cds <- setOrderingFilter(cds, ordering_genes) %>%
  reduceDimension(max_components = 2, method = 'DDRTree') %>%
  orderCells(root_state = 3)

# --- 4.1 手搓拟时序基因表达平滑热图 ---
cds_DGT <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
peu_gene <- cds_DGT %>% 
  filter(qval < 0.01 & num_cells_expressed > 100) %>%
  filter(!str_starts(gene_short_name, "^MT-|^RPL|^RPS")) %>%
  arrange(qval) %>% head(100)

newdata <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime), max(pData(cds)$Pseudotime), length.out = 100))
m <- genSmoothCurves(cds[peu_gene$gene_short_name,], trend_formula = '~sm.ns(Pseudotime, df=3)', 
                     relative_expr = TRUE, new_data = newdata)
m <- m[!apply(m, 1, sum) == 0, ]
m <- log10(m + 1)
m <- m[!apply(m, 1, sd) == 0, ]
m <- Matrix::t(scale(Matrix::t(m), center = TRUE))
m[m > 3] <- 3; m[m < -3] <- -3

callback <- function(hc, mat){
  sv <- svd(t(mat))$v[,1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

annotation_col <- data.frame(pseudotime = scales::rescale(newdata$Pseudotime, to = c(-1, 1)))
row.names(annotation_col) <- colnames(m)

pheatmap(m, useRaster = TRUE, cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = TRUE, show_colnames = FALSE, clustering_method = "ward.D2", 
         cutree_rows = 3, annotation_col = annotation_col, clustering_callback = callback,
         color = colorRampPalette(c("navy","white","firebrick3"))(100), main = "Pseudotime")

# ==============================================================================
# 5. SCENIC 转录因子网络分析 (SCENIC Analysis)
# ==============================================================================
# 提取已通过 pySCENIC 跑完的 loom 文件
sce_SCENIC <- open_loom("sce_SCENIC.loom")
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name = 'RegulonsAUC')
sub_regulonAUC <- regulonAUC[, match(colnames(Mye_traj_seob), colnames(regulonAUC))]
Mye_traj_seob@meta.data <- cbind(Mye_traj_seob@meta.data, t(assay(sub_regulonAUC[regulonAUC@NAMES,])))

# --- 5.1 RSS 特异性分析与 ggheatmap ---
cellTypes <- data.frame(row.names = colnames(Mye_traj_seob), celltype = Mye_traj_seob$cell_type)
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = cellTypes[colnames(sub_regulonAUC), "celltype"])
rss <- na.omit(rss)

rssPlot <- plotRSS(rss, zThreshold = 3, cluster_columns = FALSE, order_rows = TRUE, 
                   thr = 0.1, varName = "cellType", 
                   col.low = '#330066', col.mid = '#66CC66', col.high = '#FFCC33')

# 使用 ggheatmap 绘制特异性转录因子
rss_data <- reshape2::dcast(rssPlot$plot$data, Topic ~ cellType, value.var = 'Z')
rownames(rss_data) <- rss_data[,1]; rss_data <- rss_data[,-1]

col_ann <- data.frame(group = colnames(rss_data))
rownames(col_ann) <- colnames(rss_data)
ggheatmap(rss_data, color = colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
          cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
          annotation_cols = col_ann, legendName = "Relative value")

# --- 5.2 网络图构建 (Network Graph) ---
# 提示: CEBPB_gene.csv 和 STAT3_gene.csv 需通过 Python/R 字符串清洗从 sce.regulons.csv 中提取
CEBPB_gene <- read.csv("CEBPB_gene.csv", header = TRUE)
STAT3_gene <- read.csv("STAT3_gene.csv", header = TRUE)
TF_target <- rbind(CEBPB_gene, STAT3_gene)
TF_target$score <- as.numeric(TF_target$score)

nodes <- rbind(
  data.frame(name = unique(union(subset(TF_target, tf=="CEBPB")$tf, subset(TF_target, tf=="CEBPB")$target)), value = c(1, subset(TF_target, tf=="CEBPB")$score)),
  data.frame(name = unique(union(subset(TF_target, tf=="STAT3")$tf, subset(TF_target, tf=="STAT3")$target)), value = c(1, subset(TF_target, tf=="STAT3")$score))
)
nodes$cluster <- c(rep("CEBPB",1), rep("CEBPB_gene", nrow(subset(TF_target, tf=="CEBPB"))), 
                   rep("STAT3",1), rep("STAT3_gene", nrow(subset(TF_target, tf=="STAT3"))))
edges <- TF_target[, c("tf", "target", "score")]
edges$class <- edges$tf

layout_cir <- tbl_graph(nodes = nodes, edges = edges)

ggraph(layout_cir, layout = "linear", circular = TRUE) + 
  geom_node_point(aes(size = value, colour = cluster)) +
  geom_node_text(aes(x = 1.03 * x, y = 1.03 * y, label = name, color = cluster,
                     angle = -((-node_angle(x, y) + 90 ) %% 180) + 90), hjust = "outward") +
  geom_edge_arc(aes(colour = class)) + 
  theme_void() + theme(legend.position = "none") +
  scale_colour_manual(values = c("#407972", "#961E28", "#D46724", "#0f8096")) +
  scale_edge_colour_manual(values = c("#961E28", "#D46724", "#0f8096")) +
  scale_size_continuous(range = c(2,8)) + 
  coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.2,1.2))