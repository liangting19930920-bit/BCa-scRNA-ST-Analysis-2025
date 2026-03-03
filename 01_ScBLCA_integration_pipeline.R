# ==============================================================================
# Project: LT-ScBLCA (Bladder Cancer Single-Cell Analysis)
# Description: Data Integration, Harmony Correction, Annotation & Visualization
# ==============================================================================

# ==============================================================================
# 1. 加载依赖包 (Dependencies)
# ==============================================================================
library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(ggsci)
library(grid)
library(scRNAtoolVis)
library(cols4all)
library(tidydr)

# ==============================================================================
# 2. 批量构建 Seurat 对象 (Batch Object Construction)
# ==============================================================================
sample_meta <- tribble(
  ~sample_id, ~data_dir, ~type, ~status, ~format,
  "BC5",   "path/GSM4006648_BC5_gene_cell_exprs_table.txt", "BC_H1",  "High_Grade", "txt",
  "BC6",   "path/GSM4751267_BC6_gene_cell_exprs_table.txt", "BC_L1",  "Low_Grade",  "txt",
  "BC7",   "path/GSM4751268_BC7_gene_cell_exprs_table.txt", "BC_L2",  "Low_Grade",  "txt",
  "BCN_1", "path/GSM5329919_BCN_gene_cell_exprs_table.txt", "BCN_1",  "NAT",        "txt",  
  "BC_H11","path/BC_H11_Sc",                     "BC_H11", "High_Grade", "10X"
  # ... [此处按照上方格式补充完你的 25 个样本] ...
  )

load_sample <- function(path, format) {
  if (format == "txt") {
    return(read.table(path, row.names = 1, header = TRUE))
  } else {
    return(Read10X(data.dir = path))
  }
}

seob_list <- apply(sample_meta, 1, function(row) {
  message("Loading: ", row["sample_id"])
  counts <- load_sample(row["data_dir"], row["format"])
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
  obj$type <- row["type"]
  obj$status <- row["status"] 
  return(obj)
})

names(seob_list) <- sample_meta$sample_id
seob <- merge(x = seob_list[[1]], y = seob_list[-1], add.cell.ids = names(seob_list))
rm(seob_list); gc()

seob$type <- ifelse(seob$type == "GC", "GC_1", seob$type)

# ==============================================================================
# 3. 质控与过滤 (Quality Control)
# ==============================================================================
seob[["percent.mt"]] <- PercentageFeatureSet(seob, pattern = "^MT-")
seob <- CellCycleScoring(seob, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
seob <- subset(seob, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 50)

# ==============================================================================
# 4. 降维、聚类与批次校正 (Dimensionality Reduction & Harmony)
# ==============================================================================
seob <- seob %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunHarmony(reduction = "pca", group.by.vars = "type", reduction.save = "harmony") %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4, random.seed = 1)

# ==============================================================================
# 5. 细胞注释 (Cell Type Annotation)
# ==============================================================================
cluster_map <- list(
  "B cells"           = c('25','5','27'),
  "Epithelial cells"  = c('29','8','18','20','19','4','30','9','13','16','14','26','1','0','15','21'),
  "T/NK cells"        = c('7','2','28'),
  "Myeloid cells"     = c('12'),
  "Fibroblast cells"  = c('6','10','24'),
  "Endothelial cells" = c('3','23'),
  "Plasma cells"      = c('11','22','17')
)

annotation_vector <- unlist(lapply(names(cluster_map), function(x) setNames(rep(x, length(cluster_map[[x]])), cluster_map[[x]])))
seob$cell_type <- annotation_vector[as.character(seob$seurat_clusters)]
Idents(seob) <- "cell_type"

# ==============================================================================
# 6. 数据可视化 (Advanced Visualization)
# ==============================================================================

# --- 6.1 经典 Marker 基因气泡图 ---
list_genes <- list(
  T_NK_cells        = c("CD3D",'CD3E',"CD2",'CD4','CD8A','KLRB1','NCR1'),
  B_cells           = c("MS4A1","CD79A","CD79B","CD19"),
  Plasma_cells      = c("MZB1","IGHG1","IGHG3","IGLC2","IGKC"),
  Myeloid_cells     = c("LYZ","CD14","FCN1","CD68","C1QB","CD163","CD1E","LAMP3"),
  Fibroblast_cells  = c("COL1A1","COL3A1","DCN","FGF7"),
  Epithelial_cells  = c("EPCAM","KRT13","KRT18","KRT7"),
  Endothelial_cells = c("PECAM1","PLVAP","VWF","CLDN5")
)

p1 <- DotPlot(seob, features=list_genes, cols = c("grey", "red"), cluster.idents = TRUE) +
  RotatedAxis() +
  theme(
    panel.border = element_rect(color="black"),
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outside',
    axis.line = element_blank()
  ) + labs(x="", y="")

# 气泡图增加分界线
p_grob <- ggplotGrob(p1)
lg <- linesGrob(x=unit(c(0,1),"npc"), y=unit(c(0,0)+0.2,"npc"), gp=gpar(col="black", lwd=4))
grid.newpage()
for (k in grep("strip-t", p_grob$layout$name)) {
  p_grob$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}
grid.draw(p_grob)

# --- 6.2 提取 UMAP 坐标用于高级 ggplot 可视化 ---
UMAP_df <- cbind(
  as.data.frame(Embeddings(seob, "umap")),
  seob@meta.data[, c("seurat_clusters", "cell_type", "type", "status")]
)

mytheme <- theme_void() + theme(plot.margin = margin(5.5,15,5.5,5.5)) 
mycol <- c4a('wright25', 15)

# 高级 UMAP
p_umap_custom <- ggplot(data = UMAP_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cell_type), size = 0.5, alpha = 0.8) + 
  stat_ellipse(aes(color = cell_type, fill = cell_type),
               level = 0.95, linetype = 2, show.legend = FALSE, geom = 'polygon', alpha = 0.1) +
  mytheme + theme_dr(xlength = 0.2, ylength = 0.2, 
                     arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) + 
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol)

# --- 6.3 细胞比例频率分布 ---
Cellratio <- prop.table(table(Idents(seob), seob@meta.data$type), margin = 2) %>% as.data.frame()
ggplot(Cellratio) + 
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = 'Type', y = 'Ratio') +
  coord_flip() +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"))

# --- 6.4 自定义 Grid 布局散点图+条形图 (完整无删减版) ---
reduc <- data.frame(Seurat::Embeddings(seob, reduction = "umap"))
meta <- seob@meta.data
pc12 <- cbind(reduc, meta)
pc12$idents <- Idents(seob)
celltype <- unique(pc12$idents)

# 确保配色严格对应到具体的细胞亚群
cols <- c(
  'Epithelial cells'  = '#6EB467',
  'T/NK cells'        = '#80689C',
  'Fibroblast cells'  = '#CEB481',
  'Endothelial cells' = '#D86F6F',
  'B cells'           = '#1965B0',
  'Plasma cells'      = '#A6761D',
  'Myeloid cells'     = '#E31A1C'
)
cols <- cols[celltype]

# 统计各个细胞群的数量用于绘制条形图
cell_num <- pc12 %>%
  dplyr::group_by(idents) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::arrange(n)

# 绘制条形图框架
grid.newpage()
pushViewport(viewport(x = unit(0.61, "npc"), y = unit(0.5, "npc"),
                      width = unit(0.2, "npc"), height = unit(0.7, "npc"),
                      just = "left",
                      yscale = c(0, nrow(cell_num) + 0.5),
                      xscale = c(0, max(cell_num$n) + 0.1 * max(cell_num$n))))

grid.xaxis()
grid.rect(x = rep(0, nrow(cell_num)), y = 1:nrow(cell_num),
          width = cell_num$n, height = unit(0.08, "npc"),
          just = "left", gp = gpar(fill = cols, col = NA), default.units = "native")
grid.rect(gp = gpar(fill = "transparent"))
grid.text(label = "Number of cells", x = 0.5, y = unit(-2.5, "lines"))
popViewport()

# 绘制图例和文字标签
pushViewport(viewport(x = unit(0.81, "npc"), y = unit(0.5, "npc"),
                      width = unit(0.2, "npc"), height = unit(0.7, "npc"),
                      just = "left", yscale = c(0, nrow(cell_num) + 0.5)))

grid.points(x = rep(0.1, nrow(cell_num)), y = 1:nrow(cell_num), pch = 19,
            gp = gpar(col = cols), size = unit(1.5, "char"))

grid.text(label = cell_num$idents, x = 0.2, y = 1:nrow(cell_num),
          just = "left", default.units = "native")
popViewport()

