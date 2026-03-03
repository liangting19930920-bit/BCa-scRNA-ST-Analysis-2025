# ==============================================================================
# Project: Bladder Cancer Single-Cell Analysis (ScBLCA)
# Description: Virtual Knockout Analysis of ITGA5 and ITGB1 using scTenifoldKnk
# ==============================================================================

# ==============================================================================
# 0. 环境准备与依赖加载 (Dependencies)
# ==============================================================================
# 提示: scTenifoldKnk 可通过 devtools::install_github("cailab-tamu/scTenifoldKnk") 安装
library(Seurat)
library(scTenifoldKnk)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

# 全局路径配置
WORK_DIR <- "data/"
OUT_DIR <- "results/Virtual_KO/"
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
setwd(WORK_DIR)

# ==============================================================================
# 1. 细胞下采样与基因严格过滤 (Data Preparation & Filtering)
# ==============================================================================
message(">>> Preparing Data for Virtual Knockout...")


# 为了控制计算资源和结果可重复性，进行随机下采样
set.seed(123) 
cells_to_keep <- sample(colnames(seob2), min(6000, ncol(seob2)))
target_seob <- subset(seob2, cells = cells_to_keep)

# 基因过滤策略：
# 1. 至少在 10 个细胞中表达 (counts > 0)
# 2. 细胞群体中的平均表达量 > 0.1
counts_mat <- target_seob@assays$RNA@counts
genes_in_min_cells <- rowSums(counts_mat > 0) >= 10
genes_avg_exp <- Matrix::rowMeans(counts_mat) > 0.1

genes_to_keep <- genes_in_min_cells & genes_avg_exp
target_seob <- subset(target_seob, features = rownames(target_seob)[genes_to_keep])

# 检查目标基因是否通过过滤并在矩阵中
target_genes <- c("ITGA5", "ITGB1")
if(!all(target_genes %in% rownames(target_seob))) {
  stop("Target genes for KO are not present after filtering!")
}

# 提取最终 Count 矩阵
countMatrix <- as.matrix(target_seob@assays$RNA@counts)

# ==============================================================================
# 2. ITGA5 / ITGB1 联合虚拟敲除 (scTenifoldKnk Execution)
# ==============================================================================
message(paste(">>> Running scTenifoldKnk for:", paste(target_genes, collapse = " & ")))



# 构建单细胞基因调控网络 (scGRN) 并执行流形对齐评分
result_ITG_KO <- scTenifoldKnk(
  countMatrix = countMatrix, 
  gKO = target_genes,          # 联合敲除 ITGA5 和 ITGB1
  qc = FALSE,                  # 已经在 Seurat 中做过 QC
  nc_nNet = 10,                # 核心参数: 子网络数量
  nc_nCells = 500,             # 核心参数: 每个网络中随机抽取的细胞数
  nc_nComp = 3                 # 核心参数: PCA 的主成分数量
)

# 保存虚拟敲除原始结果
saveRDS(result_ITG_KO, file = paste0(OUT_DIR, "scTenifoldKnk_ITGA5_ITGB1_Result.rds"))

# ==============================================================================
# 3. 敲除效应扰动分析可视化 (Perturbation Visualization)
# ==============================================================================
message(">>> Visualizing Differential Regulation...")

df_res <- result_ITG_KO$diffRegulation
df_res$log_pval <- -log10(df_res$p.adj)

# --- 3.1 Top 20 受影响基因柱状图 ---
# 提取 FC 排名靠前的基因，并剔除被敲除的基因本身
top_genes <- head(df_res[order(-df_res$FC), ], 30)
top_genes_filtered <- top_genes[!top_genes$gene %in% target_genes, ] %>% head(20)

p_bar <- ggplot(top_genes_filtered, aes(x = reorder(gene, FC), y = FC)) +
  geom_bar(stat = 'identity', fill = '#4b5cc4', width = 0.7) +
  coord_flip() +
  labs(title = "Top 20 Up-Regulated Genes Post-KO", x = "Gene", y = "Fold Change (FC)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# --- 3.2 Z-score vs p-value 火山图 ---
# 标记显著扰动的基因 (Z-score 绝对值 > 2 且 p.adj < 0.01)
label_genes <- subset(df_res, abs(Z) > 2 & p.adj < 0.01 & !gene %in% target_genes)



p_volcano <- ggplot(df_res, aes(x = Z, y = log_pval)) +
  geom_point(alpha = 0.5, color = "grey50") +
  geom_point(data = label_genes, aes(x = Z, y = log_pval), color = "#d83215", alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") + 
  geom_text_repel(data = label_genes, aes(label = gene), size = 4, max.overlaps = 50) +
  labs(title = "Virtual KO: ITGA5 & ITGB1", x = "Z-score", y = "-log10(p-value)") +
  theme_classic(base_size = 14) + 
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 14, color = "black"))

# ==============================================================================
# 4. 扰动下游通路 GO 富集分析 (Downstream GO Enrichment)
# ==============================================================================
message(">>> Running GO Enrichment on highly perturbed genes...")

# 提取差异最显著的前 500 个基因进行富集
top_genes500 <- head(df_res[order(-df_res$FC), "gene"], 500)

go_res <- enrichGO(gene = top_genes500, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                   ont = "ALL", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

dt_go <- as.data.frame(go_res)
dt_go$LogP <- -log10(as.numeric(dt_go$p.adjust))
dt_go$Description <- factor(dt_go$Description, levels = dt_go$Description[order(dt_go$LogP)])

# 保存富集结果表
# write.csv(dt_go, file = paste0(OUT_DIR, "GO_Enrichment_ITG_KO.csv"))



# --- 4.1 顶级 GO 梯级颜色条形图 ---
p_go <- ggplot(head(dt_go, 20), aes(x = Description, y = Count, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +  
  labs(title = "Enriched Pathways Altered by KO", x = "", y = "Gene Counts", fill = "-Log10(P)") +
  theme_test(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, vjust = 1),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_fill_gradient2(low = "#486b98", mid = "#f5f2b1", high = "#b93735", 
                       midpoint = median(head(dt_go, 20)$LogP))  

print(p_bar | p_volcano)
print(p_go)