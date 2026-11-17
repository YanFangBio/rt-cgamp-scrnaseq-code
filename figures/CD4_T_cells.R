# Tumor-infiltrating CD4 T cells (scRNA-seq): UMAP, cluster composition per sample, marker heatmap
#
# Input : IR_CD4T_tumor.RData  (tumor CD4 T cells; clustering at resolution 0.4 and the UMAP are
#         stored in the object, see processing/03_CD4_T_cells.R)
# Output: CD4_UMAP.pdf, CD4_cluster_composition.csv / .pdf, CD4_marker_heatmap.pdf
# Figure panels produced by each block: docs/figure_map.csv

source("config.R")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

IR_CD4T_tumor <- load_object("IR_CD4T_tumor")
DefaultAssay(IR_CD4T_tumor) <- "RNA"
Idents(IR_CD4T_tumor) = IR_CD4T_tumor@meta.data$integrated_snn_res.0.4
IR_CD4T_tumor$Yan_annotate = Idents(IR_CD4T_tumor)

# Names given to the six clusters in the figure
cluster_map <- c(
  "0" = "C0-Tnaive",
  "1" = "C1-Treg",
  "2" = "C2-Tint",
  "3" = "C3-Treg",
  "4" = "C4-Th1",
  "5" = "C5-Unknown"
)
IR_CD4T_tumor$state_by_cluster <- factor(
  unname(cluster_map[as.character(IR_CD4T_tumor$Yan_annotate)]),
  levels = unname(cluster_map)
)
Idents(IR_CD4T_tumor) <- "state_by_cluster"

# -----------------------------
# UMAP
# -----------------------------
umap_plot <- DimPlot(IR_CD4T_tumor, reduction = "umap", pt.size = 2, label = TRUE, label.size = 5)
ggsave(filename = file.path(out_dir, "CD4_UMAP.pdf"), plot = umap_plot, width = 10, height = 8)

# -----------------------------
# Cluster composition per tumor sample
# -----------------------------
composition <- as.data.frame(table(cluster = IR_CD4T_tumor$state_by_cluster, sample = IR_CD4T_tumor$orig.ident)) %>%
  group_by(sample) %>%
  mutate(total_CD4T = sum(Freq), percent = 100 * Freq / total_CD4T) %>%
  ungroup() %>%
  rename(n_cells = Freq)
composition$sample <- factor(composition$sample, levels = c("local_ctrl_T", "local_IR_T", "distal_ctrl_T", "distal_IR_T"))
write.csv(composition, file.path(out_dir, "CD4_cluster_composition.csv"), row.names = FALSE)

p_comp <- ggplot(composition, aes(x = sample, y = percent, fill = cluster)) +
  geom_col(width = 0.75) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.title.x = element_blank()) +
  labs(y = "% of CD4 T", fill = NULL)
ggsave(file.path(out_dir, "CD4_cluster_composition.pdf"), p_comp, width = 5, height = 4)

# -----------------------------
# Average expression of marker genes per cluster
# -----------------------------
# Calculate the average expression per cluster
avg_expression <- AverageExpression(IR_CD4T_tumor, return.seurat = TRUE)

cd4gene = c("Ccr7","Sell","Ifng", "Tnf","Il2","Ly6c1","Gm19951", "Igfbp4","Nme1" ,"Rtp4" ,"Isg15","Irf7", "Tgtp1","Igtp","Bach2","Aff3" , "Maml2", "Smc4","Tnfrsf4","Tnfrsf9" ,"Ctla4","S100a6", "Tigit" ,"Il2ra","Foxp3","Ikzf2","Il2rb", "Il10", "Izumo1r","Id3" , "Dusp2","Cd69"  ,"Ifit1" ,"Ifit3" ,"Ifit3b","Bst2" ,"S100a4" , "Itgb1", "Id2" ,"Nkg7", "Tcf7" )
avg_exp_matrix <- as.matrix(get_assay_matrix(avg_expression, "RNA", "data")[cd4gene, ])
new_cluster_order <- unname(cluster_map[c("0", "1", "2", "3", "5", "4")])
avg_exp_matrix_ordered <- avg_exp_matrix[, new_cluster_order]

pdf(file.path(out_dir, "CD4_marker_heatmap.pdf"))
draw(Heatmap(avg_exp_matrix_ordered,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = TRUE,
        column_names_rot = 45,
        column_gap = unit(0, "mm")))
dev.off()
