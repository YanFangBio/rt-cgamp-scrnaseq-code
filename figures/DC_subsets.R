# Dendritic-cell subsets in local and distal tumors (scRNA-seq): UMAP, top-10 marker heatmap,
# activation / maturation marker ridge plots, cDC1- and cDC2-associated genes
#
# Input : IR_myeloid_DC.RData  (dendritic cells; annotation cDC1 / cDC2 / CCR7+ in meta.data$Yan_annotate,
#         see processing/04_myeloid_DC.R). The "CCR7+" cluster is named "Mature DC" in the manuscript.
# Output: DC_UMAP.pdf, DC_top10_marker_heatmap.pdf, DC_markers.csv, DC_ridgeplots.pdf,
#         DC_cDC1_features.pdf, DC_cDC2_features.pdf
# Figure panels produced by each block: docs/figure_map.csv

source("config.R")
library(Seurat)
library(dplyr)
library(ggplot2)

IR_myeloid_DC <- load_object("IR_myeloid_DC")
DefaultAssay(IR_myeloid_DC) <- "RNA"
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$Yan_annotate
print(table(IR_myeloid_DC@meta.data$Yan_annotate, IR_myeloid_DC@meta.data$orig.ident))

# -----------------------------
# UMAP
# -----------------------------
pdf(file.path(out_dir, "DC_UMAP.pdf"), width = 4, height = 3)
print(DimPlot(IR_myeloid_DC, label.size = 4, repel = TRUE, label = FALSE))
dev.off()

# -----------------------------
# Top 10 marker genes per DC population
# -----------------------------
# The original call was FindAllMarkers(IR_myeloid_DC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
# under Seurat 4.1.0. `mean.fxn = mean_fxn_v4` (config.R) keeps that version's fold-change definition,
# so the top-10 genes of the published panel are returned by any Seurat version.
all.markers <- FindAllMarkers(IR_myeloid_DC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5, mean.fxn = mean_fxn_v4)
write.csv(all.markers, file.path(out_dir, "DC_markers.csv"), row.names = FALSE)
top10_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
print(top10_markers)

pdf(file.path(out_dir, "DC_top10_marker_heatmap.pdf"), width = 8, height = 6)
print(DoHeatmap(IR_myeloid_DC, features = top10_markers$gene) + NoLegend())
dev.off()

# -----------------------------
# Activation and maturation markers
# -----------------------------
pdf(file.path(out_dir, "DC_ridgeplots.pdf"), width = 12, height = 3)
print(RidgePlot(object = IR_myeloid_DC, features = c("Cd80","Cd274","Ccr7","Il12b","Cd40","Cd86"), ncol=6))
dev.off()

# -----------------------------
# cDC1- and cDC2-associated genes
# -----------------------------
#cDC1
p_cdc1 <- DotPlot(IR_myeloid_DC, features = c("Batf3","Id2","Irf8","Zbtb46", "Btla", "Cadm1","Cd8a","Itgae","Itgax","Ly75","Thbd","Xcr1","Clec9a")) + RotatedAxis()
ggsave(file.path(out_dir, "DC_cDC1_features.pdf"), p_cdc1, width = 6, height = 4)
#cDC2  ("Cd1" was also requested in the original call; it is not a gene in the matrix and Seurat dropped it)
p_cdc2 <- DotPlot(IR_myeloid_DC, features = c("Irf4","Klf4", "Cd14", "Cd163","Clec10a","Notch2","Itgam","Sirpa","Cx3cr1","Cd2")) + RotatedAxis()
ggsave(file.path(out_dir, "DC_cDC2_features.pdf"), p_cdc2, width = 6, height = 4)
