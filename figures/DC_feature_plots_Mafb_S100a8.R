# Mafb and S100a8 expression in the DC-enriched scRNA-seq dataset (macrophage / neutrophil contamination check)
#
# Input : IR_myeloid_DC.RData
# Output: DC_Mafb_S100a8.pdf
# Figure panel: docs/figure_map.csv

source("config.R")
library(Seurat)
library(ggplot2)

IR_myeloid_DC <- load_object("IR_myeloid_DC")
DefaultAssay(IR_myeloid_DC) <- "RNA"
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$Yan_annotate

p <- FeaturePlot(IR_myeloid_DC, features = c("Mafb", "S100a8"))
ggsave(file.path(out_dir, "DC_Mafb_S100a8.pdf"), p, width = 8, height = 3.5)
