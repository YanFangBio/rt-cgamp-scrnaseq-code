# 02. CD8 T cells: clustering, removal of low-quality cells, tumor subset
#
# Processing record, as run interactively on the institutional HPC (R 4.1.1, Seurat 4.1.0) in 2023.
# Source: "20230512 radiation analysis.R", lines 285-349, 957-1048, 2469-2508, 2570-2570. The code is unchanged except that
#   (i) absolute HPC paths are replaced by file.path(data_dir, ...) / file.path(upstream_dir, ...),
#   (ii) setwd() and install.packages() calls are commented out,
#   (iii) console output and shell lines that had been pasted into the script are commented out.
# It documents how the objects used by the figure scripts were derived; it was not re-executed
# for this deposit. Exploratory plotting calls between the processing steps are kept so that the
# record of each decision (resolution, clusters removed) stays next to the step it justifies.
# Input : IR_CD8T.RData
# Output: IR_CD8T.RData (updated), IR_CD8T_tumor.RData (input of figures/CD8_T_cell_states.R and figures/CD8_composition_TCR_expansion.R)
# Note  : lines 350-956 and 1049-2468 of the source (annotation of the LN + tumor object, TCR clone
#         analyses, other lineages) are not part of the path to IR_CD8T_tumor and are omitted.

source("config.R")
# upstream_dir: folder with the inputs produced before this script (see README, "Upstream inputs")
upstream_dir <- Sys.getenv("SCRNA_UPSTREAM_DIR", unset = "upstream")

# ---- original lines 285-349 ----------------------------------------
load(file.path(data_dir, "IR_CD8T.RData"))
#CD8 T cells clustering and annotation
DefaultAssay(IR_CD8T) <- "integrated"
ElbowPlot(IR_CD8T)
IR_CD8T <- FindNeighbors(IR_CD8T, dims = 1:20)
IR_CD8T <- FindClusters(object = IR_CD8T,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_CD8T <- RunUMAP(IR_CD8T, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_CD8T) <- "RNA"
Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.0.6
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.0.8
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.1
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.1.4
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

#Use 0.8 as resolution
Idents(object = IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.1
IR_CD8T[["Yan_annotate"]] <- Idents(object = IR_CD8T)
IR_CD8T$integrated_snn_res.1 <- NULL
IR_CD8T$integrated_snn_res.1.4 <- NULL
IR_CD8T$integrated_snn_res.0.4 <- NULL
IR_CD8T$integrated_snn_res.0.6 <- NULL
IR_CD8T$integrated_snn_res.0.8 <- NULL
IR_CD8T$seurat_clusters <- NULL

#Use resolution 0.8 ans double check the clustering
Idents(IR_CD8T) = IR_CD8T@meta.data$Yan_annotate
DimPlot(subset(IR_CD8T, idents = c("0","1","2")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD8T, idents = c("3","4", "5")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD8T, idents = c("6","7","8")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD8T, idents = c("9","10", "11","12")), label.size = 4,repel = T,label = T)
VlnPlot(IR_CD8T, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
DotPlot(IR_CD8T, features = c("Clec10a", "Folr2", "Adgre1")) + RotatedAxis()
DotPlot(IR_CD8T, features = c("Cxcr3", "Cd28", "Cd44")) + RotatedAxis()
#Find markers
DefaultAssay(IR_CD8T) <- "RNA"
IR_CD8T <- NormalizeData(IR_CD8T)
IR_CD8T <- FindVariableFeatures(IR_CD8T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_CD8T)
IR_CD8T <- ScaleData(IR_CD8T, features = all.genes)
all.markers <- FindAllMarkers(IR_CD8T, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "Tcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_CD8T, features = top5_markers$gene) + NoLegend()

library(ggplot2)
# install.packages("svglite")

umap_plot <- DimPlot(IR_CD8T, reduction = "umap", pt.size =2)
ggsave(filename = file.path(data_dir, "CD8umap_plot.svg"), plot = umap_plot, width = 10, height = 8)
umap_plot <- DimPlot(IR_CD8T, reduction = "umap", pt.size =2, label = TRUE, label.size = 5, label.box = TRUE)
ggsave(filename = file.path(data_dir, "CD8umap_plot_label.svg"), plot = umap_plot, width = 10, height = 8)



# ---- original lines 957-1048 ----------------------------------------
#-----------------------------2023/05/30---------------------------------
#remove low quality dying cells
Idents(IR_CD8T) = IR_CD8T@meta.data$Yan_annotate
IR_CD8T = subset(IR_CD8T, idents = c(0,2,3,4,5,6,8,9,11,12))
#CD8 T cells clustering and annotation
DefaultAssay(IR_CD8T) <- "integrated"
ElbowPlot(IR_CD8T)
IR_CD8T <- FindNeighbors(IR_CD8T, dims = 1:20)
IR_CD8T <- FindClusters(object = IR_CD8T,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_CD8T <- RunUMAP(IR_CD8T, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_CD8T) <- "RNA"
Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.0.6
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.0.8
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.1
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.1.4
DimPlot(subset(IR_CD8T, downsample = 2000), label.size = 4,repel = T,label = T)
#Use 0.8 as resolution

Idents(object = IR_CD8T) = IR_CD8T@meta.data$integrated_snn_res.0.8
IR_CD8T[["Yan_annotate"]] <- Idents(object = IR_CD8T)
IR_CD8T$integrated_snn_res.1 <- NULL
IR_CD8T$integrated_snn_res.1.4 <- NULL
IR_CD8T$integrated_snn_res.0.4 <- NULL
IR_CD8T$integrated_snn_res.0.6 <- NULL
IR_CD8T$integrated_snn_res.0.8 <- NULL
IR_CD8T$seurat_clusters <- NULL

#Use resolution 0.8 ans double check the clustering
Idents(IR_CD8T) = IR_CD8T@meta.data$Yan_annotate
DimPlot(subset(IR_CD8T, idents = c("0","1","2")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD8T, idents = c("3","4", "5")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD8T, idents = c("6","7","8")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD8T, idents = c("9","10", "11")), label.size = 4,repel = T,label = T)
VlnPlot(IR_CD8T, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
DotPlot(IR_CD8T, features = c("Clec10a", "Folr2", "Adgre1")) + RotatedAxis()

#Find markers
DefaultAssay(IR_CD8T) <- "RNA"
IR_CD8T <- NormalizeData(IR_CD8T)
IR_CD8T <- FindVariableFeatures(IR_CD8T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_CD8T)
IR_CD8T <- ScaleData(IR_CD8T, features = all.genes)
all.markers <- FindAllMarkers(IR_CD8T, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "Tcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_CD8T, features = top5_markers$gene) + NoLegend()

#Annotate T cells
##CD8 T cells
CD8_marker = c("Tcf7","Sell","Lef1","Ccr7","Il7r", "Gzma","Gzmb","Gzmk","Cxcr5","Ltb","Cx3cr1","Prf1","Gnly","Fgfbp2", "Klrg1", "Ifng","Ccl4", "Tim2", "Itgae", "Havcr2", "Entpd1", "Layn","Pdcd1","Ctla4","Lag3", "Tigit")
DotPlot(IR_CD8T, features = CD8_marker) + RotatedAxis()
CD8_marker = c("Fasl" ,"Gzma","Gzmb","Gzmk", "Klrg1", "Ifng","Ccl4")
#Naive 0,2,3,4,5,6,10
markers0 <- FindMarkers(IR_CD8T, ident.1 = c(0), ident.2 = c(2,3,4,5,6,10),only.pos = T)
markers2 <- FindMarkers(IR_CD8T, ident.1 = c(2), ident.2 = c(0,3,4,5,6,10),only.pos = T)
markers3 <- FindMarkers(IR_CD8T, ident.1 = c(3), ident.2 = c(0,2,4,5,6,10),only.pos = T)
markers4 <- FindMarkers(IR_CD8T, ident.1 = c(4), ident.2 = c(0,2,3,5,6,10),only.pos = T)
markers5 <- FindMarkers(IR_CD8T, ident.1 = c(5), ident.2 = c(0,2,3,4,6,10),only.pos = T)
markers6 <- FindMarkers(IR_CD8T, ident.1 = c(6), ident.2 = c(0,2,3,4,5,10),only.pos = T)
markers10 <- FindMarkers(IR_CD8T, ident.1 = c(10), ident.2 = c(0,2,3,4,5,6),only.pos = T)
markers1 <- FindMarkers(IR_CD8T, ident.1 = c(1), ident.2 = c(0,2,3,4,5,6,10),only.pos = T)
markers0$gene <- rownames(markers0)
markers2$gene <- rownames(markers2)
markers3$gene <- rownames(markers3)
markers4$gene <- rownames(markers4)
markers5$gene <- rownames(markers5)
markers6$gene <- rownames(markers6)
markers10$gene <- rownames(markers10)
markers1$gene <- rownames(markers1)
library("writexl")
write_xlsx(markers0, file.path(data_dir, "markers0.xlsx"))
write_xlsx(markers2, file.path(data_dir, "markers2.xlsx"))
write_xlsx(markers3, file.path(data_dir, "markers3.xlsx"))
write_xlsx(markers4, file.path(data_dir, "markers4.xlsx"))
write_xlsx(markers5, file.path(data_dir, "markers5.xlsx"))
write_xlsx(markers6, file.path(data_dir, "markers6.xlsx"))
write_xlsx(markers10, file.path(data_dir, "markers10.xlsx"))
write_xlsx(markers1, file.path(data_dir, "markers1.xlsx"))

save(IR_CD8T, file = file.path(data_dir, "IR_CD8T.RData"))

# ---- original lines 2469-2508 ----------------------------------------
load(file.path(data_dir, "IR_CD8T.RData"))
Idents(IR_CD8T) = IR_CD8T@meta.data$orig.ident
IR_CD8T_tumor <- subset(IR_CD8T, idents = c("local_ctrl_T","distal_ctrl_T","local_IR_T","distal_IR_T"))
IR_CD8T_LN <- subset(IR_CD8T, idents = c("local_ctrl_LN","distal_ctrl_LN","local_IR_LN" ,"distal_IR_LN"))

DefaultAssay(IR_CD8T_tumor) <- "integrated"
ElbowPlot(IR_CD8T_tumor)
IR_CD8T_tumor <- FindNeighbors(IR_CD8T_tumor, dims = 1:20)
IR_CD8T_tumor <- FindClusters(object = IR_CD8T_tumor,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_CD8T_tumor <- RunUMAP(IR_CD8T_tumor, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_CD8T_tumor) <- "RNA"
Idents(IR_CD8T_tumor) = IR_CD8T_tumor@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_CD8T_tumor, downsample = 2000), label.size = 4,repel = T,label = T)

CD8_marker = c("Tcf7","Sell","Lef1","Ccr7","Il7r", "Itgae","Klrg1","Gzma","Gzmb","Gzmk","Prf1", "Ifng","Cd69","Mki67","Cd44","Havcr2", "Entpd1", "Pdcd1","Tox","Ctla4","Lag3", "Tigit","Zfp683","Prdm1","Runx3","Klf2","Cxcr3","Cxcr6")

FeaturePlot(IR_CD8T_tumor, features = CD8_marker, label=T) 
DotPlot(IR_CD8T_tumor, features = CD8_marker) + RotatedAxis()

#Find markers
DefaultAssay(IR_CD8T_tumor) <- "RNA"
IR_CD8T_tumor <- NormalizeData(IR_CD8T_tumor)
IR_CD8T_tumor <- FindVariableFeatures(IR_CD8T_tumor, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_CD8T_tumor)
IR_CD8T_tumor <- ScaleData(IR_CD8T_tumor, features = all.genes)
all.markers <- FindAllMarkers(IR_CD8T_tumor, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "IR_CD8T_tumorTcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_CD8T_tumor, features = top5_markers$gene) + NoLegend()

library(ggplot2)
library(svglite)

umap_plot <- DimPlot(IR_CD8T_tumor, reduction = "umap", pt.size =2)
ggsave(filename = file.path(data_dir, "IR_CD8T_tumorumap_plot.svg"), plot = umap_plot, width = 10, height = 8)
umap_plot <- DimPlot(IR_CD8T_tumor, reduction = "umap", pt.size =2, label = TRUE, label.size = 5, label.box = TRUE)
ggsave(filename = file.path(data_dir, "IR_CD8T_tumorumap_plot_label.svg"), plot = umap_plot, width = 10, height = 8)

# ---- original lines 2570-2570 ----------------------------------------
save(IR_CD8T_tumor, file = file.path(data_dir, "IR_CD8T_tumor.RData"))
