# 03. CD4 T cells: removal of low-quality clusters, clustering, tumor subset
#
# Processing record, as run interactively on the institutional HPC (R 4.1.1, Seurat 4.1.0) in 2023.
# Source: "20230512 radiation analysis.R", lines 1393-1474, 2717-2731. The code is unchanged except that
#   (i) absolute HPC paths are replaced by file.path(data_dir, ...) / file.path(upstream_dir, ...),
#   (ii) setwd() and install.packages() calls are commented out,
#   (iii) console output and shell lines that had been pasted into the script are commented out.
# It documents how the objects used by the figure scripts were derived; it was not re-executed
# for this deposit. Exploratory plotting calls between the processing steps are kept so that the
# record of each decision (resolution, clusters removed) stays next to the step it justifies.
# Input : IR_CD4T.RData
# Output: IR_CD4T.RData (updated), IR_CD4T_tumor.RData (input of figures/CD4_T_cells.R)
# Note  : in the source, save(IR_CD4T_tumor) is written before the tumor-level re-clustering. The stored
#         IR_CD4T_tumor.RData does contain that re-clustering (integrated_snn_res.0.4 with six clusters
#         and the tumor UMAP), so the object was saved again after it; that save is added at the end.

source("config.R")
# upstream_dir: folder with the inputs produced before this script (see README, "Upstream inputs")
upstream_dir <- Sys.getenv("SCRNA_UPSTREAM_DIR", unset = "upstream")

# ---- original lines 1393-1474 ----------------------------------------
#CD4 T cells
load(file.path(data_dir, "IR_CD4T.RData"))
#remove low quality clusters
DimPlot(IR_CD4T, reduction = "umap", label = T)
VlnPlot(IR_CD4T, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
IR_CD4T = subset(IR_CD4T, idents = c('cluster27_CD4', 'cluster18_cd4', '10', '14', '2', '24', '9'))

#clustering and annotation
DefaultAssay(IR_CD4T) <- "integrated"
ElbowPlot(IR_CD4T)
IR_CD4T <- FindNeighbors(IR_CD4T, dims = 1:20)
IR_CD4T <- FindClusters(object = IR_CD4T,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_CD4T <- RunUMAP(IR_CD4T, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_CD4T) <- "RNA"
Idents(IR_CD4T) = IR_CD4T@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_CD4T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD4T) = IR_CD4T@meta.data$integrated_snn_res.0.6
DimPlot(subset(IR_CD4T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD4T) = IR_CD4T@meta.data$integrated_snn_res.0.8
DimPlot(subset(IR_CD4T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD4T) = IR_CD4T@meta.data$integrated_snn_res.1
DimPlot(subset(IR_CD4T, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(IR_CD4T) = IR_CD4T@meta.data$integrated_snn_res.1.4
DimPlot(subset(IR_CD4T, downsample = 2000), label.size = 4,repel = T,label = T)

#Use 1 as resolution
Idents(object = IR_CD4T) = IR_CD4T@meta.data$integrated_snn_res.1
IR_CD4T[["Yan_annotate"]] <- Idents(object = IR_CD4T)
IR_CD4T$integrated_snn_res.1 <- NULL
IR_CD4T$integrated_snn_res.1.4 <- NULL
IR_CD4T$integrated_snn_res.0.4 <- NULL
IR_CD4T$integrated_snn_res.0.6 <- NULL
IR_CD4T$integrated_snn_res.0.8 <- NULL
IR_CD4T$seurat_clusters <- NULL

Idents(IR_CD4T) = IR_CD4T@meta.data$Yan_annotate
DimPlot(subset(IR_CD4T, idents = c("0","1","2")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD4T, idents = c("3","4", "5")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD4T, idents = c("6","7","8")), label.size = 4,repel = T,label = T)
DimPlot(subset(IR_CD4T, idents = c("9","10", "11")), label.size = 4,repel = T,label = T)
VlnPlot(IR_CD4T, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))

IR_CD4T = subset(IR_CD4T, idents = c(0,1,2,3,4,5,6,7,8,9,10))

#Find markers
DefaultAssay(IR_CD4T) <- "RNA"
IR_CD4T <- NormalizeData(IR_CD4T)
IR_CD4T <- FindVariableFeatures(IR_CD4T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_CD4T)
IR_CD4T <- ScaleData(IR_CD4T, features = all.genes)
all.markers <- FindAllMarkers(IR_CD4T, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "cd4Tcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_CD4T, features = top5_markers$gene) + NoLegend()

#Annotate T cells
DotPlot(IR_CD4T, features = c("Cd3d","Cd4","Cd8b1","Il2ra","Foxp3","Ikzf2")) + RotatedAxis()
DotPlot(IR_CD4T, features = c("Cd69","Cd44","Sell","Cd320","Il2ra","Ctla4","Foxp3","Ikzf2", "Tbx21", "Ifng", "Tnf","Il2", "Gata3", "Il4", "Il5", "Il13","Spi1","Il9","Rorc", "Il17a","Il22a","Bcl6","Il21")) + RotatedAxis()
DotPlot(IR_CD4T, features = c("Ifng", "Tnf","Il2")) + RotatedAxis()
DotPlot(IR_CD4T, features = c("Tnfrsf4", "Tnfrsf9","Ctla4","Nkg7")) + RotatedAxis()
DimPlot(IR_CD4T, reduction = "umap", label = T, split.by = 'orig.ident' )
## Find unique markers for each Treg cluster
markers <- FindMarkers(IR_CD4T, ident.1 = c(6), ident.2 = c(7,10),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(7), ident.2 = c(6,10),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(10), ident.2 = c(6,7),only.pos = T)
## Find markers for naive cd4 t
markers <- FindMarkers(IR_CD4T, ident.1 = c(0), ident.2 = c(1,2,3,4,5,9),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(1), ident.2 = c(0,2,3,4,5,9),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(2), ident.2 = c(0,1,3,4,5,9),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(3), ident.2 = c(0,1,2,4,5,9),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(4), ident.2 = c(0,1,2,3,5,9),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(5), ident.2 = c(0,1,2,3,4,9),only.pos = T)
markers <- FindMarkers(IR_CD4T, ident.1 = c(9), ident.2 = c(0,1,2,3,4,5),only.pos = T)

save(IR_CD4T, file = file.path(data_dir, "IR_CD4T.RData"))

# ---- original lines 2717-2731 ----------------------------------------
load(file.path(data_dir, "IR_CD4T.RData"))
Idents(IR_CD4T) = IR_CD4T@meta.data$orig.ident
IR_CD4T_tumor <- subset(IR_CD4T, idents = c("local_ctrl_T","distal_ctrl_T","local_IR_T","distal_IR_T"))
save(IR_CD4T_tumor, file = file.path(data_dir, "IR_CD4T_tumor.RData"))

DefaultAssay(IR_CD4T_tumor) <- "integrated"
ElbowPlot(IR_CD4T_tumor)
IR_CD4T_tumor <- FindNeighbors(IR_CD4T_tumor, dims = 1:20)
IR_CD4T_tumor <- FindClusters(object = IR_CD4T_tumor,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_CD4T_tumor <- RunUMAP(IR_CD4T_tumor, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_CD4T_tumor) <- "RNA"
Idents(IR_CD4T_tumor) = IR_CD4T_tumor@meta.data$integrated_snn_res.0.4
DimPlot(IR_CD4T_tumor, label.size = 4,repel = T,label = T)
IR_CD4T_tumor$Yan_annotate = Idents(IR_CD4T_tumor)

# [added for the deposit: see Note in the header]
save(IR_CD4T_tumor, file = file.path(data_dir, "IR_CD4T_tumor.RData"))
