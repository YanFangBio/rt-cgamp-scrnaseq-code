# 05. NK cells: removal of low-quality cells and of a CD8 T-cell cluster, clustering
#
# Processing record, as run interactively on the institutional HPC (R 4.1.1, Seurat 4.1.0) in 2023.
# Source: "20230512 radiation analysis.R", lines 1914-2011. The code is unchanged except that
#   (i) absolute HPC paths are replaced by file.path(data_dir, ...) / file.path(upstream_dir, ...),
#   (ii) setwd() and install.packages() calls are commented out,
#   (iii) console output and shell lines that had been pasted into the script are commented out.
# It documents how the objects used by the figure scripts were derived; it was not re-executed
# for this deposit. Exploratory plotting calls between the processing steps are kept so that the
# record of each decision (resolution, clusters removed) stays next to the step it justifies.
# Input : IR_NK.RData
# Output: IR_NK.RData (updated; read by figures/ligand_receptor_dotplot.R)

source("config.R")
# upstream_dir: folder with the inputs produced before this script (see README, "Upstream inputs")
upstream_dir <- Sys.getenv("SCRNA_UPSTREAM_DIR", unset = "upstream")

# ---- original lines 1914-2011 ----------------------------------------
#-------------------------------------2023/06/07-----------------------------------NK cells
load(file.path(data_dir, "IR_NK.RData"))
#clustering and annotation
DefaultAssay(IR_NK) <- "integrated"
ElbowPlot(IR_NK)
IR_NK <- FindNeighbors(IR_NK, dims = 1:20)
IR_NK <- FindClusters(object = IR_NK,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_NK <- RunUMAP(IR_NK, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_NK) <- "RNA"
Idents(IR_NK) = IR_NK@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_NK, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_NK, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
#remove low quality cells
IR_NK = subset(IR_NK, idents = c(0,1,2,3))

#clustering and annotation
DefaultAssay(IR_NK) <- "integrated"
ElbowPlot(IR_NK)
IR_NK <- FindNeighbors(IR_NK, dims = 1:20)
IR_NK <- FindClusters(object = IR_NK,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_NK <- RunUMAP(IR_NK, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_NK) <- "RNA"
Idents(IR_NK) = IR_NK@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_NK, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_NK) = IR_NK@meta.data$integrated_snn_res.0.6
DimPlot(subset(IR_NK, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_NK) = IR_NK@meta.data$integrated_snn_res.0.8
DimPlot(subset(IR_NK, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_NK) = IR_NK@meta.data$integrated_snn_res.1
DimPlot(subset(IR_NK, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_NK) = IR_NK@meta.data$integrated_snn_res.1.4
DimPlot(subset(IR_NK, downsample = 2000), label.size = 4,repel = T,label = T)

#Use 0.4 as resolution
Idents(object = IR_NK) = IR_NK@meta.data$integrated_snn_res.0.4
IR_NK[["Yan_annotate"]] <- Idents(object = IR_NK)
IR_NK$integrated_snn_res.1 <- NULL
IR_NK$integrated_snn_res.1.4 <- NULL
IR_NK$integrated_snn_res.0.4 <- NULL
IR_NK$integrated_snn_res.0.6 <- NULL
IR_NK$integrated_snn_res.0.8 <- NULL
IR_NK$seurat_clusters <- NULL

Idents(IR_NK) = IR_NK@meta.data$Yan_annotate
DimPlot(subset(IR_NK, downsample = 2000), label.size = 4,repel = T,label = T)

#Find markers
DefaultAssay(IR_NK) <- "RNA"
IR_NK <- NormalizeData(IR_NK)
IR_NK <- FindVariableFeatures(IR_NK, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_NK)
IR_NK <- ScaleData(IR_NK, features = all.genes)
all.markers <- FindAllMarkers(IR_NK, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "nkcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_NK, features = top5_markers$gene) + NoLegend()
DotPlot(IR_NK, features = c("Klrb1c","Cd3d","Cd4","Cd8b1")) + RotatedAxis()

#merge c0 and c4, remove c1(cd8T cell)
IR_NK = subset(IR_NK, idents = c(0,2,3,4))
DotPlot(IR_NK, features = c("Klrb1c","Cd3d","Cd4","Cd8b1")) + RotatedAxis()
DimPlot(IR_NK, label.size = 4,repel = T,label = T)
t = colnames(IR_NK)[which(IR_NK@meta.data$Yan_annotate %in% c(0,4))]
Idents(IR_NK, cells = t) = "0"
IR_NK@meta.data$Yan_annotate = Idents(IR_NK)

#Find markers
DefaultAssay(IR_NK) <- "RNA"
IR_NK <- NormalizeData(IR_NK)
IR_NK <- FindVariableFeatures(IR_NK, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_NK)
IR_NK <- ScaleData(IR_NK, features = all.genes)
all.markers <- FindAllMarkers(IR_NK, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "nkcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_NK, features = top5_markers$gene) + NoLegend()
DotPlot(IR_NK, features = c("Klrb1c","Cd3d","Cd4","Cd8b1")) + RotatedAxis()
DotPlot(IR_NK, features = c("Tbr1")) + RotatedAxis()
c3markers <- FindMarkers(IR_NK, ident.1 = 3, ident.2 = 2)
DotPlot(IR_NK, features = c("Tbr1","Itga2","Ifng","Gzma","Gzmb","Il7r","Itga1","Gata3","Rorc")) + RotatedAxis()
DotPlot(IR_NK, features = c("Tbr1","Gata3","Rorc")) + RotatedAxis()
## cluster 2 and 3 has very similar feature, merge
t = WhichCells(IR_NK, idents = c(2,3))
Idents(IR_NK, cells = t) = "Gzm+ NK"
t = WhichCells(IR_NK, idents = c(0))
Idents(IR_NK, Icells = t) = "ILC1"
IR_NK@meta.data$Yan_annotate = Idents(IR_NK)
table(IR_NK@meta.data$Yan_annotate,IR_NK@meta.data$orig.ident)
save(IR_NK, file = file.path(data_dir, "IR_NK.RData"))
table(IR_NK@meta.data$Yan_annotate,IR_NK@meta.data$orig.ident)


