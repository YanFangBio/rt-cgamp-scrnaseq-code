# 04. Myeloid cells: removal of contaminants, DC subset, annotation of cDC1 / cDC2 / CCR7+ (mature DC)
#
# Processing record, as run interactively on the institutional HPC (R 4.1.1, Seurat 4.1.0) in 2023.
# Source: "20230512 radiation analysis.R", lines 2012-2277. The code is unchanged except that
#   (i) absolute HPC paths are replaced by file.path(data_dir, ...) / file.path(upstream_dir, ...),
#   (ii) setwd() and install.packages() calls are commented out,
#   (iii) console output and shell lines that had been pasted into the script are commented out.
# It documents how the objects used by the figure scripts were derived; it was not re-executed
# for this deposit. Exploratory plotting calls between the processing steps are kept so that the
# record of each decision (resolution, clusters removed) stays next to the step it justifies.
# Input : IR_myeloid.RData
# Output: IR_myeloid.RData (updated), IR_myeloid_DC.RData (input of figures/DC_subsets.R, figures/DC_feature_plots_Mafb_S100a8.R
#         and figures/ligand_receptor_dotplot.R)

source("config.R")
# upstream_dir: folder with the inputs produced before this script (see README, "Upstream inputs")
upstream_dir <- Sys.getenv("SCRNA_UPSTREAM_DIR", unset = "upstream")

# ---- original lines 2012-2277 ----------------------------------------
#DC
load(file.path(data_dir, "IR_myeloid.RData"))
DefaultAssay(IR_myeloid) <- "integrated"
ElbowPlot(IR_myeloid)
IR_myeloid <- FindNeighbors(IR_myeloid, dims = 1:20)
IR_myeloid <- FindClusters(object = IR_myeloid,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_myeloid <- RunUMAP(IR_myeloid, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_myeloid) <- "RNA"
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_myeloid, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
#remove low quality cells
IR_myeloid = subset(IR_myeloid, idents = c(0,1,2,3,4,6,7,8,9,10))
DefaultAssay(IR_myeloid) <- "integrated"
ElbowPlot(IR_myeloid)
IR_myeloid <- FindNeighbors(IR_myeloid, dims = 1:20)
IR_myeloid <- FindClusters(object = IR_myeloid,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_myeloid <- RunUMAP(IR_myeloid, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_myeloid) <- "RNA"
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_myeloid, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
IR_myeloid@meta.data$Yan_annotate = Idents(IR_myeloid)
IR_myeloid@meta.data$integrated_snn_res.1 = NULL
IR_myeloid@meta.data$integrated_snn_res.0.4 = NULL
IR_myeloid@meta.data$integrated_snn_res.0.6 = NULL
IR_myeloid@meta.data$integrated_snn_res.0.8 = NULL
IR_myeloid@meta.data$integrated_snn_res.1.4 = NULL
IR_myeloid@meta.data$seurat_clusters = NULL
save(IR_myeloid, file = file.path(data_dir, "IR_myeloid.RData"))

#Find markers
DefaultAssay(IR_myeloid) <- "RNA"
IR_myeloid <- NormalizeData(IR_myeloid)
IR_myeloid <- FindVariableFeatures(IR_myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_myeloid)
IR_myeloid <- ScaleData(IR_myeloid, features = all.genes)
all.markers <- FindAllMarkers(IR_myeloid, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "IR_myeloid_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_myeloid, features = top5_markers$gene) + NoLegend()
DotPlot(IR_myeloid, features = c("Klrb1c","Cd3d","Cd4","Cd8b1")) + RotatedAxis()
#5 and 10 are T cells,
table(IR_myeloid@meta.data$Yan_annotate[which(IR_myeloid@meta.data$tcr_cluster != 'NA')])
IR_myeloid = subset(IR_myeloid, idents = c(0,1,2,3,4,6,7,8,9))

#reclustering
DefaultAssay(IR_myeloid) <- "integrated"
ElbowPlot(IR_myeloid)
IR_myeloid <- FindNeighbors(IR_myeloid, dims = 1:20)
IR_myeloid <- FindClusters(object = IR_myeloid,
                               resolution = c(0.4))
IR_myeloid <- RunUMAP(IR_myeloid, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_myeloid) <- "RNA"
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_myeloid, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
IR_myeloid@meta.data$Yan_annotate = Idents(IR_myeloid)
IR_myeloid@meta.data$integrated_snn_res.0.4 = NULL
IR_myeloid@meta.data$seurat_clusters = NULL
save(IR_myeloid, file = file.path(data_dir, "IR_myeloid.RData"))
#Find markers
DefaultAssay(IR_myeloid) <- "RNA"
IR_myeloid <- NormalizeData(IR_myeloid)
IR_myeloid <- FindVariableFeatures(IR_myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_myeloid)
IR_myeloid <- ScaleData(IR_myeloid, features = all.genes)
all.markers <- FindAllMarkers(IR_myeloid, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "IR_myeloid_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_myeloid, features = top5_markers$gene) + NoLegend()

# cluster 9 is B cell, remove
IR_myeloid = subset(IR_myeloid, idents = c(0,1,2,3,4,5,6,7,8))
DefaultAssay(IR_myeloid) <- "RNA"
IR_myeloid <- NormalizeData(IR_myeloid)
IR_myeloid <- FindVariableFeatures(IR_myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_myeloid)
IR_myeloid <- ScaleData(IR_myeloid, features = all.genes)
all.markers <- FindAllMarkers(IR_myeloid, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "IR_myeloid_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_myeloid, features = top5_markers$gene) + NoLegend()
save(IR_myeloid, file = file.path(data_dir, "IR_myeloid.RData"))

DotPlot(IR_myeloid, features = c("Itgax","Adgre1","Cd8a1","Itgam", "Xcr1", "Batf3")) + RotatedAxis()
#pDC
DotPlot(IR_myeloid, features = c("Irf4","Irf7","Irf8","Ccr7", "Cd209a", "Lilrb4a","Nrp1","Ifnb1","Siglech")) + RotatedAxis()
#cDC1
DotPlot(IR_myeloid, features = c("Batf3","Id2","Irf8","Zbtb46", "Btla", "Cadm1","Cd8a","Itgae","Itgax","Ly75","Thbd","Xcr1","Clec9a")) + RotatedAxis()
#cDC2
DotPlot(IR_myeloid, features = c("Id2","Irf4","Klf4","Zbtb46", "Cd14", "Cd163","Clec10a","Notch2","Itgam","Sirpa","Cx3cr1","Cd1","Cd2")) + RotatedAxis()
DotPlot(IR_myeloid, features = c("Itgam","Tbx21","Runx3","Srebf2","Cd86","Rorc","Cebpa", "Clec12a","Clec10a","Tnf","Il6")) + RotatedAxis()
save(IR_myeloid, file = file.path(data_dir, "IR_myeloid.RData"))

#By checking the F4/80 and the Batf3, I think the granularity is not enough. So reclustering based on filtered cells.
load(file.path(data_dir, "IR_myeloid.RData"))
#reclustering
DefaultAssay(IR_myeloid) <- "integrated"
ElbowPlot(IR_myeloid)
IR_myeloid <- FindNeighbors(IR_myeloid, dims = 1:20)
IR_myeloid <- FindClusters(object = IR_myeloid,
                               resolution = c(0.4, 0.8, 1, 1.2, 1.4, 1.6, 1.8))
IR_myeloid <- RunUMAP(IR_myeloid, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_myeloid) <- "RNA"
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.0.8
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.1
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.1.2
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.1.4
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.1.6
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid) = IR_myeloid@meta.data$integrated_snn_res.1.8
DimPlot(subset(IR_myeloid, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_myeloid, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
FeaturePlot(IR_myeloid, features = c("Itgax","Adgre1","Cd8a","Itgam", "Xcr1", "Batf3"), label=T) 

IR_myeloid@meta.data$Yan_annotate = Idents(IR_myeloid)
IR_myeloid@meta.data$integrated_snn_res.0.4 = NULL
IR_myeloid@meta.data$integrated_snn_res.0.8 = NULL
IR_myeloid@meta.data$integrated_snn_res.1 = NULL
IR_myeloid@meta.data$integrated_snn_res.1.2 = NULL
IR_myeloid@meta.data$integrated_snn_res.1.4 = NULL
IR_myeloid@meta.data$integrated_snn_res.1.6 = NULL
IR_myeloid@meta.data$integrated_snn_res.1.8 = NULL
IR_myeloid@meta.data$seurat_clusters = NULL
DotPlot(IR_myeloid, features = c("Itgax","Adgre1","Cd8a","Itgam", "Xcr1", "Batf3","H2-Ab1")) + RotatedAxis()
#pDC
DotPlot(IR_myeloid, features = c("Irf4","Irf7","Irf8","Ccr7", "Cd209a", "Lilrb4a","Nrp1","Ifnb1","Siglech")) + RotatedAxis()
#cDC1
DotPlot(IR_myeloid, features = c("Batf3","Id2","Irf8","Zbtb46", "Btla", "Cadm1","Cd8a","Itgae","Itgax","Ly75","Thbd","Xcr1","Clec9a")) + RotatedAxis()
#cDC2
DotPlot(IR_myeloid, features = c("Id2","Irf4","Klf4","Zbtb46", "Cd14", "Cd163","Clec10a","Notch2","Itgam","Sirpa","Cx3cr1","Cd1","Cd2")) + RotatedAxis()
DotPlot(IR_myeloid, features = c("Itgam","Tbx21","Runx3","Srebf2","Cd86","Rorc","Cebpa", "Clec12a","Clec10a","Tnf","Il6")) + RotatedAxis()

save(IR_myeloid, file = file.path(data_dir, "IR_myeloid.RData"))

#Find markers
DefaultAssay(IR_myeloid) <- "RNA"
IR_myeloid <- NormalizeData(IR_myeloid)
IR_myeloid <- FindVariableFeatures(IR_myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_myeloid)
IR_myeloid <- ScaleData(IR_myeloid, features = all.genes)
all.markers <- FindAllMarkers(IR_myeloid, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "IR_myeloid_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_myeloid, features = top5_markers$gene) + NoLegend()
table(IR_myeloid@meta.data$Yan_annotate,IR_myeloid@meta.data$orig.ident)
RidgePlot(object = IR_myeloid, features = c("Cd86"))
RidgePlot(object = IR_myeloid, features = c("Cd80"))
RidgePlot(object = IR_myeloid, features = c("Cd274"))
Idents(IR_myeloid) = IR_myeloid@meta.data$orig.ident

#Define the DC/macrophage/monocytes
Idents(IR_myeloid) = IR_myeloid@meta.data$Yan_annotate
DotPlot(IR_myeloid, features = c("Ptprc","Adgre1","Itgam","Itgax","Ly6g", "Ly6c2","Zbtb46","H2-Aa", "H2-Eb1","H2-Ab1")) + RotatedAxis()
FeaturePlot(IR_myeloid, features = c("Ly6g", "Ly6c2"), label=T) 

## Subset DC
IR_myeloid_DC = subset(IR_myeloid, idents = c(11,13,18,20,21,22)) 
DotPlot(IR_myeloid_DC, features = c("Ptprc","Adgre1","Itgam","Itgax","Ly6g", "Ly6c2","Zbtb46","H2-Aa", "H2-Eb1","H2-Ab1")) + RotatedAxis()
#reclustering
DefaultAssay(IR_myeloid_DC) <- "integrated"
ElbowPlot(IR_myeloid_DC)
IR_myeloid_DC <- FindNeighbors(IR_myeloid_DC, dims = 1:20)
IR_myeloid_DC <- FindClusters(object = IR_myeloid_DC,
                               resolution = c(0.4, 0.8, 1, 1.2, 1.4, 1.6, 1.8))
IR_myeloid_DC <- RunUMAP(IR_myeloid_DC, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_myeloid_DC) <- "RNA"
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_myeloid_DC, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$integrated_snn_res.0.8
DimPlot(subset(IR_myeloid_DC, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$integrated_snn_res.1
DimPlot(subset(IR_myeloid_DC, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_myeloid_DC, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
#remove 4 and 8 which have high mitochondria
IR_myeloid_DC = subset(IR_myeloid_DC, idents = c(0,1,2,3,5,6,7,9)) 
#reclustering
DefaultAssay(IR_myeloid_DC) <- "integrated"
ElbowPlot(IR_myeloid_DC)
IR_myeloid_DC <- FindNeighbors(IR_myeloid_DC, dims = 1:20)
IR_myeloid_DC <- FindClusters(object = IR_myeloid_DC,
                               resolution = c(0.4, 0.8, 1, 1.2, 1.4, 1.6, 1.8))
IR_myeloid_DC <- RunUMAP(IR_myeloid_DC, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_myeloid_DC) <- "RNA"
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_myeloid_DC, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$integrated_snn_res.0.8
DimPlot(subset(IR_myeloid_DC, downsample = 2000), label.size = 4,repel = T,label = T)
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$integrated_snn_res.1
DimPlot(subset(IR_myeloid_DC, downsample = 2000), label.size = 4,repel = T,label = T)

DotPlot(IR_myeloid_DC, features = c("Ptprc","Itgam","Itgax","Zbtb46","H2-Aa", "H2-Eb1","H2-Ab1","H2-K1","H2-D1")) + RotatedAxis()

IR_myeloid_DC@meta.data$Yan_annotate = Idents(IR_myeloid_DC)
IR_myeloid_DC@meta.data$integrated_snn_res.0.4 = NULL
IR_myeloid_DC@meta.data$integrated_snn_res.0.8 = NULL
IR_myeloid_DC@meta.data$integrated_snn_res.1 = NULL
IR_myeloid_DC@meta.data$integrated_snn_res.1.2 = NULL
IR_myeloid_DC@meta.data$integrated_snn_res.1.4 = NULL
IR_myeloid_DC@meta.data$integrated_snn_res.1.6 = NULL
IR_myeloid_DC@meta.data$integrated_snn_res.1.8 = NULL
IR_myeloid_DC@meta.data$seurat_clusters = NULL

DotPlot(IR_myeloid_DC, features = c("Itgax","Adgre1","Cd8a","Itgam", "Xcr1", "Batf3","H2-Ab1")) + RotatedAxis()
#pDC
DotPlot(IR_myeloid_DC, features = c("Irf4","Irf7","Irf8","Ccr7", "Cd209a", "Lilrb4a","Nrp1","Ifnb1","Siglech")) + RotatedAxis()
#cDC1
DotPlot(IR_myeloid_DC, features = c("Batf3","Id2","Irf8","Zbtb46", "Btla", "Cadm1","Cd8a","Itgae","Itgax","Ly75","Thbd","Xcr1","Clec9a")) + RotatedAxis()
#cDC2
DotPlot(IR_myeloid_DC, features = c("Irf4","Klf4", "Cd14", "Cd163","Clec10a","Notch2","Itgam","Sirpa","Cx3cr1","Cd1","Cd2")) + RotatedAxis()
DotPlot(IR_myeloid_DC, features = c("Itgam","Tbx21","Runx3","Srebf2","Cd86","Rorc","Cebpa", "Clec12a","Clec10a","Tnf","Il6")) + RotatedAxis()
RidgePlot(object = IR_myeloid_DC, features = c("Cd86"))
RidgePlot(object = IR_myeloid_DC, features = c("Cd80"))
RidgePlot(object = IR_myeloid_DC, features = c("Cd274"))
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$orig.ident
Idents(IR_myeloid_DC) = IR_myeloid_DC@meta.data$Yan_annotate
#Find markers
DefaultAssay(IR_myeloid_DC) <- "RNA"
IR_myeloid_DC <- NormalizeData(IR_myeloid_DC)
IR_myeloid_DC <- FindVariableFeatures(IR_myeloid_DC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_myeloid_DC)
IR_myeloid_DC <- ScaleData(IR_myeloid_DC, features = all.genes)
all.markers <- FindAllMarkers(IR_myeloid_DC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "IR_myeloid_DC_markers.xlsx"))
top10_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
top10_markers
DoHeatmap(IR_myeloid_DC, features = top10_markers$gene) + NoLegend()
DotPlot(IR_myeloid_DC, features = c("Ccr2","Ccr7")) + RotatedAxis()
table(IR_myeloid_DC@meta.data$Yan_annotate,IR_myeloid_DC@meta.data$orig.ident)
t = WhichCells(IR_myeloid_DC, idents = c(2,4,5))
Idents(IR_myeloid_DC, cells = t) = "CCR7+"
t = WhichCells(IR_myeloid_DC, idents = c(0,1))
Idents(IR_myeloid_DC, cells = t) = "cDC2"
t = WhichCells(IR_myeloid_DC, idents = c(3))
Idents(IR_myeloid_DC, cells = t) = "cDC1"
IR_myeloid_DC@meta.data$Yan_annotate =Idents(IR_myeloid_DC)
FeaturePlot(IR_myeloid_DC, features = c("Ccr7", "Il12b"), label=T) 







DotPlot(IR_myeloid_DC, features = c("Cd209a","Il4ra")) + RotatedAxis()
DotPlot(IR_myeloid_DC, features = c("Fas","Aldh1a2")) + RotatedAxis()
DotPlot(IR_myeloid_DC, features = c("Sirpa","H2-DMb2","Itgam")) + RotatedAxis()
save(IR_myeloid_DC, file = file.path(data_dir, "IR_myeloid_DC.RData"))
