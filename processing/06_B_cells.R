# 06. B cells: removal of dead cells and of cells bound to macrophages, clustering and annotation
#
# Processing record, as run interactively on the institutional HPC (R 4.1.1, Seurat 4.1.0) in 2023.
# Source: "20230512 radiation analysis.R", lines 2367-2466. The code is unchanged except that
#   (i) absolute HPC paths are replaced by file.path(data_dir, ...) / file.path(upstream_dir, ...),
#   (ii) setwd() and install.packages() calls are commented out,
#   (iii) console output and shell lines that had been pasted into the script are commented out.
# It documents how the objects used by the figure scripts were derived; it was not re-executed
# for this deposit. Exploratory plotting calls between the processing steps are kept so that the
# record of each decision (resolution, clusters removed) stays next to the step it justifies.
# Input : IR_B.RData
# Output: IR_B.RData (updated; read by figures/ligand_receptor_dotplot.R)

source("config.R")
# upstream_dir: folder with the inputs produced before this script (see README, "Upstream inputs")
upstream_dir <- Sys.getenv("SCRNA_UPSTREAM_DIR", unset = "upstream")

# ---- original lines 2367-2466 ----------------------------------------
#-----------------------B cells-----------------------------------------------------------
load(file.path(data_dir, "IR_B.RData"))
#CD8 T cells clustering and annotation
DefaultAssay(IR_B) <- "integrated"
ElbowPlot(IR_B)
IR_B <- FindNeighbors(IR_B, dims = 1:20)
IR_B <- FindClusters(object = IR_B,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_B <- RunUMAP(IR_B, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_B) <- "RNA"
Idents(IR_B) = IR_B@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_B, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_B, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
IR_B = subset(IR_B, idents = c(0,2,3,4,5))#remove dead cells
DefaultAssay(IR_B) <- "integrated"
IR_B <- FindNeighbors(IR_B, dims = 1:20)
IR_B <- FindClusters(object = IR_B,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
IR_B <- RunUMAP(IR_B, reduction = "pca", dims = 1:20, verbose = F)
DefaultAssay(IR_B) <- "RNA"
Idents(IR_B) = IR_B@meta.data$integrated_snn_res.0.4
DimPlot(subset(IR_B, downsample = 2000), label.size = 4,repel = T,label = T)
VlnPlot(IR_B, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
Idents(IR_B) = IR_B@meta.data$integrated_snn_res.0.6
DimPlot(subset(IR_B, downsample = 2000), label.size = 4,repel = T,label = T)

#Use 0.6 as resolution
Idents(object = IR_B) = IR_B@meta.data$integrated_snn_res.0.6
IR_B[["Yan_annotate"]] <- Idents(object = IR_B)
IR_B$integrated_snn_res.1 <- NULL
IR_B$integrated_snn_res.1.4 <- NULL
IR_B$integrated_snn_res.0.4 <- NULL
IR_B$integrated_snn_res.0.6 <- NULL
IR_B$integrated_snn_res.0.8 <- NULL
IR_B$seurat_clusters <- NULL

#Find markers
DefaultAssay(IR_B) <- "RNA"
IR_B <- NormalizeData(IR_B)
IR_B <- FindVariableFeatures(IR_B, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_B)
IR_B <- ScaleData(IR_B, features = all.genes)
all.markers <- FindAllMarkers(IR_B, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "Bcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_B, features = top5_markers$gene) + NoLegend()
FeaturePlot(IR_B, features = c("C1qb", "C1qc","Fcer1g","Rgs1","Ctsb"), label=T) 
a = WhichCells(IR_B, slot = 'counts', expression = Fcer1g == 0)
IR_B = subset(IR_B, cells = a)#remove cell binds to macrophages
a = WhichCells(IR_B, slot = 'counts', expression = C1qb == 0)
IR_B = subset(IR_B, cells = a)#remove cell binds to macrophages
FeaturePlot(IR_B, features = c("C1qb", "C1qc","Fcer1g","Rgs1","Ctsb"), label=T) 
DefaultAssay(IR_B) <- "RNA"
IR_B <- NormalizeData(IR_B)
IR_B <- FindVariableFeatures(IR_B, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_B)
IR_B <- ScaleData(IR_B, features = all.genes)
all.markers <- FindAllMarkers(IR_B, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "Bcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_B, features = top5_markers$gene) + NoLegend()


#Annotate B cells
#Group 4&6
#Group 0&1&2&3&8
a = WhichCells(IR_B, idents = c(4, 6))
Idents(IR_B, cells = a) = "4"
a = WhichCells(IR_B, idents = c(0,1,2,3,8))
Idents(IR_B, cells = a) = "0"
table(Idents(IR_B))

# [console output] 0    4    5    7
# [console output] 8669 1428  858  285
DefaultAssay(IR_B) <- "RNA"
IR_B <- NormalizeData(IR_B)
IR_B <- FindVariableFeatures(IR_B, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(IR_B)
IR_B <- ScaleData(IR_B, features = all.genes)
all.markers <- FindAllMarkers(IR_B, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library("writexl")
write_xlsx(all.markers, file.path(data_dir, "Bcell_markers.xlsx"))
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
DoHeatmap(IR_B, features = top5_markers$gene) + NoLegend()
#Premature B cell marker
DotPlot(IR_B, features = c("Bub1","Ckap21","Dlgap5","Cep55","Neil3","Cdca5","Sgo1","Ska1","Ncapg","Pax5","Nut2")) + RotatedAxis()
DimPlot(subset(IR_B, downsample = 2000), label.size = 4,repel = T,label = T)
#Mature B cell markers
DotPlot(IR_B, features = c("Ms4a1","Ly6a","Faim3","Fcer2a","Bank1","Bcl2","Ltb","Fcrl1","Ctsh","H2-DM-a")) + RotatedAxis()

DotPlot(IR_B, features = c("Ighg1")) + RotatedAxis()
FeaturePlot(IR_B, features = c("Ms4a1","Ly6a","Faim3","Fcer2a","Bank1","Bcl2","Ltb","Fcrl1","Ctsh","H2-DM-a","Ighg1"), label=T) 
IR_B[["Yan_annotate"]] <- Idents(object = IR_B)
save(IR_B, file = file.path(data_dir, "IR_B.RData"))
table(IR_B@meta.data$Yan_annotate, IR_B@meta.data$orig.ident)
