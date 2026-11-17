# 01. Doublet removal, clustering of all cells, TCR annotation and split into lineages
#
# Processing record, as run interactively on the institutional HPC (R 4.1.1, Seurat 4.1.0) in 2023.
# Source: "20230512 radiation analysis.R", lines 1-283. The code is unchanged except that
#   (i) absolute HPC paths are replaced by file.path(data_dir, ...) / file.path(upstream_dir, ...),
#   (ii) setwd() and install.packages() calls are commented out,
#   (iii) console output and shell lines that had been pasted into the script are commented out.
# It documents how the objects used by the figure scripts were derived; it was not re-executed
# for this deposit. Exploratory plotting calls between the processing steps are kept so that the
# record of each decision (resolution, clusters removed) stays next to the step it justifies.
# Input : TIL_seurat_bf_cluster.RData (QC-filtered, DoubletFinder-labelled, integrated object; produced upstream)
#         TCR_cdr3b_RotationEncodingBL62.txt (TCR beta CDR3 clusters per cell barcode; produced upstream)
# Output: TIL_seurat_IR.RData, IR_CD8T.RData, IR_CD4T.RData, IR_NK.RData, IR_B.RData, IR_myeloid.RData
# Note  : clusters 18, 21 and 27 were split by hand with Seurat::CellSelector() (interactive). The cells
#         chosen are recorded in the lineage objects, so this step cannot be replayed from code alone.

source("config.R")
# upstream_dir: folder with the inputs produced before this script (see README, "Upstream inputs")
upstream_dir <- Sys.getenv("SCRNA_UPSTREAM_DIR", unset = "upstream")

# ---- original lines 1-283 ----------------------------------------
# module load R/4.1.1-gccmkl 
# R

#install.packages('Seurat')
#load all the libraries that will be needed for the tutorial.
library(Seurat)
packageVersion("Seurat")#4.1.0
library(ggplot2)
#library(SingleR)
library(dplyr)
#library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

load(file.path(upstream_dir, "TIL_seurat_bf_cluster.RData"))
DefaultAssay(TIL_seurat) <- "integrated"
DimPlot(subset(TIL_seurat, downsample = 2000), reduction = "pca")
ElbowPlot(TIL_seurat,ndims=30)
TIL_seurat <- FindNeighbors(TIL_seurat, dims = 1:30)
TIL_seurat <- FindClusters(TIL_seurat, resolution = 0.8)
TIL_seurat <- RunUMAP(TIL_seurat, reduction = "pca", dims = 1:30, verbose = F)
table(TIL_seurat@meta.data$seurat_clusters)
## View doublets clusters
Idents(TIL_seurat) = TIL_seurat@meta.data$seurat_clusters
DimPlot(subset(TIL_seurat, idents = c("8","28","31","35")), label.size = 4,repel = T,label = T)
t = WhichCells(TIL_seurat, idents = c("8","28","31","35"))
Idents(TIL_seurat) = TIL_seurat[['doublet.status']]
DimPlot(subset(TIL_seurat, downsample = 2000), label.size =4, repel=T, label=F)
DefaultAssay(TIL_seurat) <- "RNA"
DotPlot(TIL_seurat, features = c("Clec10a","Cd9","Adgre1","Folr2","Itgax","Cd3d","Cd4","Cd8b1","Klrb1c")) + RotatedAxis()
FeaturePlot(subset(TIL_seurat, downsample = 2000), features = c("Clec10a","Cd9","Adgre1","Folr2","Itgax"))
FeaturePlot(subset(TIL_seurat, downsample = 2000), features = c("Cd3d","Cd4","Cd8b1","Klrb1c"))

#Split sample into doublets or singlets
## Singlets
## remove doublet clusters
Idents(TIL_seurat) = TIL_seurat@meta.data$seurat_clusters
t = WhichCells(TIL_seurat, idents = c("0","1","2","3","4", "5", "6","7","9","10", "11", "12","13","14","15","16", "17","18","19","20","21","22", "23", "24","25","26","27", "29", "30","32","33","34"))
TIL_seurat_single = subset(TIL_seurat, cells = t)
## remove labeled doublets
Idents(TIL_seurat_single) = TIL_seurat_single@meta.data$doublet.status
t = WhichCells(TIL_seurat_single, idents = c("Singlet"))
TIL_seurat_single = subset(TIL_seurat_single, cells = t) ##87460 cells singlets
##doublets
Idents(TIL_seurat) = TIL_seurat@meta.data$doublet.status
t1 = WhichCells(TIL_seurat, idents = c("Doublet"))
Idents(TIL_seurat) = TIL_seurat@meta.data$seurat_clusters
t2 = WhichCells(TIL_seurat, idents = c("8","28","31","35"))
t = c(t1[which(!(t1 %in% t2))],t2)
TIL_seurat_double = subset(TIL_seurat, cells = t) ##7825 cells doublets
FeaturePlot(TIL_seurat_double, features = c("Clec10a","Cd9","Adgre1","Folr2","Itgax"))
FeaturePlot(TIL_seurat_double, features = c("Cd3d","Cd4","Cd8b1","Klrb1c"))
FeaturePlot(TIL_seurat_double, features = c("Clec10a","Adgre1","Cd9","Cd4","Cd8b1","Klrb1c"))

##Only take control and IR group for further analysis
Idents(TIL_seurat) = TIL_seurat@meta.data$orig.ident
t = WhichCells(TIL_seurat, idents = c("local_ctrl_T", "distal_ctrl_T", "local_IR_T", "distal_IR_T", "local_ctrl_LN", "distal_ctrl_LN", "local_IR_LN", "distal_IR_LN"))
TIL_seurat_IR = subset(TIL_seurat, cells = t)
# setwd(file.path(upstream_dir, "RadiationAnalysis"))   # working directory on the HPC; files are addressed through data_dir / upstream_dir below
save(TIL_seurat_IR, file = file.path(data_dir, "TIL_seurat_IR.RData"))
#------------------05122023--------------------
load(file.path(data_dir, "TIL_seurat_IR.RData"))
#remove labeled doublets and clusters:
Idents(TIL_seurat_IR) = TIL_seurat_IR@meta.data$seurat_clusters
t = WhichCells(TIL_seurat_IR, idents = c("0","1","2","3","4", "5", "6","7","9","10", "11", "12","13","14","15","16", "17","18","19","20","21","22", "23", "24","25","26","27", "29", "30","32","33","34"))
TIL_seurat_IR = subset(TIL_seurat_IR, cells = t)
## remove labeled doublets
Idents(TIL_seurat_IR) = TIL_seurat_IR@meta.data$doublet.status
t = WhichCells(TIL_seurat_IR, idents = c("Singlet"))
TIL_seurat_IR = subset(TIL_seurat_IR, cells = t) ##60661 cells singlets
# setwd(file.path(upstream_dir, "RadiationAnalysis"))   # working directory on the HPC; files are addressed through data_dir / upstream_dir below
save(TIL_seurat_IR, file = file.path(data_dir, "TIL_seurat_IR.RData"))
#
# reclustering the filtered cells with resolution 0.4, 0.6, 0.8, 1.0, 1.4. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.
DefaultAssay(TIL_seurat_IR) <- "integrated"
TIL_seurat_IR <- FindNeighbors(TIL_seurat_IR, dims = 1:30)
TIL_seurat_IR <- FindClusters(object = TIL_seurat_IR,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
TIL_seurat_IR <- RunUMAP(TIL_seurat_IR, reduction = "pca", dims = 1:30, verbose = F)
DefaultAssay(TIL_seurat_IR) <- "RNA"
Idents(TIL_seurat_IR) = TIL_seurat_IR@meta.data$integrated_snn_res.0.4
DimPlot(subset(TIL_seurat_IR, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(TIL_seurat_IR) = TIL_seurat_IR@meta.data$integrated_snn_res.0.6
DimPlot(subset(TIL_seurat_IR, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(TIL_seurat_IR) = TIL_seurat_IR@meta.data$integrated_snn_res.0.8
DimPlot(subset(TIL_seurat_IR, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(TIL_seurat_IR) = TIL_seurat_IR@meta.data$integrated_snn_res.1
DimPlot(subset(TIL_seurat_IR, downsample = 2000), label.size = 4,repel = T,label = T)

Idents(TIL_seurat_IR) = TIL_seurat_IR@meta.data$integrated_snn_res.1.4
DimPlot(subset(TIL_seurat_IR, downsample = 2000), label.size = 4,repel = T,label = T)

FeaturePlot(TIL_seurat_IR, features = c("Cd3d","Cd4","Cd8b1","Klrb1c", "Cd19"))

Idents(TIL_seurat_IR) = TIL_seurat_IR[['tissue']]
t = WhichCells(TIL_seurat_IR, idents=c('LN isolated cells'))
Idents(TIL_seurat_IR) = TIL_seurat_IR[['integrated_snn_res.1.4']]
DimPlot(subset(TIL_seurat_IR, cells = t), label.size = 4,repel = T,label = T)
Idents(TIL_seurat_IR) = TIL_seurat_IR[['tissue']]
t = WhichCells(TIL_seurat_IR, idents=c('tumor isolated cells'))
Idents(TIL_seurat_IR) = TIL_seurat_IR[['integrated_snn_res.1.4']]
DimPlot(subset(TIL_seurat_IR, cells = t), label.size = 4,repel = T,label = T)

#Use resolution 1.0 ans double check the clustering
DimPlot(subset(TIL_seurat_IR, idents = c("0","1","2","3","4", "5")), label.size = 4,repel = T,label = T)
DimPlot(subset(TIL_seurat_IR, idents = c("6","7","8","9","10", "11")), label.size = 4,repel = T,label = T)
DimPlot(subset(TIL_seurat_IR, idents = c("12","13","14","15","16", "17")), label.size = 4,repel = T,label = T)
DimPlot(subset(TIL_seurat_IR, idents = c("18","19","20","21","22", "23")), label.size = 4,repel = T,label = T)
DimPlot(subset(TIL_seurat_IR, idents = c("24","25","26","27","28", "29")), label.size = 4,repel = T,label = T)
DimPlot(subset(TIL_seurat_IR, idents = c("30","31","32")), label.size = 4,repel = T,label = T)
# > TIL_seurat_IR$DF.classifications_0.25_0.22_866 <- NULL
# > TIL_seurat_IR$DF.classifications_0.25_0.28_450 <- NULL
# > TIL_seurat_IR$DF.classifications_0.25_0.3_117 <- NULL
# > TIL_seurat_IR$pANN_0.25_0.27_294 <- NULL
# > TIL_seurat_IR$pANN_0.25_0.3_117 <- NULL
# > TIL_seurat_IR$pANN_0.25_0.28_450 <- NULL
# > TIL_seurat_IR$pANN_0.25_0.22_866 <- NULL
# .... to clean the data

####--------note: no18, 21,27 don't clustered well

DefaultAssay(TIL_seurat_IR) <- "RNA"
Idents(TIL_seurat_IR) = TIL_seurat_IR[['integrated_snn_res.1']]
DotPlot(TIL_seurat_IR, features = c("Ptprc","Itgax","Itgam","Cd3d","Cd4","Cd8b1","Klrb1c", "Cd19")) + RotatedAxis()
DotPlot(TIL_seurat_IR, features = c("Ptprc","Itgax","Itgam","Cd3e","Cd3d","Cd4","Cd8b1","Klrb1c", "Nkg7","Cd19")) + RotatedAxis()

#8.4 Differential expression and marker selection
DefaultAssay(TIL_seurat_IR) <- "RNA"
TIL_seurat_IR <- NormalizeData(TIL_seurat_IR)
TIL_seurat_IR <- FindVariableFeatures(TIL_seurat_IR, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(TIL_seurat_IR)
TIL_seurat_IR <- ScaleData(TIL_seurat_IR, features = all.genes)
all.markers <- FindAllMarkers(TIL_seurat_IR, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

DoHeatmap(TIL_seurat_IR, features = top3_markers$gene) + NoLegend()
# setwd(file.path(upstream_dir, "RadiationAnalysis"))   # working directory on the HPC; files are addressed through data_dir / upstream_dir below
save(TIL_seurat_IR, file = file.path(data_dir, "TIL_seurat_IR.RData"))
#--------------------------------20230515-----------------------------------------------------------------------------------------
load(file.path(data_dir, "TIL_seurat_IR.RData"))

#integrate TCR information to define T cells (TCR clustered by CDR3b)
# Clonotype-centric info.
clono = read.delim(file.path(upstream_dir, "TCR_cdr3b_RotationEncodingBL62.txt"), skip = 2, header=F)
colnames (clono) = c("cdr3", "tcr_cluster", "v_gene","d_gene","j_gene","tissue","cell_barcode")
clono <- clono[!duplicated(clono$cell_barcode), ]#14740cells
clono <- clono[,c("cell_barcode", "tcr_cluster","cdr3")]
rownames(clono) <- clono[,1]
clono[,1] <- NULL



#change the cell name in TIL_seurat_IR

original_cell_name = colnames(TIL_seurat_IR)
t = grep("1_1", original_cell_name)
t2 = grep("1_10", original_cell_name)
t3 = grep("1_11", original_cell_name)
t4 = grep("1_12", original_cell_name)
cellname <- gsub("1_1", "1", original_cell_name[which(!(t %in% c(t2, t3, t4)))])
newname = cellname

t = grep("1_2", original_cell_name)
cellname <- gsub("1_2", "2", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_3", original_cell_name)
cellname <- gsub("1_3", "3", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_4", original_cell_name)
cellname <- gsub("1_4", "4", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_5", original_cell_name)
cellname <- gsub("1_5", "5", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_6", original_cell_name)
cellname <- gsub("1_6", "6", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_7", original_cell_name)
cellname <- gsub("1_7", "7", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_8", original_cell_name)
cellname <- gsub("1_8", "8", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_9", original_cell_name)
cellname <- gsub("1_9", "9", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_10", original_cell_name)
cellname <- gsub("1_10", "10", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_11", original_cell_name)
cellname <- gsub("1_11", "11", original_cell_name[t])
newname = c(newname,cellname)

t = grep("1_12", original_cell_name)
cellname <- gsub("1_12", "12", original_cell_name[t])
newname = c(newname,cellname)

TIL_seurat_IR = RenameCells(TIL_seurat_IR,new.names=newname)


TIL_seurat_IR <- AddMetaData(object=TIL_seurat_IR, metadata=clono)
a = TIL_seurat_IR@meta.data[which(TIL_seurat_IR$tcr_cluster != "<NA>"),]
t_cell_name = rownames(a)
DimPlot(TIL_seurat_IR, reduction = "umap", label = T, cells.highlight = t_cell_name)


ob.list <- SplitObject(TIL_seurat_IR, split.by = "orig.ident")
plot.list <- lapply(X = ob.list, FUN = function(x) {
    DimPlot(x, reduction = "umap", label = T, label.size = 4)
})

#Split 18,21,27
Idents(TIL_seurat_IR) = TIL_seurat_IR@meta

#18
plot <- DimPlot(subset(TIL_seurat_IR,idents = 18), reduction = "umap")
select.cells <- CellSelector(plot = plot)
cluster_group1 = select.cells
clustercell = WhichCells(TIL_seurat_IR, idents = c(18))
cluster_group2 = clustercell[which(!(clustercell%in%cluster_group1))]
TIL_seurat_IR[["Yan_annotate"]] <- Idents(object = TIL_seurat_IR)
Idents(TIL_seurat_IR, cells = cluster_group1) = "cluster18_cd8"
Idents(TIL_seurat_IR, cells = cluster_group2) = "cluster18_cd4"
TIL_seurat_IR[["Yan_annotate"]] <- Idents(object = TIL_seurat_IR)

#21
plot <- DimPlot(subset(TIL_seurat_IR,idents = 21), reduction = "umap")
select.cells <- CellSelector(plot = plot)
cluster_group1 = select.cells
clustercell = WhichCells(TIL_seurat_IR, idents = c(21))
cluster_group2 = clustercell[which(!(clustercell%in%cluster_group1))]
Idents(TIL_seurat_IR, cells = cluster_group1) = "cluster21_unknow"
Idents(TIL_seurat_IR, cells = cluster_group2) = "cluster21_B"
TIL_seurat_IR[["Yan_annotate"]] <- Idents(object = TIL_seurat_IR)

#27
plot <- DimPlot(subset(TIL_seurat_IR,idents = 27), reduction = "umap")
select.cells <- CellSelector(plot = plot)
cluster_group1 = select.cells
clustercell = WhichCells(TIL_seurat_IR, idents = c(27))
cluster_group2 = clustercell[which(!(clustercell%in%cluster_group1))]
Idents(TIL_seurat_IR, cells = cluster_group1) = "cluster27_CD8"
Idents(TIL_seurat_IR, cells = cluster_group2) = "cluster27_CD4"
TIL_seurat_IR[["Yan_annotate"]] <- Idents(object = TIL_seurat_IR)


#Split cells into CD8, CD4, NK,  B, myeloid and unknown.
##CD8 cells 1,5,7,11,16,cluster18_cd8,cluster27_CD8
##CD4 cells 2, 4, 9, 10, 14, 19, 24, cluster18_cd4, cluster27_CD4
##NK cells 8,23
##B cell 0, 3, 12, 13, 29, cluster21_B
##Myeloid 6, 15, 17,20, 22, 25, 26, cluster21_unknow
IR_CD8T <- subset(TIL_seurat_IR, idents = c(1,5,7,11,16,'cluster18_cd8','cluster27_CD8'))
IR_CD4T <- subset(TIL_seurat_IR, idents = c(2, 4, 9, 10, 14, 19, 24, 'cluster18_cd4', 'cluster27_CD4'))
IR_NK <- subset(TIL_seurat_IR, idents = c(8,23))
IR_B <- subset(TIL_seurat_IR, idents = c(0, 3, 12, 13, 29, 'cluster21_B'))
IR_myeloid <- subset(TIL_seurat_IR, idents = c(6, 15, 17,20, 22, 25, 26, 'cluster21_unknow'))

DimPlot(IR_CD8T, reduction = "umap", label = T, label.size = 4)
DimPlot(IR_CD4T, reduction = "umap", label = T, label.size = 4)
DimPlot(IR_NK, reduction = "umap", label = T, label.size = 4)
DimPlot(IR_B, reduction = "umap", label = T, label.size = 4)
DimPlot(IR_myeloid, reduction = "umap", label = T, label.size = 4)
FeaturePlot(IR_CD4T, features = c("Cd3d","Foxp3","Cd4"))
save(IR_CD8T, file = file.path(data_dir, "IR_CD8T.RData"))
save(IR_CD4T, file = file.path(data_dir, "IR_CD4T.RData"))
save(IR_NK, file = file.path(data_dir, "IR_NK.RData"))
save(IR_B, file = file.path(data_dir, "IR_B.RData"))
save(IR_myeloid, file = file.path(data_dir, "IR_myeloid.RData"))
