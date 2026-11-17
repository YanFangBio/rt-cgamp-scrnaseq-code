# Shared settings for all scripts in this repository.
# Run every script from the repository root, e.g.  Rscript figures/CD8_T_cell_states.R
#
# data_dir : folder that holds the processed Seurat objects (*.RData)
# out_dir  : folder that receives the figure panels (PDF) and tables (CSV)
data_dir <- Sys.getenv("SCRNA_DATA_DIR", unset = "data")
out_dir  <- Sys.getenv("SCRNA_OUT_DIR",  unset = "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load one Seurat object saved with save(<name>, file = "<name>.RData")
load_object <- function(name) {
  e <- new.env()
  load(file.path(data_dir, paste0(name, ".RData")), envir = e)
  get(name, envir = e)
}

# Read an assay matrix with either SeuratObject >= 5 (layer) or < 5 (slot)
get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  if (packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(obj, assay = assay, layer = layer)
  } else {
    GetAssayData(obj, assay = assay, slot = layer)
  }
}

# Fold-change definition of Seurat 4.1.0 (the version the 2023 analysis was run with):
# avg_log2FC = log2(mean(expm1(x)) + 1) of group 1 minus the same for group 2.
# Seurat >= 5 changed its default, which re-orders markers ranked by avg_log2FC. Passing this
# function as `mean.fxn` gives the Seurat 4.1.0 values in any Seurat version.
mean_fxn_v4 <- function(x) log(x = rowMeans(x = expm1(x = x)) + 1, base = 2)

# Fold-change definition of Seurat >= 5 (the default there, written out so that it can be passed
# explicitly): avg_log2FC = log2((sum(expm1(x)) + 1) / number of cells) of group 1 minus the same
# for group 2. Used for the Reactome GSEA ranking, which was run with Seurat 5.
mean_fxn_v5 <- function(x) log(x = (rowSums(x = expm1(x = x)) + 1) / NCOL(x), base = 2)
