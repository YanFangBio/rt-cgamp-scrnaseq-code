# Candidate recruitment / retention signals (Xcl1-Xcr1, Ccl5-Ccr5, Flt3l-Flt3) by cell population,
# site, tissue and treatment
#
# Input : IR_B.RData, IR_CD4T.RData, IR_CD8T.RData, IR_myeloid_DC.RData, IR_NK.RData
# Output: ligand_receptor/ (summary tables and the dot plot, PDF and PNG)
# Figure panel: docs/figure_map.csv
#
# The published panel shows the two tumor rows of this plot (Tumor_ctrl, Tumor_IR). In the figure,
# "Tumor_IR" is labelled "Tumor_Tx" and "CCR7+" is labelled "Mature DC"; the red connector lines
# between ligand and receptor were drawn by hand.
#
# Objects are loaded one at a time, summarized, and removed, so the objects are never merged.

source("config.R")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ----------------------------
# 1. Paths
# ----------------------------

base_dir <- data_dir

rdata_files <- c(
  IR_B          = file.path(base_dir, "IR_B.RData"),
  IR_CD4T       = file.path(base_dir, "IR_CD4T.RData"),
  IR_CD8T       = file.path(base_dir, "IR_CD8T.RData"),
  IR_myeloid_DC = file.path(base_dir, "IR_myeloid_DC.RData"),
  IR_NK         = file.path(base_dir, "IR_NK.RData")
)

out_dir <- file.path(out_dir, "ligand_receptor")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genes_to_show <- c(
  "Xcl1",
  "Xcl2",
  "Ccl5",
  "Flt3l",
  "Xcr1",
  "Ccr5",
  "Flt3"
)

celltype_order <- c("NK", "CD8 T", "CD4 T", "B", "cDC1", "cDC2", "CCR7+")

sample_group_order <- c(
  "local_Tumor_ctrl",
  "local_Tumor_IR",
  "distal_Tumor_ctrl",
  "distal_Tumor_IR",
  "local_LN_ctrl",
  "local_LN_IR",
  "distal_LN_ctrl",
  "distal_LN_IR"
)

# ============================================================
# Helper functions
# ============================================================

load_seurat_from_rdata <- function(path) {
  e <- new.env()
  load(path, envir = e)
  
  obj_names <- ls(e)
  seurat_names <- obj_names[
    sapply(obj_names, function(nm) inherits(get(nm, envir = e), "Seurat"))
  ]
  
  if (length(seurat_names) == 0) {
    stop("No Seurat object found in: ", path)
  }
  
  obj <- get(seurat_names[1], envir = e)
  
  rm(e)
  gc()
  
  return(obj)
}


get_assay_data_compat <- function(obj, assay = "RNA", layer = "data") {
  DefaultAssay(obj) <- assay
  
  if (packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(obj, assay = assay, layer = layer))
  } else {
    return(GetAssayData(obj, assay = assay, slot = layer))
  }
}


ensure_rna_data_layer <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  
  data_mat <- tryCatch(
    get_assay_data_compat(obj, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  
  if (is.null(data_mat) || nrow(data_mat) == 0 || ncol(data_mat) == 0) {
    message("RNA data layer not found. Running NormalizeData().")
    obj <- NormalizeData(obj, assay = "RNA")
  }
  
  return(obj)
}


parse_sample_metadata <- function(meta) {
  meta %>%
    mutate(
      site = case_when(
        str_detect(orig.ident, "^local") ~ "local",
        str_detect(orig.ident, "^distal") ~ "distal",
        TRUE ~ "unknown"
      ),
      treatment = case_when(
        str_detect(orig.ident, "_ctrl_") ~ "ctrl",
        str_detect(orig.ident, "_IR_") ~ "IR",
        TRUE ~ "unknown"
      ),
      tissue_type = case_when(
        str_detect(orig.ident, "_T$") ~ "Tumor",
        str_detect(orig.ident, "_LN$") ~ "LN",
        TRUE ~ "unknown"
      ),
      sample_group = paste(site, tissue_type, treatment, sep = "_")
    )
}
# ============================================================
# Process one Seurat object
# ============================================================

process_one_object <- function(object_name, file_path, genes_to_show) {
  
  message("Processing: ", object_name)
  
  obj <- load_seurat_from_rdata(file_path)
  obj <- ensure_rna_data_layer(obj)
  
  DefaultAssay(obj) <- "RNA"
  
  # ----------------------------
  # Assign cell type
  # ----------------------------
  
  if (object_name == "IR_NK") {
    obj$plot_celltype <- "NK"
  } else if (object_name == "IR_CD8T") {
    obj$plot_celltype <- "CD8 T"
  } else if (object_name == "IR_CD4T") {
    obj$plot_celltype <- "CD4 T"
  } else if (object_name == "IR_B") {
    obj$plot_celltype <- "B"
  } else if (object_name == "IR_myeloid_DC") {
    
    if (!"Yan_annotate" %in% colnames(obj@meta.data)) {
      stop("IR_myeloid_DC does not have Yan_annotate.")
    }
    
    obj$plot_celltype <- as.character(obj$Yan_annotate)
    
    obj <- subset(
      obj,
      subset = plot_celltype %in% c("cDC1", "cDC2", "CCR7+")
    )
  } else {
    stop("Unknown object name: ", object_name)
  }
  
  # ----------------------------
  # Gene presence
  # ----------------------------
  
  gene_presence <- data.frame(
    object = object_name,
    celltype = paste(unique(obj$plot_celltype), collapse = ";"),
    gene = genes_to_show,
    present = genes_to_show %in% rownames(obj),
    n_cells_object = ncol(obj),
    stringsAsFactors = FALSE
  )
  
  genes_present <- intersect(genes_to_show, rownames(obj))
  
  if (length(genes_present) == 0) {
    warning("No requested genes present in ", object_name)
    
    rm(obj)
    gc()
    
    return(list(
      summary = data.frame(),
      gene_presence = gene_presence
    ))
  }
  
  # ----------------------------
  # Extract only requested genes
  # This is tiny: genes x cells
  # ----------------------------
  
  expr_mat <- get_assay_data_compat(
    obj,
    assay = "RNA",
    layer = "data"
  )[genes_present, , drop = FALSE]
  
  # ----------------------------
  # Metadata
  # ----------------------------
  
  meta <- obj@meta.data %>%
    mutate(cell = rownames(.)) %>%
    parse_sample_metadata() %>%
    select(
      cell,
      orig.ident,
      site,
      tissue_type,
      treatment,
      sample_group,
      plot_celltype
    )
  
  # ----------------------------
  # Summarize
  # ----------------------------
  
  summary_df <- lapply(genes_present, function(gene) {
    data.frame(
      cell = colnames(expr_mat),
      gene = gene,
      expr = as.numeric(expr_mat[gene, ]),
      stringsAsFactors = FALSE
    ) %>%
      left_join(meta, by = "cell") %>%
      filter(
        !is.na(plot_celltype),
        !is.na(sample_group),
        sample_group %in% sample_group_order
      ) %>%
      group_by(
        sample_group,
        site,
        tissue_type,
        treatment,
        plot_celltype,
        gene
      ) %>%
      summarise(
        n_cells = n(),
        n_expr = sum(expr > 0),
        pct_expr = 100 * mean(expr > 0),
        avg_expr = mean(expr),
        avg_expr_positive = ifelse(sum(expr > 0) > 0, mean(expr[expr > 0]), 0),
        .groups = "drop"
      )
  }) %>%
    bind_rows() %>%
    mutate(object = object_name) %>%
    select(object, everything())
  
  # Free memory before returning
  rm(obj, expr_mat, meta)
  gc()
  
  return(list(
    summary = summary_df,
    gene_presence = gene_presence
  ))
}

# ============================================================
# Run one object at a time
# ============================================================

all_results <- list()
all_gene_presence <- list()

for (nm in names(rdata_files)) {
  result <- process_one_object(
    object_name = nm,
    file_path = rdata_files[[nm]],
    genes_to_show = genes_to_show
  )
  
  all_results[[nm]] <- result$summary
  all_gene_presence[[nm]] <- result$gene_presence
  
  gc()
}

plot_summary <- bind_rows(all_results)
gene_presence <- bind_rows(all_gene_presence)

# Clean intermediate lists
rm(all_results, all_gene_presence)
gc()

print(gene_presence)
print(plot_summary)

write.csv(
  gene_presence,
  file.path(out_dir, "gene_presence_memory_safe_no_merge.csv"),
  row.names = FALSE
)

write.csv(
  plot_summary,
  file.path(out_dir, "expression_summary_memory_safe_no_merge.csv"),
  row.names = FALSE
)

# ============================================================
# Complete missing combinations for plotting
# ============================================================

plot_summary <- plot_summary %>%
  mutate(
    gene = factor(gene, levels = genes_to_show),
    plot_celltype = factor(plot_celltype, levels = celltype_order),
    sample_group = factor(sample_group, levels = sample_group_order),
    facet_row = paste(tissue_type, treatment, sep = "_"),
    facet_row = factor(
      facet_row,
      levels = c("Tumor_ctrl", "Tumor_IR", "LN_ctrl", "LN_IR")
    ),
    site = factor(site, levels = c("local", "distal"))
  )

# Only complete combinations that are biologically possible from observed groups
plot_summary_complete <- plot_summary %>%
  complete(
    sample_group,
    site,
    facet_row,
    plot_celltype,
    gene,
    fill = list(
      n_cells = 0,
      n_expr = 0,
      pct_expr = 0,
      avg_expr = 0,
      avg_expr_positive = 0
    )
  ) %>%
  filter(!is.na(site), !is.na(facet_row))

print(plot_summary_complete)

# ============================================================
# Save small objects only
# ============================================================

save(
  plot_summary,
  plot_summary_complete,
  gene_presence,
  file = file.path(out_dir, "single_dotplot_axis_summary_tables_only.RData")
)

message("Done. Results saved to: ", out_dir)

# ============================================================
# Reorder genes and remove Xcl2 from existing plot_summary
# ============================================================

gene_order <- c("Xcl1", "Xcr1", "Ccl5", "Ccr5", "Flt3l", "Flt3")

plot_summary2 <- plot_summary %>%
  filter(gene %in% gene_order) %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    plot_celltype = factor(
      plot_celltype,
      levels = c("NK", "CD8 T", "CD4 T", "B", "cDC1", "cDC2", "CCR7+")
    ),
    facet_row = paste(tissue_type, treatment, sep = "_"),
    facet_row = factor(
      facet_row,
      levels = c("Tumor_ctrl", "Tumor_IR", "LN_ctrl", "LN_IR")
    ),
    site = factor(site, levels = c("local", "distal"))
  )

# Optional: complete missing combinations so empty groups still show as blank/zero
plot_summary2_complete <- plot_summary2 %>%
  tidyr::complete(
    site,
    facet_row,
    plot_celltype,
    gene,
    fill = list(
      n_cells = 0,
      n_expr = 0,
      pct_expr = 0,
      avg_expr = 0,
      avg_expr_positive = 0
    )
  ) %>%
  filter(!is.na(site), !is.na(facet_row))

# ============================================================
# Single dot plot with requested gene order
# ============================================================

p_grid_ordered <- ggplot(
  plot_summary2_complete,
  aes(
    x = gene,
    y = plot_celltype,
    size = pct_expr,
    color = avg_expr
  )
) +
  geom_point() +
  facet_grid(
    facet_row ~ site
  ) +
  scale_size(
    range = c(0, 6),
    limits = c(0, max(plot_summary2_complete$pct_expr, na.rm = TRUE))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "cDC1 recruitment / retention axis by site, tissue, and treatment",
    subtitle = "Dot size = % expressing cells; color = average normalized expression",
    x = "Gene",
    y = "Cell population",
    size = "% positive",
    color = "Avg expression"
  )

print(p_grid_ordered)

ggsave(
  filename = file.path(out_dir, "ligand_receptor_dotplot.pdf"),
  plot = p_grid_ordered,
  width = 11,
  height = 10
)

ggsave(
  filename = file.path(out_dir, "ligand_receptor_dotplot.png"),
  plot = p_grid_ordered,
  width = 11,
  height = 10,
  dpi = 300
)
