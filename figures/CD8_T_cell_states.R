# Tumor-infiltrating CD8 T-cell states (scRNA-seq): UMAP, marker heatmap, density of control vs
# treated cells, fraction of Tcf7+ / Pdcd1+ / Tox+ cells
#
# Input : IR_CD8T_tumor.RData  (tumor CD8 T cells; clustering at resolution 0.4 and the UMAP
#         are stored in the object, see processing/02_CD8_T_cells.R)
# Output: CD8_UMAP.pdf, CD8_marker_heatmap_all_genes.pdf, CD8_marker_heatmap.pdf,
#         CD8_density_ctrl_vs_treated.pdf, CD8_marker_positive_fraction.pdf / .csv
# Figure panels produced by each block: docs/figure_map.csv

source("config.R")
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

IR_CD8T_tumor <- load_object("IR_CD8T_tumor")
DefaultAssay(IR_CD8T_tumor) <- "RNA"

# Clusters at resolution 0.4 are the six CD8 T-cell states
IR_CD8T_tumor$Yan_annotate <- IR_CD8T_tumor$integrated_snn_res.0.4
Idents(IR_CD8T_tumor) <- IR_CD8T_tumor@meta.data$integrated_snn_res.0.4

# -----------------------------
# State annotation of the six clusters
# -----------------------------
IR_CD8T_tumor$Yan_annotate_chr <- as.character(IR_CD8T_tumor$Yan_annotate)

cluster_map <- c(
  "0" = "C0-Tpex",
  "1" = "C1-Tmem",
  "2" = "C2-Tint",
  "3" = "C3-Ttex",
  "4" = "C4-Teff",
  "5" = "C5-Tpex"
)

IR_CD8T_tumor$state_annot <- unname(cluster_map[IR_CD8T_tumor$Yan_annotate_chr])

IR_CD8T_tumor$state_annot <- factor(
  IR_CD8T_tumor$state_annot,
  levels = c(
    "C1-Tmem",
    "C0-Tpex",
    "C5-Tpex",
    "C2-Tint",
    "C3-Ttex",
    "C4-Teff"
  )
)
print(table(IR_CD8T_tumor$Yan_annotate, IR_CD8T_tumor$state_annot))
print(table(IR_CD8T_tumor$orig.ident, IR_CD8T_tumor$state_annot))

# Same six clusters, named, in cluster-number order (keeps the default cluster colors)
IR_CD8T_tumor$state_by_cluster <- factor(
  unname(cluster_map[IR_CD8T_tumor$Yan_annotate_chr]),
  levels = unname(cluster_map)
)
Idents(IR_CD8T_tumor) <- "state_by_cluster"

# -----------------------------
# UMAP of tumor CD8 T cells
# -----------------------------
umap_plot <- DimPlot(IR_CD8T_tumor, reduction = "umap", pt.size = 2, label = TRUE, label.size = 5)
ggsave(filename = file.path(out_dir, "CD8_UMAP.pdf"), plot = umap_plot, width = 10, height = 8)

# -----------------------------
# Average expression of functional markers per cluster
# -----------------------------
CD8_marker = c("Tcf7","Sell","Lef1","Ccr7","Il7r", "Itgae","Klrg1","Gzma","Gzmb","Gzmk","Prf1", "Ifng","Cd69","Mki67","Cd44","Havcr2", "Entpd1", "Pdcd1","Tox","Ctla4","Lag3", "Tigit","Zfp683","Prdm1","Runx3","Klf2","Cxcr3","Cxcr6")

# Calculate the average expression per cluster
avg_expression <- AverageExpression(IR_CD8T_tumor, return.seurat = TRUE)
avg_exp_matrix <- as.matrix(get_assay_matrix(avg_expression, "RNA", "data")[unique(CD8_marker), ])

# Define breakpoints and corresponding colors
breaks <- c(min(avg_exp_matrix), 1, max(avg_exp_matrix))
colors <- c("blue", "white", "red")

# Create color mapping function
color_mapping <- colorRamp2(breaks, colors)

# All 28 markers, rows and columns clustered
pdf(file.path(out_dir, "CD8_marker_heatmap_all_genes.pdf"), width = 6, height = 6)
draw(Heatmap(avg_exp_matrix,
        cluster_rows = TRUE,
        show_column_names = TRUE,
        column_names_rot = 0,
        column_gap = unit(0, "mm"),
        col = color_mapping,
        name = "Expression"))
dev.off()

# Panel as published: the 19 markers shown in the figure, grouped by function,
# with the cluster order of the figure and the same color mapping as above
panel_genes <- c(
  Memory     = "Klf2", Memory     = "Il7r", Memory     = "Tcf7", Memory = "Lef1", Memory = "Sell", Memory = "Ccr7",
  Activated  = "Cd69", Activated  = "Pdcd1", Activated = "Mki67",
  Effector   = "Gzmb", Effector   = "Gzma", Effector   = "Ifng", Effector = "Prf1",
  Exhaustion = "Tox", Exhaustion  = "Tigit", Exhaustion = "Lag3", Exhaustion = "Ctla4", Exhaustion = "Havcr2", Exhaustion = "Entpd1"
)
panel_cluster_order <- c("C3-Ttex", "C2-Tint", "C4-Teff", "C5-Tpex", "C0-Tpex", "C1-Tmem")

pdf(file.path(out_dir, "CD8_marker_heatmap.pdf"), width = 5, height = 5)
draw(Heatmap(avg_exp_matrix[unname(panel_genes), panel_cluster_order],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = factor(names(panel_genes), levels = unique(names(panel_genes))),
        show_column_names = TRUE,
        column_names_rot = 45,
        col = color_mapping,
        name = "Expression"))
dev.off()

# -----------------------------
# Density of control and treated cells on the UMAP, local and distal tumors
# -----------------------------
# Extract UMAP coordinates from IR_CD8T_tumor
umap_df <- Embeddings(IR_CD8T_tumor, reduction = "umap") %>%
  as.data.frame()

umap_df$cell <- rownames(umap_df)
umap_df$orig.ident <- IR_CD8T_tumor$orig.ident
umap_df$state_annot <- IR_CD8T_tumor$state_annot

# Keep only tumor groups of interest
umap_df <- umap_df %>%
  filter(orig.ident %in% c(
    "local_ctrl_T", "local_IR_T",
    "distal_ctrl_T", "distal_IR_T"
  ))
local_df <- umap_df %>%
  filter(orig.ident %in% c("local_ctrl_T", "local_IR_T")) %>%
  mutate(group = factor(
    orig.ident,
    levels = c("local_ctrl_T", "local_IR_T"),
    labels = c("Local Ctrl", "Local IR")
  ))

p_local <- ggplot(local_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    color = "grey88",
    size = 0.35,
    alpha = 0.35
  ) +
  geom_density_2d(
    data = subset(local_df, group == "Local Ctrl"),
    aes(color = group),
    linewidth = 0.7,
    bins = 4
  ) +
  geom_density_2d(
    data = subset(local_df, group == "Local IR"),
    aes(color = group),
    linewidth = 0.7,
    bins = 4
  ) +
  scale_color_manual(
    values = c(
      "Local Ctrl" = "#A0A0A0",
      "Local IR" = "#D6604D"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Local tumor CD8 T cells",
    x = "UMAP_1",
    y = "UMAP_2"
  )

distal_df <- umap_df %>%
  filter(orig.ident %in% c("distal_ctrl_T", "distal_IR_T")) %>%
  mutate(group = factor(
    orig.ident,
    levels = c("distal_ctrl_T", "distal_IR_T"),
    labels = c("Distal Ctrl", "Distal IR")
  ))

p_distal <- ggplot(distal_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    color = "grey88",
    size = 0.35,
    alpha = 0.35
  ) +
  geom_density_2d(
    data = subset(distal_df, group == "Distal Ctrl"),
    aes(color = group),
    linewidth = 0.7,
    bins = 4
  ) +
  geom_density_2d(
    data = subset(distal_df, group == "Distal IR"),
    aes(color = group),
    linewidth = 0.7,
    bins = 4
  ) +
  scale_color_manual(
    values = c(
      "Distal Ctrl" = "#A0A0A0",
      "Distal IR" = "#D6604D"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Distal tumor CD8 T cells",
    x = "UMAP_1",
    y = "UMAP_2"
  )

label_df <- umap_df %>%
  group_by(state_annot) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )
p_local_label <- p_local +
  geom_label(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = state_annot),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold",
    fill = "white",
    alpha = 0.85,
    label.size = 0.2
  )

p_distal_label <- p_distal +
  geom_label(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = state_annot),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold",
    fill = "white",
    alpha = 0.85,
    label.size = 0.2
  )

p_density_two_label <- p_local_label + p_distal_label + plot_layout(ncol = 2)
ggsave(file.path(out_dir, "CD8_density_ctrl_vs_treated.pdf"), p_density_two_label, width = 10.4, height = 3.6)

# -----------------------------
# Fraction of Tcf7+, Pdcd1+ and Tox+ cells among tumor CD8 T cells
# -----------------------------
markers_use <- c("Tcf7", "Pdcd1", "Tox")

marker_df <- FetchData(
  IR_CD8T_tumor,
  vars = c(markers_use, "orig.ident")
)

marker_df$cell <- rownames(marker_df)
threshold <- 0.5

marker_long <- marker_df %>%
  mutate(
    site = case_when(
      grepl("^local", orig.ident) ~ "Local",
      grepl("^distal", orig.ident) ~ "Distal",
      TRUE ~ NA_character_
    ),
    treatment = case_when(
      grepl("_ctrl_", orig.ident) ~ "Ctrl",
      grepl("_IR_", orig.ident) ~ "Tx",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(site), !is.na(treatment)) %>%
  pivot_longer(
    cols = all_of(markers_use),
    names_to = "marker",
    values_to = "expression"
  ) %>%
  mutate(
    marker_pos = expression > threshold,
    site = factor(site, levels = c("Local", "Distal")),
    treatment = factor(treatment, levels = c("Ctrl", "Tx")),
    marker = factor(marker, levels = c("Tcf7", "Pdcd1", "Tox"))
  )
marker_plot_df <- marker_long %>%
  group_by(site, treatment, marker) %>%
  summarise(
    n_cells = n(),
    n_marker_pos = sum(marker_pos),
    pct_marker_pos = 100 * n_marker_pos / n_cells,
    .groups = "drop"
  )

print(marker_plot_df, n = Inf)
write.csv(marker_plot_df, file.path(out_dir, "CD8_marker_positive_fraction.csv"), row.names = FALSE)

p_marker_pos <- ggplot(
  marker_plot_df,
  aes(x = treatment, y = pct_marker_pos, fill = treatment)
) +
  geom_col(
    width = 0.65,
    color = "black",
    linewidth = 0.3
  ) +
  facet_grid(marker ~ site) +
  scale_fill_manual(
    values = c(
      "Ctrl" = "#A0A0A0",
      "Tx" = "#D6604D"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    y = "% marker+ tumor CD8 T cells"
  )

ggsave(file.path(out_dir, "CD8_marker_positive_fraction.pdf"), p_marker_pos, width = 3.6, height = 3.3)
