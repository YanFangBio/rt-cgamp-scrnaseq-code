# Tumor-infiltrating CD8 T cells: cluster composition per sample, TCR clonal expansion score (CS)
# per sample, and fraction of clonally expanded cells per state and sample
#
# Input : IR_CD8T_tumor.RData  (meta.data columns used: orig.ident, integrated_snn_res.0.4, tcr_cluster,
#         Expand_status)
# Output: CD8_cluster_composition.csv / .pdf, CD8_clonal_expansion_score.csv / .pdf,
#         CD8_expanded_fraction.csv / .pdf
# Figure panels produced by each block: docs/figure_map.csv
#
# tcr_cluster   : TCR-beta CDR3 cluster of the cell (NA = no TCR detected)
# Expand_status : "Single" = clonotype seen in one CD8 T cell of the dataset (tumor + lymph node);
#                 every other level = clonotype shared by two or more CD8 T cells (expanded)

source("config.R")
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

IR_CD8T_tumor <- load_object("IR_CD8T_tumor")

IR_CD8T_tumor$Yan_annotate <- IR_CD8T_tumor$integrated_snn_res.0.4
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

# -----------------------------
# Cluster composition per sample
# -----------------------------
composition_df <- IR_CD8T_tumor@meta.data %>%
  as.data.frame() %>%
  mutate(
    sample_group = orig.ident,
    site = case_when(
      grepl("^local", sample_group) ~ "Local",
      grepl("^distal", sample_group) ~ "Distal",
      TRUE ~ NA_character_
    ),
    treatment = case_when(
      grepl("_ctrl_", sample_group) ~ "Ctrl",
      grepl("_IR_", sample_group) ~ "IR",
      TRUE ~ NA_character_
    )
  ) %>%
  count(sample_group, site, treatment, state_annot, name = "n_cells") %>%
  group_by(sample_group, site, treatment) %>%
  mutate(
    total_CD8T = sum(n_cells),
    percent = 100 * n_cells / total_CD8T
  ) %>%
  ungroup()

print(composition_df, n = Inf)
write.csv(
  composition_df,
  file.path(out_dir, "CD8_cluster_composition.csv"),
  row.names = FALSE
)
state_colors <- c(
  "C1-Tmem" = "#B89B00",
  "C0-Tpex" = "#FF6B6B",
  "C5-Tpex" = "#C85DBA",
  "C2-Tint" = "#00B050",
  "C3-Ttex" = "#00BFC4",
  "C4-Teff" = "#6A9BD8"
)

composition_df$sample_group <- factor(
  composition_df$sample_group,
  levels = c(
    "local_ctrl_T",
    "local_IR_T",
    "distal_ctrl_T",
    "distal_IR_T"
  )
)

p_stack <- ggplot(
  composition_df,
  aes(x = sample_group, y = percent, fill = state_annot)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    y = "% of tumor CD8 T cells",
    fill = "CD8 T state"
  )
ggsave(file.path(out_dir, "CD8_cluster_composition.pdf"), p_stack, width = 5, height = 4)

# -----------------------------
# Clonal expansion score (CS) = 1 - normalized Shannon entropy of the clonotype frequencies
# -----------------------------
# Function to calculate Shannon entropy
shannon_entropy <- function(probabilities) {
  # Filter out probabilities that are 0 because log(0) is undefined
  valid_probs <- probabilities[probabilities > 0]

  # Calculate entropy
  entropy <- -sum(valid_probs * log2(valid_probs))

  return(entropy)
}

# Function to calculate normalized Shannon entropy
normalized_shannon_entropy <- function(probabilities) {
  n <- length(probabilities)
  entropy <- shannon_entropy(probabilities)

  return(entropy / log2(n))
}

sample_order <- c("local_ctrl_T", "local_IR_T", "distal_ctrl_T", "distal_IR_T")
cs_df <- do.call(rbind, lapply(sample_order, function(x) {
  tcr_cluster = IR_CD8T_tumor@meta.data[which(IR_CD8T_tumor@meta.data$tcr_cluster != 'NA' & IR_CD8T_tumor@meta.data$orig.ident == x),]
  a = table(tcr_cluster$tcr_cluster)
  a = as.data.frame(a)
  colnames(a) = c("tcr_cluster", "freq")
  sorted_freq=a[order(a$freq,decreasing=T),]
  total_Tcell = sum(sorted_freq$freq)
  probabilities = (sorted_freq$freq)/total_Tcell
  data.frame(sample_group = x, n_cells_with_TCR = total_Tcell, n_clonotypes = nrow(sorted_freq),
             clonal_expansion_score = 1-normalized_shannon_entropy(probabilities))
}))
cs_df$site <- ifelse(grepl("^local", cs_df$sample_group), "Local", "Distal")
cs_df$treatment <- ifelse(grepl("_ctrl_", cs_df$sample_group), "Ctrl", "Tx")
print(cs_df)
write.csv(cs_df, file.path(out_dir, "CD8_clonal_expansion_score.csv"), row.names = FALSE)

p_cs <- ggplot(cs_df, aes(x = treatment, y = clonal_expansion_score, fill = treatment)) +
  geom_col(width = 0.65) +
  facet_wrap(~ factor(site, levels = c("Local", "Distal")), scales = "free_y") +
  scale_fill_manual(values = c("Ctrl" = "#2E3192", "Tx" = "#ED1C24")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.title.x = element_blank(), strip.background = element_blank()) +
  labs(y = "Clonal expansion (CS)")
ggsave(file.path(out_dir, "CD8_clonal_expansion_score.pdf"), p_cs, width = 5, height = 3)

# -----------------------------
# Fraction of clonally expanded cells per state and sample (expanded T / total T)
# -----------------------------
expand_df <- IR_CD8T_tumor@meta.data %>%
  as.data.frame() %>%
  mutate(expanded = Expand_status != "Single") %>%
  group_by(orig.ident, state_annot) %>%
  summarise(n_cells = n(), n_expanded = sum(expanded), expanded_fraction = n_expanded / n_cells, .groups = "drop") %>%
  mutate(
    site = factor(ifelse(grepl("^local", orig.ident), "Local", "Distal"), levels = c("Local", "Distal")),
    treatment = factor(ifelse(grepl("_ctrl_", orig.ident), "Ctrl", "Tx"), levels = c("Ctrl", "Tx"))
  )
print(expand_df, n = Inf)
write.csv(expand_df, file.path(out_dir, "CD8_expanded_fraction.csv"), row.names = FALSE)

p_expand <- ggplot(expand_df, aes(x = treatment, y = expanded_fraction, color = state_annot, group = state_annot)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ site) +
  scale_color_manual(values = state_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic(base_size = 13) +
  theme(axis.title.x = element_blank(), strip.background = element_blank()) +
  labs(y = "Expanded T / total T", color = NULL, title = "TCR expand")
ggsave(file.path(out_dir, "CD8_expanded_fraction.pdf"), p_expand, width = 5.5, height = 3.5)
