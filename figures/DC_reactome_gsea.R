# Dendritic cells (scRNA-seq): Reactome pathway enrichment (GSEA) in mature DCs (CCR7+) versus cDC1
#
# Genes are ranked by avg_log2FC (CCR7+ vs cDC1, all genes detected in >= 10% of cells of either group,
# no fold-change threshold) and tested with ReactomePA::gsePathway (Benjamini-Hochberg adjusted P < 0.1).
#
# Input : IR_myeloid_DC.RData
# Output: DC_reactome_gsea.pdf, DC_reactome_gsea_results.csv, DC_reactome_gsea_ranked_genes.csv
# Run from the repository root: Rscript figures/DC_reactome_gsea.R

source("config.R")
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)

IR_myeloid_DC <- load_object("IR_myeloid_DC")
DefaultAssay(IR_myeloid_DC) <- "RNA"
Idents(IR_myeloid_DC) <- IR_myeloid_DC@meta.data$Yan_annotate

# 1) Differential expression, mature DC (CCR7+) vs cDC1
#    mean.fxn = mean_fxn_v5 fixes the fold-change definition (Seurat >= 5 default, see config.R)
markers_CCR7_vs_cDC1 <- FindMarkers(
  IR_myeloid_DC,
  ident.1 = "CCR7+",
  ident.2 = "cDC1",
  logfc.threshold = 0,
  min.pct = 0.1,
  test.use = "wilcox",
  mean.fxn = mean_fxn_v5
)

# 2) Ranked gene list (SYMBOL -> avg_log2FC)
gene_list <- markers_CCR7_vs_cDC1$avg_log2FC
names(gene_list) <- rownames(markers_CCR7_vs_cDC1)
gene_list <- sort(gene_list, decreasing = TRUE)

# 3) Map SYMBOL -> ENTREZID
gene_map <- bitr(
  names(gene_list),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

# 4) ENTREZ-named ranked vector; duplicated IDs collapsed to the value with the largest abs(log2FC)
gene_list_mapped <- gene_list[gene_map$SYMBOL]
names(gene_list_mapped) <- gene_map$ENTREZID

gene_list2 <- tapply(
  gene_list_mapped,
  INDEX = names(gene_list_mapped),
  FUN = function(x) x[which.max(abs(x))]
)
gene_list2 <- sort(setNames(as.numeric(gene_list2), names(gene_list2)), decreasing = TRUE)

ranked_genes <- data.frame(
  ENTREZID = names(gene_list2),
  SYMBOL = gene_map$SYMBOL[match(names(gene_list2), gene_map$ENTREZID)],
  avg_log2FC = as.numeric(gene_list2)
)
# every gene must carry its own fold change
stopifnot(isTRUE(all.equal(ranked_genes$avg_log2FC, unname(gene_list[ranked_genes$SYMBOL]))))
write.csv(ranked_genes, file.path(out_dir, "DC_reactome_gsea_ranked_genes.csv"), row.names = FALSE)

# 5) Reactome GSEA
#    nPermSimple = 100000 (fgsea default: 1000) makes the P values precise enough that the pathways
#    passing the cutoff do not depend on the random seed.
set.seed(1)
gsea_reactome <- gsePathway(
  geneList = gene_list2,
  organism = "mouse",
  minGSSize = 10,
  pvalueCutoff = 0.1,
  verbose = FALSE,
  seed = TRUE,
  nPermSimple = 100000
)
gsea_table <- as.data.frame(setReadable(gsea_reactome, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")@result)
write.csv(gsea_table, file.path(out_dir, "DC_reactome_gsea_results.csv"), row.names = FALSE)

# 6) Plot: top pathways in each direction (positive NES = enriched in mature DC, negative = in cDC1)
plot_reactome_manuscript <- function(df,
                                     title = "Reactome",
                                     n_pos = 10,
                                     n_neg = 10) {

  plot_df <- as.data.frame(df) %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    mutate(sig = -log10(p.adjust))

  top_pos <- plot_df %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    slice_head(n = n_pos)

  top_neg <- plot_df %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    slice_head(n = n_neg)

  plot_df <- bind_rows(top_pos, top_neg) %>%
    mutate(
      short_label = str_trunc(Description, width = 42),
      direction = ifelse(NES > 0, "positive", "negative")
    ) %>%
    arrange(NES) %>%
    mutate(short_label = factor(short_label, levels = short_label))

  ggplot(plot_df, aes(x = NES, y = short_label)) +
    geom_segment(
      aes(x = 0, xend = NES, yend = short_label),
      linewidth = 0.4,
      color = "black"
    ) +
    geom_point(
      aes(size = sig, fill = direction),
      shape = 21,
      color = "black",
      stroke = 0.2
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "positive" = "#E64B35",
        "negative" = "#6C5CE7"
      )
    ) +
    scale_size_continuous(
      range = c(2.5, 6),
      limits = c(1, 5),
      breaks = c(1, 2, 4)
    ) +
    labs(
      title = title,
      x = "NES",
      y = NULL,
      size = expression(-log[10]("FDR"))
    ) +
    coord_cartesian(xlim = c(-2.4, 2.4)) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.title.x = element_text(size = 11),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    ) +
    guides(fill = "none")
}

p_reactome <- plot_reactome_manuscript(gsea_table, title = "Reactome", n_pos = 10, n_neg = 10)

ggsave(file.path(out_dir, "DC_reactome_gsea.pdf"), p_reactome, width = 6.8, height = 4.2)
