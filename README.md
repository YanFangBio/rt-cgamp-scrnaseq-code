# scRNA-seq analysis code

Code for the single-cell RNA-seq and TCR analyses in

> cGAMP-Driven cDC1 Maturation and Migration Promote Systemic Antitumor Immunity

The repository holds (i) the scripts that draw the scRNA-seq figure panels from the processed
Seurat objects and (ii) the processing record that produced those objects from the integrated data set.

## Layout

```
config.R                  data / output folders and shared helper functions
figures/                  one script per analysis; each runs on its own
  CD8_T_cell_states.R               tumor CD8 T cells: UMAP, marker heatmap, density, marker-positive fraction
  CD8_composition_TCR_expansion.R   tumor CD8 T cells: cluster composition, clonal expansion score, expanded fraction
  DC_subsets.R                      dendritic cells: UMAP, top-10 markers, ridge plots, cDC1 / cDC2 genes
  DC_feature_plots_Mafb_S100a8.R    dendritic cells: Mafb and S100a8 feature plots
  DC_reactome_gsea.R                dendritic cells: Reactome pathway GSEA, mature DC vs cDC1
  CD4_T_cells.R                     tumor CD4 T cells: UMAP, cluster composition, marker heatmap
  ligand_receptor_dotplot.R         Xcl1-Xcr1, Ccl5-Ccr5, Flt3l-Flt3 across cell populations
processing/               processing record (as run in 2023; see "Processing record")
  01_filter_cluster_split_lineages.R
  02_CD8_T_cells.R
  03_CD4_T_cells.R
  04_myeloid_DC.R
  05_NK_cells.R
  06_B_cells.R
docs/figure_map.csv       which script and output file correspond to which figure panel
LICENSE                   MIT License
```

## Data

Raw sequencing data: NCBI BioProject PRJNA1469973.

The figure scripts read processed Seurat objects saved as `<name>.RData` (one object per file,
object name = file name):

| Object | Content | Read by |
|---|---|---|
| `IR_CD8T_tumor` | tumor CD8 T cells (3,139 cells) | `CD8_T_cell_states.R`, `CD8_composition_TCR_expansion.R` |
| `IR_CD4T_tumor` | tumor CD4 T cells (1,900 cells) | `CD4_T_cells.R` |
| `IR_myeloid_DC` | dendritic cells (858 cells) | `DC_subsets.R`, `DC_feature_plots_Mafb_S100a8.R`, `DC_reactome_gsea.R`, `ligand_receptor_dotplot.R` |
| `IR_CD8T`, `IR_CD4T`, `IR_B`, `IR_NK` | lineage objects, tumor and lymph node | `ligand_receptor_dotplot.R` |

The objects are not stored in this repository and are not deposited elsewhere; they can be obtained
from the lead contact, as stated in the Data availability section of the article. Put them in one
folder and point `SCRNA_DATA_DIR` to it (default: `data/`).

Sample names (`orig.ident`): `local_ctrl_T`, `local_IR_T`, `distal_ctrl_T`, `distal_IR_T` (tumors) and the
matching `_LN` samples (lymph nodes). `IR` denotes the treated group (RT + cGAMP + anti-PD-L1; "Tx" in the figures).

Cluster names used in the figures:

| Object | Cluster (resolution 0.4) | Name |
|---|---|---|
| `IR_CD8T_tumor` | 0, 1, 2, 3, 4, 5 | C0-Tpex, C1-Tmem, C2-Tint, C3-Ttex, C4-Teff, C5-Tpex |
| `IR_CD4T_tumor` | 0, 1, 2, 3, 4, 5 | C0-Tnaive, C1-Treg, C2-Tint, C3-Treg, C4-Th1, C5-Unknown |
| `IR_myeloid_DC` | `Yan_annotate` | cDC1, cDC2, CCR7+ (shown as "Mature DC") |

## Running

```
export SCRNA_DATA_DIR=/path/to/objects     # folder with the .RData files
export SCRNA_OUT_DIR=output                # created if missing
Rscript figures/CD8_T_cell_states.R
Rscript figures/CD8_composition_TCR_expansion.R
Rscript figures/DC_subsets.R
Rscript figures/DC_feature_plots_Mafb_S100a8.R
Rscript figures/DC_reactome_gsea.R
Rscript figures/CD4_T_cells.R
Rscript figures/ligand_receptor_dotplot.R
```

Run the scripts from the repository root. Each writes PDF panels and, where a panel is a summary,
the underlying table as CSV. `ligand_receptor_dotplot.R` loads the lineage objects (1.8-2.4 GB each
on disk) one at a time and removes each before loading the next.

## Software

R packages: Seurat, dplyr, tidyr, ggplot2, patchwork, stringr, ComplexHeatmap, circlize; for
`DC_reactome_gsea.R` also clusterProfiler, ReactomePA (with reactome.db) and org.Mm.eg.db.

| | Original analysis | Scripts re-run for this repository |
|---|---|---|
| R | 4.1.1 | 4.5.3 |
| Seurat | 4.1.0 | 5.5.1 (SeuratObject 5.4.0) |
| ComplexHeatmap / circlize | | 2.26.1 / 0.4.18 |
| ggplot2 / dplyr / tidyr / patchwork | | 4.0.3 / 1.2.1 / 1.3.2 / 1.3.2 |
| clusterProfiler / ReactomePA / fgsea | | 4.18.4 / 1.54.0 / 1.36.2 |
| reactome.db / org.Mm.eg.db | | 1.95.0 / 3.22.0 |

Three points matter when different package versions are used:

- **Marker ranking.** Seurat 5 changed the default definition of `avg_log2FC`. `DC_subsets.R` passes
  `mean.fxn = mean_fxn_v4` (defined in `config.R`), the Seurat 4.1.0 definition, so the top-10 marker
  genes are the same in any version.
- **Pathway enrichment.** `DC_reactome_gsea.R` ranks genes by `avg_log2FC` computed with
  `mean.fxn = mean_fxn_v5` (the Seurat 5 definition, written out in `config.R`), so the ranking does not
  depend on the Seurat version. The pathways tested come from the installed `reactome.db` release; a
  different release can change which pathways pass the cutoff.
- **Stored results.** Clustering and UMAP coordinates are read from the objects and are not recomputed
  by the figure scripts.

## Figure panels that are not produced by code

Cluster labels on the UMAPs, the "Effector-like / Stem-like" brackets, the renaming of "CCR7+" to
"Mature DC" and of "Tumor_IR" to "Tumor_Tx", and the connector lines in the ligand-receptor plot were
added during figure assembly. The cluster-composition bar charts and the TCR summary plots in the
figures were drawn separately from the values these scripts compute and write as CSV.

## Processing record

`processing/01-06` are the commands that were run interactively on an HPC cluster to go from the
integrated object to the lineage objects and the tumor subsets. They are kept as run, with absolute
paths replaced by `data_dir` / `upstream_dir`, and were not re-executed for this repository.
Exploratory plotting calls are left in place next to the step they informed.

- Three clusters (18, 21, 27) were split by hand with `Seurat::CellSelector()`; this step is
  interactive and its result is stored in the objects.
- Upstream inputs that are not created here: `TIL_seurat_bf_cluster.RData` (cells after quality
  filtering, DoubletFinder labelling and integration, as described in the Methods) and
  `TCR_cdr3b_RotationEncodingBL62.txt` (TCR-beta CDR3 clusters per cell barcode).

## License

MIT License, see `LICENSE`.
