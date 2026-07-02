#' scitargets: Single-cell analysis with targets
#'
#' This package provides target factories and utilities for analyzing
#' single-cell data with the \pkg{Seurat} and \pkg{targets} packages.
#'
#' # scitargets Functions
#'
#' ## Function used by the targets factories:
#' 1. [load_seurat_data_10X]: load single-cell data from 10X cellranger output
#' 2. [filter_cell_and_run_reduction]: To be called on the loaded seurat data to filter cells based on their number of genes and percentage of mitochondrial counts.
#' 3. [demultiplex_cell]: demultiplex the cells using the HTO assay
#' 4. [extract_non_negative_cells]: filter the cells to remove the negative cells (keep both singlets and doublets).
#' 5. [extract_singlets]: filter the cells to keep only the singlets. Also allow to remove some features (genes) by names before recomputing UMAP, t-SNE, and clustering the cells.
#' 6. [azimuth_annot_pbmc]: Use azimuth to annotate celltypes using a PBMC reference
#'
#' ## Targets factories:
#'
#' Targets factories are low-level functions used to generate an ensemble of targets steps. For more information about targets factories check the [targetopia contribution section](https://wlandau.github.io/targetopia/contributing.html):
#' 1. [tar_demultiplex_hto]: Used to create a pipeline using the functions 1 to 5 from the previous steps.
#'
#' ## Other targets related functions:
#' 1. [tar_visnetwork_enhanced]: Similar to [targets::tar_visnetwork] but the graph was modified to easier to read.
#'
#' ## Differential expression / GO / GSEA ([scitargets_dea] S7 class):
#' 1. [dea_comparisons]: enumerate the clinical-group comparisons to run (with `groups` and `clusters` selectors).
#' 2. [run_dea]: run one comparison (single-cell `FindMarkers` and/or pseudobulk `DESeq2`) with GO and GSEA.
#' 3. [get_msigdbr_pathways]: fetch MSigDB collections (Hallmark, GO:BP, C7 ImmuneSigDB) for GSEA.
#' 4. Accessors: [markers_table], [go_table], [go_barplot], [go_cnetplot], [go_genes_html], [volcano_plot], [gsea_table], [gsea_barplot].
#' 5. [write_dea_xlsx] / [dea_write_xlsx]: write an xlsx report; [dea_report_lines]: generate Quarto report sections.
#'
#' ## Exported functions:
#' - magrittr pipes: \link[magrittr]{%>%}, \link[magrittr]{%<>%}
#'
#' @docType package
#' @name scitargets
#' @importFrom ggiraph girafe
#' @importFrom ggpubr annotate_figure
#' @importFrom ggrepel geom_text_repel
"_PACKAGE"
