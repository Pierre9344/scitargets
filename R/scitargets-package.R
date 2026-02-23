#' scitargets: Single-cell CITE-seq analysis with targets
#'
#' This package provides target factories and utilities for analyzing
#' single-cell [CITE-seq](https://en.wikipedia.org/wiki/CITE-Seq) data with the \pkg{Seurat} and \pkg{targets} packages.
#'
#' # scitargets Functions
#'
#' ## Function used by the targets factories:
#' 1. [load_seurat_data_10X]: load CITE-seq data from 10X cellranger output
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
#' ## Exported functions:
#' - magrittr pipes: \link[magrittr]{%>%}, \link[magrittr]{%<>%}
#'
#' @docType package
#' @name scitargets
"_PACKAGE"

