#' Extract the non-negative cells (singlets and doublets)
#'
#' Filter the Seurat object to keep only the singlets and the doublets. Recompute UMAP and t-SNE before returning the object.
#'
#' @param obj a Seurat object.
#' @param hto_ident The name of the column in obj metadata that contains the HTO demultiplexing results. Default to `HTO_classification.global`.
#' @param default_seed Default seed used for reduction methods that need it. Not used if the function run inside a targets pipeline as in this case the steps seed (defined based on its name and the pipeline global seed) is used.
#'
#' @returns A seurat object filtered to contains only singlets and doublets
#' @export
#'
#' @examples
#' \dontrun{
#' extract_non_negative_cells(obj = Seurat_obj)
#' }
extract_non_negative_cells <- function(obj = NULL,
                                   hto_ident = "HTO_classification.global",
                                   default_seed = 1234) {
  if (is.null(obj)) {
    stop("obj must be a Seurat object")
  }
  if (! hto_ident %in% colnames(obj[[]])) {
    stop("hto_ident must be part of obj metadata")
  }
  local_seed <- base::ifelse(targets::tar_active(), targets::tar_seed_get(), default_seed)
  SeuratObject::DefaultAssay(obj) <- "RNA"
  SeuratObject::Idents(obj) <- hto_ident
  obj %<>%
    base::subset(idents = "Negative", invert = TRUE) %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures(selection.method = "mean.var.plot") %>%
    Seurat::ScaleData(features = SeuratObject::VariableFeatures(.)) %>%
    Seurat::RunPCA(seed.use = local_seed)
  elbow <- Seurat::ElbowPlot(obj)
  max_dim <- utils::tail(elbow$data$dims[elbow$data$stdev > stats::median(elbow$data$stdev)], 1)
  message(paste0("Computing UMAP and t-SNE on non-negative cells (singlets and doublets). Determining the number of PC to use (1:", max_dim, ") by their standard deviation."))
  obj %<>%
    Seurat::FindNeighbors(
      reduction = "pca",
      dims = 1:max_dim
    ) %>%
    Seurat::FindClusters(
      resolution = 0.5,
      algorithm = 4,
      random.seed = abs(local_seed)
    ) %>%
    Seurat::RunUMAP(
      dims = 1:max_dim,
      seed.use = local_seed
    ) %>%
    Seurat::RunTSNE(
      dims = 1:max_dim,
      perplexity = 100,
      check_duplicates = FALSE,
      seed.use = local_seed
    )
  obj$RNA_snn_res.0.5 <- NULL
  return(obj)
}
