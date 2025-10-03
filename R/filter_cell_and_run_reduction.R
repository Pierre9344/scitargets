utils::globalVariables(c("nFeature_RNA", "percent.mt", "."))

#' Filter cell and run reduction
#'
#' Filter cells based on the number of genes they contains and the percentage of count represented by mitochondrial genes.
#' After the filtering, PCA, UMAP, and t-SNE are computed. The UMAP and t-SNE are based on the PCA (PC used are those whose standard deviation is superior to the median of all PC standard deviation).
#'
#' @param obj a Seurat object.
#' @param min_nFeature_RNA Minimum number of genes in the filtered cells.
#' @param max_nFeature_RNA Maximum number of genes in the filtered cells.
#' @param cutoff_percent_mt Percentage of count corresponding to mitochondrial genes used to filter the cells. Cells with that percentage or more will be removed.
#' @param default_seed Default seed used for reduction methods that need it. Not used if the function run inside a targets pipeline as in this case the steps seed (defined based on its name and the pipeline global seed) is used.
#'
#' @returns A seurat object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' filter_cell_and_run_reduction(obj = Seurat_obj)
#' }
filter_cell_and_run_reduction <- function(obj = NULL,
                                          min_nFeature_RNA = 200,
                                          max_nFeature_RNA = 4000,
                                          cutoff_percent_mt = 5,
                                          default_seed = 1234) {
  if (is.null(obj)) {
    stop("obj must be a Seurat object")
  }
  local_seed <- base::ifelse(targets::tar_active(), targets::tar_seed_get(), default_seed)
  obj %<>%
    base::subset(
      subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < cutoff_percent_mt
    ) %>%
    Seurat::SCTransform(
      vars.to.regress = "percent.mt",
      verbose = FALSE,
      vst.flavor = "v2",
      seed = local_seed
    ) %>%
    Seurat::RunPCA(
      features = SeuratObject::VariableFeatures(.),
      seed.use = local_seed
    )
  elbow <- Seurat::ElbowPlot(obj)
  max_dim <- utils::tail(elbow$data$dims[elbow$data$stdev > stats::median(elbow$data$stdev)], 1)
  message(paste0("Computing UMAP and t-SNE. Determining the number of PC to use (1:", max_dim, ") by their standard deviation."))
  obj %<>%
    Seurat::FindNeighbors(
      reduction = "pca",
      dims = 1:max_dim
    ) %>%
    Seurat::FindClusters(
      resolution = 0.5,
      algorithm = 4,
      random.seed = base::abs(local_seed)
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
  # Not necessary to keep the clusters at this step
  obj$SCT_snn_res.0.5 <- NULL
  obj$seurat_clusters <- NULL
  return(obj)
}
