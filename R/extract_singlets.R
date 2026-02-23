#' Extract the singlets
#'
#' Filter the Seurat object to keep only the singlets identified using the HTO and removes features (genes) if necessary before recomputing the UMAP and t-SNE and identify cluster for different resolution (0.2 to 1).
#'
#' @param obj obj a Seurat object.
#' @param dims_to_use Number of principal components to use when computing the UMAP and t-SNE. Default to 1:15. If a NULL or non-integer value is set, the number will be determined by the stdev of the PC.
#' @param feat_to_remove Names of features (genes) to remove. Can be used to remove genes located on gonosomes.
#' @param default_seed default_seed Default seed used for reduction methods that need it. Not used if the function run inside a targets pipeline as in this case the steps seed (defined based on its name and the pipeline global seed) is used.
#'
#' @returns A Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' extract_singletss(obj = Seurat_obj)
#' }
extract_singlets <- function(obj = NULL,
                             dims_to_use = 1:15,
                             feat_to_remove = NULL,
                             default_seed = 1234) {
  if (is.null(obj)) {
    stop("obj must be a Seurat object")
  } else if (!is.integer(dims_to_use)) {
    message("dims_to_use is not set, it will be replaced using the PCA elbow plot")
  }
  local_seed <- base::ifelse(targets::tar_active(), targets::tar_seed_get(), default_seed)
  # Subset singlet, set default assay and idents
  obj %<>% base::subset(
    .$HTO_classification.global == "Singlet"
    )
  # Remove chr Y genes if feat_to_remove is valid.
  # UMAP were they are not removed show a difference between man and woman (see the report)
  if (!is.null(feat_to_remove) && is.character(feat_to_remove)) {
    to_keep <- dplyr::setdiff(unique(rownames(obj)), feat_to_remove)
    SeuratObject::DefaultAssay(obj) <- "HTO"
    obj[["RNA"]] <- base::subset(obj[["RNA"]], features = to_keep)
  }
  SeuratObject::DefaultAssay(obj) <- "RNA"
  SeuratObject::Idents(obj) <- "HTO_classification"
  # Remove metadata that are no longer needed
  obj@meta.data %<>%
    dplyr::select(-tidyselect::any_of(c("RNA_snn_res.0.5",
                                        "SCT_snn_res.0.5",
                                        "seurat_clusters")))
  # Recompute dimensions reductions and clusters
  obj %<>%
    Seurat::SCTransform(
      vars.to.regress = "percent.mt",
      verbose = FALSE,
      vst.flavor = "v2",
      seed.use = local_seed
    ) %>%
    Seurat::RunPCA(
      seed.use = local_seed
    ) %>%
    Seurat::RunUMAP(
      dims = dims_to_use,
      seed.use = local_seed
    ) %>%
    Seurat::FindNeighbors(
      reduction = "pca",
      dims = dims_to_use
    )
  # identify clusters with different resolutions
  for (k in base::seq(from = 0.2, to = 1, by = 0.1)) {
    obj %<>%
      Seurat::FindClusters(
        resolution = k,
        algorithm = 4,
        cluster.name = paste0("clusters_", k),
        random.seed = abs(local_seed)
      )
  }
  return(obj)
}
