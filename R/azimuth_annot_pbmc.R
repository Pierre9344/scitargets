#' Azimuth annotations for PBMC
#'
#' @param obj a Seurat object with an SCT assay and a metadata field corresponding to cell clusters
#' @param cluster_to_use Name of the cells clusters
#'
#' @returns A seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' azimuth_annot_pbmc(SeuratObj, "clusters_name_in_metadata")
#' }
azimuth_annot_pbmc <- function(obj, cluster_to_use = "clusters_0.5") {
  if (is.null(obj)) {
    stop("obj must be a Seurat object")
  }
  if (!is.character(cluster_to_use)) {
    stop("cluster_to_use must br a character variable present in the meta.data of obj")
  } else if (cluster_to_use %in% obj@meta.data) {
    stop("cluster_to_use must br a character variable present in the meta.data of obj")
  }
  SeuratObject::Idents(obj) <- cluster_to_use
  SeuratObject::DefaultAssay(obj) <- "SCT"
  return(Azimuth::RunAzimuth(obj, reference = "pbmcref", assay = "SCT"))
}
