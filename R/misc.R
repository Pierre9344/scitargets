seurat_obj_target_name <- function(id, suffix) {
  sprintf("seurat_obj_%s_%s", id, suffix)
}

markers_target_name <- function(id) {
  sprintf("markers_%s", id)
}
