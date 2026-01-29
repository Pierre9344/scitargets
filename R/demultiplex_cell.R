utils::globalVariables(c("nFeature_RNA", "percent.mt", "."))

#' Demultiplex cell
#'
#' Demultiplex the cells using the cite-seq experiment HTO.
#'
#' @param obj A Seurat object with the cells HTO informations
#' @param default_seed Default seed used for reduction methods that need it. Not used if the function run inside a targets pipeline as in this case the steps seed (defined based on its name and the pipeline global seed) is used.
#'
#' @returns A seurat object with metadata containing the demultiplexing result.
#' @export
#'
#' @examples
#' \dontrun{
#' demultiplex_cell(path = ".", id = "my_project_id")
#' }
demultiplex_cell <- function(
    obj = NULL,
    default_seed = 1234) {
  if (is.null(obj)) {
    stop("obj must be a Seurat object")
  }
  if (!"HTO" %in% SeuratObject::Assays(obj)) {
    stop("obj must contain a assay named 'HTO' from a cite-seq experiment")
  }
  if (!is.numeric(default_seed)) {
    stop("default_seed must be numeric")
  }
  local_seed <- base::ifelse(targets::tar_active(), targets::tar_seed_get(), default_seed)
  cell <- Seurat::NormalizeData(
    object = obj,
    normalization.method = "CLR",
    margin = 2,
    assay = "HTO"
  )
  tryCatch(
    {
      cell %<>%
        Seurat::HTODemux(
          assay = "HTO",
          positive.quantile = 0.99,
          seed = local_seed
        )
      class_tab <- base::table(cell$HTO_classification.global)
      if (class_tab["Negative"] / class_tab["Singlet"] > 0.05) {
        message("Negatives represent more than 5% of the single. Preparing a second demultiplexing on the negatives and doublets.")
        rd_neg <- base::subset(cell, idents = "Negative") %>%
          Seurat::NormalizeData(
            assay = "HTO",
            normalization.method = "CLR"
          ) %>%
          Seurat::HTODemux(
            assay = "HTO",
            positive.quantile = 0.99,
            seed = local_seed
          )
        negsing_cells <- base::names(base::which(rd_neg$HTO_classification.global == "Singlet"))
        cell@meta.data[negsing_cells, ] <- rd_neg@meta.data[negsing_cells, ]
      }
      if (class_tab["Doublet"] / class_tab["Singlet"] > 0.05) {
        message("Doublets represent more than 5% of the single. Preparing a second demultiplexing on the negatives and doublets.")
        rd_dup <- base::subset(cell, idents = "Doublet") %>%
          Seurat::NormalizeData(
            assay = "HTO",
            normalization.method = "CLR"
          ) %>%
          Seurat::HTODemux(
            assay = "HTO",
            positive.quantile = 0.99,
            seed = local_seed
          )
        dup_cells <- base::names(base::which(rd_dup$HTO_classification.global == "Singlet"))
        cell@meta.data[dup_cells, ] <- rd_dup@meta.data[dup_cells, ]
      }
    },
    error = function(e) {
      message(e)
      message("First demultiplexing try failed. Second try with HTO counts per cell > 10% median.")
      tryCatch({
        # Second try with HTO counts > 10% median
        cell_sub <- names(which(colSums(cell@assays$HTO@counts) < stats::median(colSums(cell@assays$HTO@counts)) / 10))
        cell %<>%
          subset(cell = cell_sub, invert = TRUE) %>%
          Seurat::HTODemux(assay = "HTO", positive.quantile = 0.99)
        message(paste("Number of cells Excluded : ", length(cell_sub)))
        tryCatch(
          {
            # Second pass HTOdemux on negative cells if negative or doublet > 5% of singlet
            class_tab <- table(cell$HTO_classification.global)
            if (class_tab["Negative"] / class_tab["Singlet"] > 0.05) {
              rd_neg <- cell %>%
                subset(idents = "Negative") %>%
                Seurat::NormalizeData(assay = "HTO", normalization.method = "CLR") %>%
                Seurat::HTODemux(assay = "HTO", positive.quantile = 0.99)
              negsing_cells <- names(which(rd_neg$HTO_classification.global == "Singlet"))
              cell@meta.data[negsing_cells, ] <- rd_neg@meta.data[negsing_cells, ]
            }
            if (class_tab["Doublet"] / class_tab["Singlet"] > 0.05) {
              rd_dup <- cell %>%
                subset(idents = "Doublet") %>%
                Seurat::NormalizeData(assay = "HTO", normalization.method = "CLR") %>%
                Seurat::HTODemux(assay = "HTO", positive.quantile = 0.99)
              dup_cells <- names(which(rd_dup$HTO_classification.global == "Singlet"))
              cell@meta.data[dup_cells, ] <- rd_dup@meta.data[dup_cells, ]
            }
          },
          error = function(e) {
            message(e)
            message("Second demultiplexing try failed. Last try with MULTIseqDemux.")
            cell %<>% Seurat::MULTIseqDemux(cell, assay = "HTO", quantile = 0.99)
            cell$HTO_classification.global <- cell$MULTI_ID
            cell$HTO_classification.global[!cell$HTO_classification.global %in% c("Negative", "Doublet")] <- NA
            cell$HTO_classification.global <- droplevels(cell$HTO_classification.global)
            levels(cell$HTO_classification.global) <- c(levels(cell$HTO_classification.global), "Singlet")
            cell$HTO_classification.global[is.na(cell$HTO_classification.global)] <- "Singlet"
            cell$HTO_classification <- cell$MULTI_classification
            class_tab <- table(cell$HTO_classification.global)
            if (class_tab["Negative"] / class_tab["Singlet"] > 0.05) {
              print("Double Demux Negative")
              rd_neg <- cell %>%
                subset(idents = "Negative") %>%
                Seurat::NormalizeData(assay = "HTO", normalization.method = "CLR") %>%
                Seurat::MULTIseqDemux(assay = "HTO", quantile = 0.99)
              rd_neg$HTO_classification.global <- rd_neg$MULTI_ID
              rd_neg$HTO_classification.global[!rd_neg$HTO_classification.global %in% c("Negative", "Doublet")] <- NA
              rd_neg$HTO_classification.global <- droplevels(rd_neg$HTO_classification.global)
              levels(rd_neg$HTO_classification.global) <- c(levels(rd_neg$HTO_classification.global), "Singlet")
              rd_neg$HTO_classification.global[is.na(rd_neg$HTO_classification.global)] <- "Singlet"
              rd_neg$HTO_classification <- rd_neg$MULTI_classification
              negsing_cells <- colnames(rd_neg)[which(rd_neg$HTO_classification.global == "Singlet")]
              cell@meta.data[negsing_cells, ] <- rd_neg@meta.data[negsing_cells, ]
            }
            if (class_tab["Doublet"] / class_tab["Singlet"] > 0.05) {
              print("Double Demux Doublet")
              rd_dup <- cell %>%
                subset(idents = "Doublet") %>%
                Seurat::NormalizeData(assay = "HTO", normalization.method = "CLR") %>%
                Seurat::MULTIseqDemux(assay = "HTO", positive.quantile = 0.99)
              rd_dup$HTO_classification.global <- rd_dup$MULTI_ID
              rd_dup$HTO_classification.global[!rd_dup$HTO_classification.global %in% c("Negative", "Doublet")] <- NA
              rd_dup$HTO_classification.global <- droplevels(rd_dup$HTO_classification.global)
              levels(rd_dup$HTO_classification.global) <- c(levels(rd_dup$HTO_classification.global), "Singlet")
              rd_dup$HTO_classification.global[is.na(rd_dup$HTO_classification.global)] <- "Singlet"
              rd_dup$HTO_classification <- rd_dup$MULTI_classification
              dup_cells <- names(which(rd_dup$HTO_classification.global == "Singlet"))
              cell@meta.data[dup_cells, ] <- rd_dup@meta.data[dup_cells, ]
            }
            cell$hash.ID <- cell$MULTI_ID
            cell$MULTI_ID <- NULL
          })
      })
    })
  return(cell)
}
