utils::globalVariables(c(
  "min_RNA", "max_RNA", "mt_cutoff", "cellranger_path",
  "project_name", "r_id", "r_path", "p_id", "min_g", "max_g", "mt",
  "id", "feats", "seurat", "singlets_dim", "cluster_to_use"
))

#' Targets factory to demultiplex CITE-seq
#'
#' This function is used to create the pipeline steps necessary to:
#' - load the CITE-seq data
#' - filter the assay based on the number of detected genes and percentage of mitochondrial counts
#' - demultiplex the HTO
#' - extract non-negative (singlets + doublets) cells
#' - extract singlets cells.
#' - remove non-desired features (gonosomes genes) from singlets
#' - realize a celltype annotation on the singlets using azimuth
#'
#' Note that if we remove features, azimuth will be run in the version without those features.
#'
#' Dimension reduction (PCA, UMAP, t-SNE) are computed after each modification of the Seurat object that modify the count matrix.
#'
#' @param run_id Used to set the names of of the targets. Must be unique across all project
#' @param project_id Used as the project name for the Seurat object.
#' @param run_path Relative or absolute path to the directory that contains another directory whose name correspond to id and contains cellranger output (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz).
#' @param min_genes_detected Minimal number of genes detected in a cell. Default to 200.
#' @param max_genes_detected Maximal number of genes detected in a cell. Default to 4000.
#' @param mt_percent_cutoff Cutoff for the percentage of mitochondrial counts a cell can have.
#' @param singlets_dim_to_use PCA dimensions to use when computing UMAP and t-SNE (only used after extracting the singlets cells).
#' @param singlets_feat_to_remove If not NULL, a character variable indicating a targets step that contains the name of features (genes) to remove from the singlets.
#' @param singlets_clusters_to_use Name of cell cluster to use when identifying markers or running azimuth celltype annotation
#' @param run_azimuth If true, run azimuth using the (human) pbmc-ref of SeuratData
#'
#' @returns A list describing pipeline steps that can be interpreted by targets
#' @export
#'
#' @examples
#' \dontrun{
#' # Check the "get-started" and "targets-factories" vignettes for examples
#' }
tar_demultiplex_hto <- function(
  run_id = "ID_NAME",
  project_id = "PROJECT_ID",
  run_path = NA,
  min_genes_detected = 200,
  max_genes_detected = 4000,
  mt_percent_cutoff = 5,
  singlets_feat_to_remove = NULL,
  singlets_dim_to_use = 1:15,
  singlets_clusters_to_use = NULL,
  run_azimuth = F
) {
  if (is.na(run_path)) {
    stop("Run path must be a valid path pointing to a file inside the 'data/cellranger_output/' folder! 1")
  } else if (fs::is_absolute_path(run_path) &&
    !startsWith(
      normalizePath(run_id, winslash = "/", mustWork = F),
      normalizePath("./data/cellranger_output/", winslash = "/", mustWork = F)
    )
  ) {
    stop("Run path must be a valid path pointing to a file inside the 'data/cellranger_output/' folder! 2")
  } else if (!dir.exists(fs::path("./data/cellranger_output/", run_path))) {
    stop("Run path must be a valid path pointing to a file inside the 'data/cellranger_output/' folder! 3")
  }
  if (!is.character(run_id)) {
    stop("run_id must be a character value corresponding to the date of the sequencing run.")
  }
  if (!is.character(project_id)) {
    stop("project_id must be a character value corresponding to the date of the sequencing run.")
  }
  if (!is.numeric(min_genes_detected) || min_genes_detected <= 0) {
    stop("min_genes_detected must be a positive numeric value.")
  }
  if (!is.numeric(max_genes_detected) || max_genes_detected <= 0) {
    stop("max_genes_detected must be a positive numeric value.")
  }
  if (!is.numeric(mt_percent_cutoff) || mt_percent_cutoff < 0) {
    stop("mt_percent_cutoff_detected a numeric value superior or equal to 0.")
  }
  if (!is.null(singlets_feat_to_remove) && !is.character(singlets_feat_to_remove)) {
    stop("singlets_feat_to_remove must be NULL or a vector of character.")
  }
  if (!is.null(singlets_clusters_to_use) && !is.character(singlets_clusters_to_use)) {
    stop("singlets_clusters_to_use must be NULL or a vector of character representing a column in the metadata of the singlets object created by the pipeline (clusters_0.2, clusters_0.3, ..., clusters_1).")
  }
  hto_demux_steps <- list(
    ##################
    ##  Parameters  ##
    ##################
    targets::tar_target_raw(
      name = base::paste0("parameters_", run_id),
      command = base::substitute(
        expr = list(
          "run" = r_id, "project" = p_id,
          "path" = r_path,
          "min_genes_detected" = min_g,
          "max_genes_detected" = max_g,
          "mt_cutoff" = mt,
          "singlets_clusters_for_markers" = singlets_clusters
        ),
        env = list(
          r_id = run_id, p_id = project_id,
          r_path = run_path,
          min_g = min_genes_detected,
          max_g = max_genes_detected,
          mt = mt_percent_cutoff,
          singlets_clusters = singlets_clusters_to_use
        )
      ),
      description = base::paste0("parameters of the run ", run_id, " from the ", project_id, " project.")
    ),
    ##################
    ##  Files path  ##
    ##################
    targets::tar_target_raw(
      name = base::paste0("cellranger_output_", run_id),
      command = base::substitute(
        base::paste0("./data/cellranger_output/", id, "/"),
        list(id = run_path)
      ),
      format = "file",
      description = base::paste0(
        "Folder containing the genes and antibodies counts generated by cellranger from the ",
        run_id, "_NovaSeqX run (", project_id, ")"
      )
    ),
    ################
    ##  Analysis  ##
    ################
    targets::tar_target_raw(
      name = seurat_obj_target_name(run_id, "raw"),
      command = base::substitute(
        scitargets::load_seurat_data_10X(path = cellranger_path, id = project_name),
        list(
          cellranger_path = base::as.symbol(base::paste0("cellranger_output_", run_id)),
          project_name = project_id
        )
      ),
      description = base::paste0("Seurat object with the counts generated by cellranger ", run_id, " from the ", project_id, " project.")
    ),
    targets::tar_target_raw(
      name = seurat_obj_target_name(run_id, "filter"),
      command = base::substitute(
        scitargets::filter_cell_and_run_reduction(
          obj = seurat, min_nFeature_RNA = min_RNA,
          max_nFeature_RNA = max_RNA, cutoff_percent_mt = mt_cutoff
        ),
        list(
          seurat = base::as.symbol(seurat_obj_target_name(run_id, "raw")),
          min_RNA = min_genes_detected,
          max_RNA = max_genes_detected,
          mt_cutoff = mt_percent_cutoff
        )
      ),
      description = base::paste0(run_id, ": filtered low quality cells based on the number of detected genes and the percentage of mitochondrial genes (", min_genes_detected, " < nFeature_RNA < ", max_genes_detected, " and percent.mt < ", mt_percent_cutoff, ").")
    ),
    targets::tar_target_raw(
      name = seurat_obj_target_name(run_id, "hto_demux"),
      command = base::substitute(
        scitargets::demultiplex_cell(obj = seurat),
        list(seurat = base::as.symbol(seurat_obj_target_name(run_id, "filter")))
      ),
      description = base::paste0(run_id, ": demultiplex HTO.")
    ),
    targets::tar_target_raw(
      name = seurat_obj_target_name(run_id, "non_neg"),
      command = base::substitute(
        scitargets::extract_non_negative_cells(obj = seurat),
        list(seurat = base::as.symbol(seurat_obj_target_name(run_id, "hto_demux")))
      ),
      description = base::paste0(run_id, ": remove negative cells (keep singlets and doublets) and recompute reduction.")
    ),
    targets::tar_target_raw(
      name = seurat_obj_target_name(run_id, "singlets"),
      command = base::substitute(
        scitargets::extract_singlets(obj = seurat, dims_to_use = singlets_dim),
        list(
          seurat = base::as.symbol(seurat_obj_target_name(run_id, "non_neg")),
          singlets_dim = singlets_dim_to_use
        )
      ),
      description = base::paste0(run_id, ": singlets cells (all features of the original object).")
    )
  )

  # Removal non-desired features if they are specified
  if (!is.null(singlets_feat_to_remove)) {
    hto_demux_steps <- c(
      hto_demux_steps,
      targets::tar_target_raw(
        name = seurat_obj_target_name(run_id, "feat_removed_singlets"),
        command = base::substitute(
          scitargets::extract_singlets(
            obj = seurat,
            dims_to_use = singlets_dim,
            feat_to_remove = feats
          ),
          list(
            seurat = base::as.symbol(seurat_obj_target_name(run_id, "singlets")),
            singlets_dim = singlets_dim_to_use,
            feats = base::as.symbol(singlets_feat_to_remove)
          )
        ),
        description = base::paste0(run_id, ": singlets cell without the non desired features.")
      )
    )
  }

  # If clusters names specified:
  # - Identify markers
  # - Add azimuth annotations
  if (!is.null(singlets_clusters_to_use)) {
    if (is.null(singlets_feat_to_remove)) {
      hto_demux_steps <- c(
        hto_demux_steps,
        targets::tar_target_raw(
          name = markers_target_name(run_id),
          command = {
            base::substitute(
              Seurat::FindAllMarkers(
                seurat,
                assay = "SCT",
                group.by = cluster_to_use,
                random.seed = targets::tar_seed_get(),
                only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25
              ),
              list(
                seurat = base::as.symbol(seurat_obj_target_name(run_id, "singlets")),
                cluster_to_use = singlets_clusters_to_use
              )
            )
          },
          description = base::paste0(
            run_id, ": Markers identified for the different cells clusters (",
            singlets_clusters_to_use, ")"
          )
        )
      )
      if (run_azimuth) {
        hto_demux_steps <- c(
          hto_demux_steps,
          name = seurat_obj_target_name(run_id, "azimuth"),
          command = {
            base::substitute(
              scitargets::azimuth_annot_pbmc(seurat, cluster_to_use),
              list(
                seurat = base::as.symbol(seurat_obj_target_name(run_id, "singlets")),
                cluster_to_use = singlets_clusters_to_use
              )
            )
          },
          description = base::paste0(run_id, ": singlets cell with azimuth celltype annotation.")
        )
      }
    } else {
      hto_demux_steps <- c(
        hto_demux_steps,
        targets::tar_target_raw(
          name = markers_target_name(run_id),
          command = {
            base::substitute(
              Seurat::FindAllMarkers(
                object = seurat,
                assay = "SCT",
                group.by = cluster_to_use,
                random.seed = seed,
                only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25
              ),
              list(
                seurat = base::as.symbol(seurat_obj_target_name(run_id, "feat_removed_singlets")),
                cluster_to_use = singlets_clusters_to_use,
                seed = targets::tar_seed_get()
              )
            )
          },
          description = base::paste0(run_id, ": Markers identified for the different cells clusters (", singlets_clusters_to_use, ", after features removal)")
        )
      )
      if (run_azimuth) {
        hto_demux_steps <- c(
          targets::tar_target_raw(
            name = seurat_obj_target_name(run_id, "azimuth"),
            command = {
              base::substitute(
                scitargets::azimuth_annot_pbmc(seurat, cluster_to_use),
                list(
                  seurat = base::as.symbol(seurat_obj_target_name(run_id, "feat_removed_singlets")),
                  cluster_to_use = singlets_clusters_to_use
                )
              )
            },
            description = base::paste0(run_id, ": singlets cell with azimuth celltype annotation (after features removal).")
          )
        )
      }
    }
  }


  return(hto_demux_steps)
}
