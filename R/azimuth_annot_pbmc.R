#' Azimuth annotations for PBMC
#'
#' @param obj a Seurat object with an SCT assay and a metadata field corresponding to cell clusters
#' @param cluster_to_use Name of the cells clusters
#' @param reference Reference passed to [Azimuth::RunAzimuth()]. `NULL` (the
#'   default) uses the `pbmcref.SeuratData` package's own `azimuth/` directory
#'   when that package is installed, and falls back to the name `"pbmcref"`
#'   otherwise. Pass a directory holding `ref.Rds` and `idx.annoy` to use a
#'   reference from elsewhere.
#' @param homolog_table Path to a local copy of Azimuth's `homologs.rds`. `NULL`
#'   (the default) reads `getOption("scitargets.azimuth_homologs")`, then the
#'   `AZIMUTH_HOMOLOGS` environment variable, and leaves Azimuth to download the
#'   table when neither is set. Fetch a copy from
#'   `https://seurat.nygenome.org/azimuth/references/homologs.rds` on a machine
#'   with access.
#'
#' @details
#' Two separate downloads have to be dealt with, and only the first is
#' configurable.
#'
#' `Azimuth:::LoadReference()` treats its `path` as a URL unless it is an
#' existing local directory, so `reference = "pbmcref"` always downloads from
#' seurat.nygenome.org even when the reference is already installed. Resolving
#' the installed package's directory first takes the local branch instead.
#'
#' `RunAzimuth.Seurat()` then calls `ConvertGeneNames()` with the homolog table
#' URL **hardcoded**, with no argument to override it. `ConvertGeneNames()` does
#' accept a local file, so when `homolog_table` is given this function swaps that
#' function in Azimuth's namespace for the duration of the call, restoring it
#' immediately afterwards. That is a workaround for a hardcoded URL rather than a
#' supported extension point, and it is skipped entirely unless a table is
#' configured, so the default behaviour is untouched.
#'
#' On a host that can reach seurat.nygenome.org neither applies. On one that
#' cannot (an HPC compute node behind a proxy answering 403) both are needed:
#' fixing only the reference moves the failure from the reference download to the
#' homolog download.
#'
#' @returns A seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' azimuth_annot_pbmc(SeuratObj, "clusters_name_in_metadata")
#' }
azimuth_annot_pbmc <- function(obj, cluster_to_use = "clusters_0.5",
                               reference = NULL, homolog_table = NULL) {
  if (is.null(obj)) {
    stop("obj must be a Seurat object")
  }
  if (!is.character(cluster_to_use)) {
    stop("cluster_to_use must br a character variable present in the meta.data of obj")
  } else if (cluster_to_use %in% obj@meta.data) {
    stop("cluster_to_use must br a character variable present in the meta.data of obj")
  }
  if (is.null(reference)) {
    reference <- local_pbmcref() %||% "pbmcref"
  }
  if (is.null(homolog_table)) {
    homolog_table <- local_homolog_table()
  }
  SeuratObject::Idents(obj) <- cluster_to_use
  PreviousDefaultAssay <- SeuratObject::DefaultAssay(obj)
  SeuratObject::DefaultAssay(obj) <- "RNA"
  SeuratObject::DefaultAssay(obj) <- PreviousDefaultAssay
  with_local_homologs(
    homolog_table,
    Azimuth::RunAzimuth(obj, reference = reference, assay = "RNA")
  )
}

#' Path to a local Azimuth homolog table
#'
#' Resolved from `getOption("scitargets.azimuth_homologs")`, then the
#' `AZIMUTH_HOMOLOGS` environment variable. Returns `NULL` when neither is set,
#' or when the configured path does not exist, in which case Azimuth downloads
#' the table as usual.
#'
#' @returns A file path, or `NULL`.
#' @export
local_homolog_table <- function() {
  path <- getOption("scitargets.azimuth_homologs", default = NULL)
  if (is.null(path) || !nzchar(path)) {
    path <- Sys.getenv("AZIMUTH_HOMOLOGS", unset = "")
  }
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  path
}

#' Run an expression with Azimuth's homolog table taken from a local file
#'
#' `RunAzimuth.Seurat()` hardcodes the homolog table URL, so the only way to
#' point it at a local copy is to replace `ConvertGeneNames()` while it runs.
#' The original is restored on exit, including on error.
#'
#' @param path Local homolog table, or `NULL` to run `expr` unchanged.
#' @param expr Expression to evaluate.
#' @returns The value of `expr`.
#' @keywords internal
with_local_homologs <- function(path, expr) {
  if (is.null(path)) {
    return(expr)
  }
  ns <- asNamespace("Azimuth")
  original <- get("ConvertGeneNames", envir = ns)
  patched <- function(object, reference.names, homolog.table) {
    # homolog.table is deliberately ignored: it is the hardcoded URL.
    original(object = object, reference.names = reference.names,
             homolog.table = path)
  }
  utils::assignInNamespace("ConvertGeneNames", patched, ns = "Azimuth")
  on.exit(
    utils::assignInNamespace("ConvertGeneNames", original, ns = "Azimuth"),
    add = TRUE
  )
  # `expr` is a promise, so it is evaluated HERE, with the patch installed, and
  # not when the argument was matched. force() makes that explicit rather than
  # leaving it to a bare symbol on the last line.
  force(expr)
}

#' Path to an installed Azimuth PBMC reference
#'
#' The directory inside `pbmcref.SeuratData` holding the two files
#' [Azimuth::RunAzimuth()] needs, or `NULL` when the package is not installed or
#' does not carry both files. Exported so a pipeline can check the reference is
#' usable before spending an hour reaching the annotation step.
#'
#' @returns A directory path, or `NULL`.
#' @export
local_pbmcref <- function() {
  if (!requireNamespace("pbmcref.SeuratData", quietly = TRUE)) {
    return(NULL)
  }
  dir <- system.file("azimuth", package = "pbmcref.SeuratData")
  # Both files, not just the directory: LoadReference() checks for exactly these
  # and a partial install would fail later with a less obvious message.
  if (!nzchar(dir) || !all(file.exists(file.path(dir, c("ref.Rds", "idx.annoy"))))) {
    return(NULL)
  }
  dir
}
