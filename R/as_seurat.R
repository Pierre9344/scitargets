#' Coerce a Seurat object, or a path to one, into a Seurat object
#'
#' The target factories in this package are given the NAME of an upstream target
#' holding a Seurat object. That target does not always hold the object itself.
#' A pipeline whose merged object is too expensive to rebuild everywhere can
#' write it to disk and declare its target `format = "file"`, so that the
#' expensive normalisation runs once on a machine with enough memory and every
#' other machine reads the result. The target's VALUE is then a path, not an
#' object.
#'
#' `as_seurat()` accepts either form, so a factory no longer has to care which
#' one it was handed, and a pipeline can switch a target between the two without
#' the downstream steps changing.
#'
#' The reader is chosen from the file extension:
#'
#' | extension | reader |
#' | --- | --- |
#' | `.qs2`, `.qs` | [qs2::qs_read()] |
#' | `.rds` | [base::readRDS()] |
#' | `.RData`, `.rda` | [base::load()] |
#'
#' `.RData` must hold exactly ONE object. `save()`/`load()` carry the variable's
#' name inside the file rather than returning the value, so a file holding
#' several objects has no unambiguous answer to "which one is the Seurat
#' object?". `.qs2` and `.rds` do not have that problem and are preferred.
#'
#' No caching happens here: the file is read on every call. A pipeline whose
#' branches each need the object should memoise on its own side rather than pay
#' the read per branch.
#'
#' @param x A Seurat object, or a single path to a file holding one.
#' @param nthreads Decompression threads, passed to [qs2::qs_read()]. Ignored by
#'   the other readers.
#' @returns A Seurat object.
#' @examples
#' \dontrun{
#' # Both of these give the object:
#' as_seurat(seurat_obj)
#' as_seurat("out/seurat/merged_filtered.qs2")
#' }
#' @export
as_seurat <- function(x, nthreads = 1L) {
  if (methods::is(x, "Seurat")) {
    return(x)
  }
  if (!is.character(x) || length(x) != 1L) {
    stop("as_seurat(): expected a Seurat object or a single file path, got ",
      base::class(x)[1L],
      if (is.character(x)) base::paste0(" of length ", length(x)) else "",
      ".", call. = FALSE)
  }
  if (is.na(x) || !nzchar(x)) {
    stop("as_seurat(): the file path is ",
      if (is.na(x)) "NA" else "an empty string",
      ". The upstream target probably did not produce one.", call. = FALSE)
  }
  if (!file.exists(x)) {
    stop("as_seurat(): no file at '", x, "'.", call. = FALSE)
  }
  ext <- base::tolower(base::sub(".*\\.", "", base::basename(x)))
  obj <- switch(ext,
    qs2 = ,
    qs = {
      if (!requireNamespace("qs2", quietly = TRUE)) {
        stop("as_seurat(): reading '", base::basename(x), "' needs the qs2 ",
          "package. install.packages(\"qs2\")", call. = FALSE)
      }
      qs2::qs_read(x, nthreads = nthreads)
    },
    rds = base::readRDS(x),
    rdata = ,
    rda = {
      e <- base::new.env(parent = base::emptyenv())
      nms <- base::load(x, envir = e)
      if (length(nms) != 1L) {
        stop("as_seurat(): '", x, "' holds ", length(nms), " objects (",
          base::paste(nms, collapse = ", "), "); it must hold exactly one. ",
          "Prefer .qs2 or .rds, which store a value rather than a name.",
          call. = FALSE)
      }
      base::get(nms, envir = e)
    },
    stop("as_seurat(): do not know how to read '", base::basename(x),
      "'. Supported extensions: .qs2, .qs, .rds, .RData, .rda.", call. = FALSE)
  )
  if (!methods::is(obj, "Seurat")) {
    stop("as_seurat(): '", x, "' holds an object of class ",
      base::class(obj)[1L], ", not a Seurat object.", call. = FALSE)
  }
  obj
}

#' File extensions [as_seurat()] and [save_seurat()] understand
#'
#' Exposed so a factory can reject a badly named path when the pipeline is
#' DEFINED rather than after the work that was supposed to fill it.
#'
#' @returns A character vector of lower-case extensions, without the dot.
#' @export
seurat_file_formats <- function() {
  c("qs2", "qs", "rds", "rdata", "rda")
}

#' Write a Seurat object, choosing the writer from the file extension
#'
#' The counterpart of [as_seurat()]: the two agree on which extensions mean
#' which format, so a file written here reads back there.
#'
#' | extension | writer |
#' | --- | --- |
#' | `.qs2`, `.qs` | [qs2::qs_save()] |
#' | `.rds` | [base::saveRDS()] |
#' | `.RData`, `.rda` | [base::save()], under the name `obj` |
#'
#' Prefer `.qs2`: it is markedly smaller and faster than `.rds` for a Seurat
#' object, and unlike `.RData` it stores a value rather than a variable name.
#'
#' The parent directory is created when missing, so a caller does not have to.
#'
#' @param obj A Seurat object.
#' @param path Destination file. Its extension picks the writer.
#' @param nthreads Compression threads, passed to [qs2::qs_save()]. Ignored by
#'   the other writers.
#' @returns `path`, invisibly, so the call can end a `format = "file"` target.
#' @export
save_seurat <- function(obj, path, nthreads = 1L) {
  if (!methods::is(obj, "Seurat")) {
    stop("save_seurat(): `obj` must be a Seurat object, got ",
      base::class(obj)[1L], ".", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("save_seurat(): `path` must be a single non-empty file path.",
      call. = FALSE)
  }
  ext <- base::tolower(base::sub(".*\\.", "", base::basename(path)))
  if (!ext %in% seurat_file_formats()) {
    stop("save_seurat(): do not know how to write '", base::basename(path),
      "'. Supported extensions: .qs2, .qs, .rds, .RData, .rda.", call. = FALSE)
  }
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)
  switch(ext,
    qs2 = ,
    qs = {
      if (!requireNamespace("qs2", quietly = TRUE)) {
        stop("save_seurat(): writing '", base::basename(path), "' needs the ",
          "qs2 package. install.packages(\"qs2\")", call. = FALSE)
      }
      qs2::qs_save(obj, path, nthreads = nthreads)
    },
    rds = base::saveRDS(obj, path),
    rdata = ,
    # Exactly one object, named `obj`, so as_seurat() finds a single name to
    # restore. save() stores the NAME rather than the value, which is why a
    # multi-object .RData has no unambiguous answer on the way back in.
    rda = base::save(obj, file = path)
  )
  invisible(path)
}
