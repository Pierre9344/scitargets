# hdWGCNA report helpers -------------------------------------------------------
#
# Plotting/table helpers for the per-cell-type hdWGCNA report, plus
# `hdwgcna_report_lines()` which assembles a full Quarto section for one cell
# group (mirrors `dea_report_lines()`). The helpers hold the plotting logic so
# the generated chunks stay short; `hdwgcna_report_lines()` only builds text.

.hdwgcna_require <- function() {
  if (!requireNamespace("hdWGCNA", quietly = TRUE)) {
    stop("Package 'hdWGCNA' is required for the hdWGCNA report helpers.", call. = FALSE)
  }
}

#' Modules associated with at least one trait
#'
#' Returns the module names whose harmonised module eigengene is significantly
#' correlated (FDR < `p_cutoff`) with at least one trait in the object's
#' `ModuleTraitCorrelation` (`"all_cells"` slot). The `"grey"` module is dropped.
#'
#' @param obj An hdWGCNA Seurat object carrying a `ModuleTraitCorrelation` result.
#' @param p_cutoff FDR threshold (default 0.05).
#' @returns A character vector of module names (possibly empty).
#' @export
hdwgcna_associated_modules <- function(obj, p_cutoff = 0.05) {
  .hdwgcna_require()
  mt <- hdWGCNA::GetModuleTraitCorrelation(obj)
  fdr_m <- as.matrix(mt$fdr[["all_cells"]])
  mods <- colnames(fdr_m)[apply(fdr_m, 2, function(cl) any(cl < p_cutoff, na.rm = TRUE))]
  setdiff(mods, "grey")
}

#' Long-format module composition table
#'
#' Returns the gene-to-module assignment from [hdWGCNA::GetModules()] in a tidy
#' long format: one row per gene, with its module, the module colour and the
#' gene's intramodular connectivity (kME) within its assigned module. Rows are
#' ordered by module and then by decreasing kME, so each module's hub genes come
#' first. This is the table [tar_hdwgcna()] writes to an `.xlsx` file (one per
#' cell group) when its `module_table_file` argument is set.
#'
#' @param obj An hdWGCNA Seurat object (after `ConstructNetwork` /
#'   `ModuleEigengenes`), as produced by the `wgcna_<group>` target.
#' @param include_grey Keep the unassigned `"grey"` genes (default `TRUE`).
#' @returns A `data.frame` with columns `module`, `color`, `gene` and `kME`.
#' @export
hdwgcna_module_composition <- function(obj, include_grey = TRUE) {
  .hdwgcna_require()
  mods <- hdWGCNA::GetModules(obj)
  mod_chr <- as.character(mods$module)
  # kME of each gene within its own module (the wide kME_<module> columns held
  # by GetModules() are collapsed to a single per-gene membership value).
  kme <- vapply(seq_len(nrow(mods)), function(i) {
    col <- paste0("kME_", mod_chr[i])
    if (col %in% names(mods)) as.numeric(mods[[col]][i]) else NA_real_
  }, numeric(1))
  out <- data.frame(
    module = mod_chr,
    color  = as.character(mods$color),
    gene   = as.character(mods$gene_name),
    kME    = kme,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(include_grey)) out <- out[out$module != "grey", , drop = FALSE]
  out <- out[order(out$module, -out$kME), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Module-to-gene list for one cell group
#'
#' Returns the gene membership of every module of one hdWGCNA network object as a
#' named `list`: one element per module, holding the character vector of genes
#' assigned to that module (ordered by decreasing kME, so hub genes come first).
#' This is the per-cell-group building block of the nested module gene list
#' produced by [tar_hdwgcna()] (`wgcna_module_gene_list`), where the elements are
#' nested under each cell type, i.e. `ls[[<cell type>]][[<module>]]`.
#'
#' @param obj An hdWGCNA Seurat object (after `ConstructNetwork` /
#'   `ModuleEigengenes`), as produced by the `wgcna_<group>` target.
#' @param include_grey Keep the unassigned `"grey"` module (default `FALSE`).
#' @returns A named `list` of character vectors (module name -> genes). Empty
#'   modules are dropped; module order follows the object's module levels.
#' @export
hdwgcna_module_gene_list <- function(obj, include_grey = FALSE) {
  comp <- hdwgcna_module_composition(obj, include_grey = include_grey)
  # preserve the module factor ordering from GetModules()
  mods <- hdWGCNA::GetModules(obj)
  lvls <- levels(mods$module)
  if (!isTRUE(include_grey)) lvls <- setdiff(lvls, "grey")
  lvls <- lvls[lvls %in% comp$module]
  out <- lapply(lvls, function(m) comp$gene[comp$module == m])
  names(out) <- lvls
  out[vapply(out, length, integer(1)) > 0L]
}

# Build the traits x modules correlation / FDR matrices from the "all_cells"
# slot, dropping grey. `trait_labels` (named: one-hot column -> display label)
# is applied to the rows; NULL keeps the raw row names.
.hdwgcna_mt_long <- function(obj, trait_labels = NULL) {
  .hdwgcna_require()
  mt <- hdWGCNA::GetModuleTraitCorrelation(obj)
  cor_m <- as.matrix(mt$cor[["all_cells"]])
  fdr_m <- as.matrix(mt$fdr[["all_cells"]])[rownames(cor_m), colnames(cor_m), drop = FALSE]
  keep <- colnames(cor_m) != "grey"
  cor_m <- cor_m[, keep, drop = FALSE]
  fdr_m <- fdr_m[, keep, drop = FALSE]
  rn <- rownames(cor_m)
  labs_v <- if (is.null(trait_labels)) stats::setNames(rn, rn) else trait_labels
  list(
    cor = cor_m, fdr = fdr_m, rn = rn,
    label = unname(labs_v[rn])
  )
}

#' Module-trait correlation heatmap (FDR in the tiles)
#'
#' A ggplot tile heatmap of the `"all_cells"` module-trait correlations: tile
#' colour is the correlation `r`, tile text is the FDR with significance stars.
#' Modules are columns, traits (clinical groups) are rows.
#'
#' @param obj An hdWGCNA Seurat object with a `ModuleTraitCorrelation` result.
#' @param trait_labels Optional named character vector mapping the one-hot trait
#'   columns (e.g. `is_CR1M6M`) to display labels (e.g. `CR1M/6M`).
#' @param p_cutoff FDR threshold for the single significance star (default 0.05).
#' @returns A `ggplot` object.
#' @export
hdwgcna_module_trait_heatmap <- function(obj, trait_labels = NULL, p_cutoff = 0.05) {
  m <- .hdwgcna_mt_long(obj, trait_labels)
  cor_m <- m$cor
  fdr_m <- m$fdr
  mt_df <- data.frame(
    clinical_group = factor(rep(m$label, times = ncol(cor_m)), levels = rev(m$label)),
    module = factor(rep(colnames(cor_m), each = nrow(cor_m)), levels = colnames(cor_m)),
    r = as.vector(cor_m),
    fdr = as.vector(fdr_m),
    stringsAsFactors = FALSE
  )
  mt_df$stars <- ifelse(is.na(mt_df$fdr), "",
    ifelse(mt_df$fdr < 0.001, "***",
      ifelse(mt_df$fdr < 0.01, "**",
        ifelse(mt_df$fdr < p_cutoff, "*", ""))))
  mt_df$label <- ifelse(is.na(mt_df$fdr), "",
    paste0(formatC(mt_df$fdr, format = "g", digits = 2),
      ifelse(mt_df$stars == "", "", paste0("\n", mt_df$stars))))
  lim <- max(abs(mt_df$r), na.rm = TRUE)
  ggplot2::ggplot(mt_df, ggplot2::aes(x = module, y = clinical_group, fill = r)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3, lineheight = 0.85) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "grey95", high = "red",
      midpoint = 0, limits = c(-lim, lim), name = "Correlation (r)"
    ) +
    ggplot2::labs(x = NULL, y = "Clinical group") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

#' Module-trait correlation table
#'
#' A long data frame of the `"all_cells"` module-trait correlations (`module`,
#' `clinical_group`, `r`, `fdr`), sorted by FDR. Companion to
#' [hdwgcna_module_trait_heatmap()].
#'
#' @inheritParams hdwgcna_module_trait_heatmap
#' @returns A `data.frame`.
#' @export
hdwgcna_module_trait_table <- function(obj, trait_labels = NULL) {
  m <- .hdwgcna_mt_long(obj, trait_labels)
  cor_m <- m$cor
  fdr_m <- m$fdr
  tab <- data.frame(
    module = rep(colnames(cor_m), each = nrow(cor_m)),
    clinical_group = rep(m$label, times = ncol(cor_m)),
    r = round(as.vector(cor_m), 3),
    fdr = signif(as.vector(fdr_m), 3),
    row.names = NULL, stringsAsFactors = FALSE
  )
  tab[order(tab$fdr), ]
}

# tidytext-style within-facet ordering, inlined to avoid the dependency.
.hdwgcna_wrap_terms <- function(x, width = 45) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}
.hdwgcna_reorder_within <- function(x, by, within, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by, FUN = mean)
}
.hdwgcna_scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

#' Enrichr GO barplot for the hdWGCNA modules
#'
#' For one Enrichr database, draws the top enriched GO terms per module as a
#' faceted barplot (bar length = `-log10(adjusted p-value)`, dashed line at FDR
#' 0.05). Mirrors the GO barplot of [dea_report_lines()].
#'
#' @param enrichr_df The combined Enrichr table (a `wgcna_<group>_enrichr`
#'   target), with columns `db`, `module`, `Term`, `Adjusted.P.value`, etc.
#' @param db The database to plot (one value of `enrichr_df$db`).
#' @param modules Optional character vector restricting the modules shown (e.g.
#'   the clinically-associated ones from [hdwgcna_associated_modules()]); `NULL`
#'   shows all non-grey modules.
#' @param n_terms Top terms per module (default 6).
#' @param interactive If `TRUE` (default) return a `ggiraph` girafe with per-bar
#'   tooltips; otherwise a static `ggplot`.
#' @param width_svg girafe canvas width (interactive only).
#' @returns A `girafe` (interactive) or `ggplot` (static) object.
#' @export
hdwgcna_enrichr_barplot <- function(enrichr_df, db, modules = NULL, n_terms = 6,
                                    interactive = TRUE, width_svg = 10) {
  d <- enrichr_df[enrichr_df$db == db & enrichr_df$module != "grey", , drop = FALSE]
  if (!is.null(modules)) d <- d[d$module %in% modules, , drop = FALSE]
  if (nrow(d) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
        label = paste0("No GO terms to display for ", db,
          "\n(no module had enrichment results)")) +
      ggplot2::theme_void()
    return(if (isTRUE(interactive)) ggiraph::girafe(ggobj = p) else p)
  }
  d$neg_log10 <- -log10(pmax(d$Adjusted.P.value, .Machine$double.xmin))
  d <- do.call(rbind, lapply(split(d, d$module), function(m) {
    utils::head(m[order(m$Adjusted.P.value), , drop = FALSE], n_terms)
  }))
  d$module <- factor(d$module, levels = sort(unique(as.character(d$module))))
  d$Term_w <- .hdwgcna_wrap_terms(d$Term)
  d$row_id <- seq_len(nrow(d))
  d$tooltip <- paste0(
    "Module: ", d$module,
    "\nTerm: ", d$Term,
    "\nOverlap: ", d$Overlap,
    "\nraw p-value: ", signif(d$P.value, 3),
    "\nadj. p-value: ", signif(d$Adjusted.P.value, 3),
    "\nCombined score: ", signif(d$Combined.Score, 3),
    "\nGenes: ", d$Genes
  )
  bar <- if (isTRUE(interactive)) {
    ggiraph::geom_col_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = row_id), fill = "cornflowerblue"
    )
  } else {
    ggplot2::geom_col(fill = "cornflowerblue")
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = neg_log10, y = .hdwgcna_reorder_within(Term_w, neg_log10, module)
  )) +
    bar +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
      color = "red", linewidth = 0.4) +
    ggplot2::facet_wrap(~module, scales = "free", ncol = 2) +
    .hdwgcna_scale_y_reordered() +
    ggplot2::labs(x = "- log10 adjusted P-value", y = NULL,
      title = db, subtitle = "dashed line: adjusted p = 0.05") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (!isTRUE(interactive)) {
    return(p)
  }
  n_mod <- length(unique(d$module))
  ggiraph::girafe(
    ggobj = p, width_svg = width_svg, height_svg = max(4, ceiling(n_mod / 2) * 2.4),
    options = list(
      ggiraph::opts_hover(css = "fill:orange;stroke:black;"),
      ggiraph::opts_tooltip(css = "background-color:white;color:black;padding:6px;border:1px solid grey;")
    )
  )
}

#' Dot plot of module eigengenes across cell types
#'
#' Adds the harmonised module eigengenes (hMEs) of the analysed group's modules
#' to the metadata and draws a `Seurat::DotPlot` across `group_by`.
#'
#' @param obj An hdWGCNA Seurat object.
#' @param group_by Metadata column to group the dot plot by (e.g. the cell-type
#'   clustering column).
#' @returns A `ggplot` object.
#' @export
hdwgcna_me_dotplot <- function(obj, group_by) {
  .hdwgcna_require()
  MEs <- hdWGCNA::GetMEs(obj, harmonized = TRUE)
  modules <- hdWGCNA::GetModules(obj)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  MEs <- MEs[rownames(obj@meta.data), , drop = FALSE]
  obj@meta.data <- cbind(obj@meta.data, MEs)
  Seurat::DotPlot(obj, features = mods, group.by = group_by) +
    Seurat::RotatedAxis() +
    ggplot2::scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
}

#' Per-cell module scores via Seurat::AddModuleScore
#'
#' Scores **every cell of the object** (all cell types, not a single-cell-type
#' subset) for each module of one hdWGCNA network with [Seurat::AddModuleScore()]
#' (the average expression of a module's genes minus a matched random control
#' set). The modules are those learned in one cell group, but the scores are
#' computed across all cells so a module's activity can be compared everywhere
#' and grouped by clinical condition downstream.
#'
#' @param obj An hdWGCNA Seurat object (the `wgcna_<group>` target), which holds
#'   the full set of cells passed to [tar_hdwgcna()] plus that group's modules and
#'   the (log-normalised) `assay`.
#' @param trait_col Optional clinical-trait metadata column; copied into the
#'   returned `condition` column so the scores can be grouped by clinical
#'   condition downstream. `NULL` leaves `condition` as `NA`.
#' @param clustering_col Optional cell-type / cluster metadata column; copied into
#'   the returned `cell_type` column so each cell keeps its own cell-type label
#'   (the scores still cover all cells). `NULL` leaves `cell_type` as `NA`.
#' @param assay Assay used for scoring (default `"RNA"`, log-normalised).
#' @param include_grey Score the `"grey"` module too (default `FALSE`).
#' @param nbin,ctrl Passed to [Seurat::AddModuleScore()].
#' @returns A `data.frame` with one row per cell in the object: `cell`,
#'   `cell_type`, `condition`, then one numeric column per module (named after the
#'   module). An empty `data.frame` if there is no non-grey module.
#' @export
hdwgcna_module_scores <- function(obj, trait_col = NULL, clustering_col = NULL,
                                  assay = "RNA", include_grey = FALSE,
                                  nbin = 24, ctrl = 100) {
  .hdwgcna_require()
  gene_list <- hdwgcna_module_gene_list(obj, include_grey = include_grey)
  if (length(gene_list) == 0L) {
    warning("No module to score.")
    return(data.frame())
  }
  # Score all cells of the object (every cell type), not a single-cell-type subset.
  SeuratObject::DefaultAssay(obj) <- assay
  obj <- tryCatch(
    Seurat::AddModuleScore(obj, features = gene_list, name = "modulescore_",
                           assay = assay, nbin = nbin, ctrl = ctrl),
    error = function(e) {
      warning("AddModuleScore failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(obj)) return(data.frame())
  score_cols <- paste0("modulescore_", seq_along(gene_list))
  scores <- obj@meta.data[, score_cols, drop = FALSE]
  colnames(scores) <- names(gene_list)
  out <- data.frame(
    cell      = rownames(obj@meta.data),
    cell_type = if (!is.null(clustering_col) && clustering_col %in% colnames(obj@meta.data)) {
      as.character(obj@meta.data[[clustering_col]])
    } else {
      NA_character_
    },
    condition = if (!is.null(trait_col) && trait_col %in% colnames(obj@meta.data)) {
      as.character(obj@meta.data[[trait_col]])
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cbind(out, scores)
}

#' Heatmap of per-cell module scores grouped by clinical condition
#'
#' Tile heatmap of the per-cell module scores from [hdwgcna_module_scores()]:
#' rows are modules, columns are individual cells **grouped (faceted) by clinical
#' condition**, fill is the (row-scaled) module score. Built with `ggplot2` for
#' consistency with the other report heatmaps.
#'
#' @param score_df The per-cell score table from [hdwgcna_module_scores()] (a
#'   `wgcna_<group>_modscore` target).
#' @param trait_col Column of `score_df` holding the clinical condition used to
#'   group the cells (default `"condition"`, as set by [hdwgcna_module_scores()]).
#' @param module_cols Module columns to plot; `NULL` (default) uses every column
#'   that is not `cell`/`cell_type`/`condition`.
#' @param scale_rows Z-score each module across the displayed cells (default
#'   `TRUE`) so modules on different scales are comparable.
#' @param max_cells_per_group Down-sample to at most this many cells per clinical
#'   condition to keep the figure legible (default 300; `NULL` keeps all cells).
#' @param seed Seed for the down-sampling (default 1).
#' @param celltype_col Column of `score_df` holding each cell's cell type, used to
#'   order the cells and draw the cell-type annotation bar (default `"cell_type"`).
#' @param annotate_celltype Add a cell-type colour annotation bar above the heatmap
#'   and order the cells by cell type within each clinical condition (default
#'   `TRUE`). Ignored when `celltype_col` is missing or has a single value.
#' @returns A `ggplot` object (or a `patchwork` of the annotation bar + heatmap
#'   when `annotate_celltype = TRUE`).
#' @export
hdwgcna_module_score_heatmap <- function(score_df, trait_col = "condition",
                                         module_cols = NULL, scale_rows = TRUE,
                                         max_cells_per_group = 300, seed = 1,
                                         celltype_col = "cell_type",
                                         annotate_celltype = TRUE) {
  empty_plot <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg) +
      ggplot2::theme_void()
  }
  if (is.null(score_df) || !is.data.frame(score_df) || nrow(score_df) == 0L) {
    return(empty_plot("No module scores to display"))
  }
  meta_cols <- unique(c("cell", "cell_type", "condition", trait_col, celltype_col))
  if (is.null(module_cols)) module_cols <- setdiff(colnames(score_df), meta_cols)
  if (length(module_cols) == 0L) return(empty_plot("No module-score columns found"))

  df <- score_df[!is.na(score_df[[trait_col]]), , drop = FALSE]
  if (nrow(df) == 0L) return(empty_plot("No cell has a clinical condition"))

  # down-sample per clinical condition for legibility
  if (!is.null(max_cells_per_group)) {
    set.seed(seed)
    keep <- unlist(lapply(split(seq_len(nrow(df)), df[[trait_col]]), function(idx) {
      if (length(idx) > max_cells_per_group) sample(idx, max_cells_per_group) else idx
    }), use.names = FALSE)
    df <- df[sort(keep), , drop = FALSE]
  }

  mat <- as.matrix(df[, module_cols, drop = FALSE])
  if (isTRUE(scale_rows)) {
    mat <- apply(mat, 2, function(x) {
      s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
    })
  }

  cond <- as.character(df[[trait_col]])
  has_ct <- isTRUE(annotate_celltype) && celltype_col %in% colnames(df) &&
    length(unique(stats::na.omit(df[[celltype_col]]))) > 1L
  ct <- if (has_ct) as.character(df[[celltype_col]]) else NULL

  # order cells by condition, then by cell type within each condition so that
  # the facet holds contiguous, cell-type-sorted blocks of cells
  ord <- if (has_ct) order(cond, ct) else order(cond)
  cells_ordered <- df$cell[ord]

  long <- data.frame(
    cell      = factor(rep(df$cell, times = length(module_cols)), levels = cells_ordered),
    module    = factor(rep(module_cols, each = nrow(df)), levels = rev(module_cols)),
    score     = as.vector(mat),
    condition = factor(rep(cond, times = length(module_cols))),
    stringsAsFactors = FALSE
  )
  lim <- max(abs(long$score), na.rm = TRUE)

  heat <- ggplot2::ggplot(long, ggplot2::aes(x = cell, y = module, fill = score)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(cols = ggplot2::vars(condition), scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "grey95", high = "red", midpoint = 0,
      limits = c(-lim, lim),
      name = if (isTRUE(scale_rows)) "Module score\n(row z-score)" else "Module score"
    ) +
    ggplot2::labs(x = "Cells (grouped by clinical condition, ordered by cell type)", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(2, "pt")
    )

  if (!has_ct) {
    return(heat + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold")))
  }

  # cell-type annotation bar (kept divided by clinical condition); the condition
  # strip labels sit on this top track, so they are hidden on the heatmap below.
  ann_df <- data.frame(
    cell      = factor(df$cell, levels = cells_ordered),
    cell_type = ct,
    condition = factor(cond),
    stringsAsFactors = FALSE
  )
  ann <- ggplot2::ggplot(ann_df, ggplot2::aes(x = cell, y = "Cell type", fill = cell_type)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(cols = ggplot2::vars(condition), scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_discrete(name = "Cell type") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(2, "pt"),
      strip.text = ggplot2::element_text(face = "bold")
    )
  heat <- heat + ggplot2::theme(strip.text = ggplot2::element_blank())

  patchwork::wrap_plots(ann, heat, ncol = 1, heights = c(1, 10), guides = "collect")
}

#' Violin plot of one module's per-cell scores across clinical groups
#'
#' Draws the per-cell `AddModuleScore` of a **single module** as violins across
#' the clinical groups (x-axis, one violin + inner boxplot per group, filled by
#' group), **faceted by cell type**. Pairwise comparisons between the clinical
#' groups are annotated with [ggpubr::stat_compare_means()] (all pairs, unpaired
#' Wilcoxon rank-sum test by default -- an appropriate non-parametric test for
#' the per-cell scores). Companion to [hdwgcna_module_score_heatmap()]: the
#' heatmap shows every module at once, this drills into one module cell-type by
#' cell-type. Meant to be looped over the module columns (one tab per module).
#'
#' Only cell types where **every** displayed clinical group has at least
#' `min_cells` cells are kept, so every pairwise test is computable in every
#' facet (rare cell types are dropped).
#'
#' @param score_df The per-cell score table from [hdwgcna_module_scores()] (a
#'   `wgcna_<group>_modscore` target).
#' @param module Name of the single module column of `score_df` to plot.
#' @param trait_col Column of `score_df` holding the clinical group (default
#'   `"condition"`, as set by [hdwgcna_module_scores()]).
#' @param celltype_col Column of `score_df` holding each cell's cell type
#'   (default `"cell_type"`); used for the facets.
#' @param trait_groups Optional ordered clinical-group levels (as in
#'   [tar_hdwgcna()]); sets the x-axis order and the pairwise comparisons. `NULL`
#'   uses the sorted unique conditions.
#' @param test Test passed to [ggpubr::stat_compare_means()]: `"wilcox.test"`
#'   (default, unpaired Wilcoxon) or `"t.test"`.
#' @param ncol Number of facet columns (default 2).
#' @param min_cells Minimum cells per clinical group required for a cell type to
#'   be shown (default 10).
#' @param scales Facet scales (default `"free_y"`, so each cell type is legible).
#' @returns A `ggplot` object (a placeholder plot with a message if there is
#'   nothing to show).
#' @export
hdwgcna_module_score_violin <- function(score_df, module,
                                        trait_col = "condition",
                                        celltype_col = "cell_type",
                                        trait_groups = NULL,
                                        test = c("wilcox.test", "t.test"),
                                        ncol = 2, min_cells = 10,
                                        scales = "free_y") {
  empty_plot <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg) +
      ggplot2::theme_void()
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required for hdwgcna_module_score_violin().", call. = FALSE)
  }
  test <- match.arg(test)
  if (is.null(score_df) || !is.data.frame(score_df) || nrow(score_df) == 0L) {
    return(empty_plot("No module scores to display"))
  }
  if (!module %in% colnames(score_df)) {
    return(empty_plot(sprintf("Module '%s' not found", module)))
  }

  df <- data.frame(
    score     = as.numeric(score_df[[module]]),
    group     = as.character(score_df[[trait_col]]),
    cell_type = if (celltype_col %in% colnames(score_df)) {
      as.character(score_df[[celltype_col]])
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$score) & !is.na(df$group), , drop = FALSE]
  if (nrow(df) == 0L) return(empty_plot("No cell has a clinical group"))

  grp_levels <- if (is.null(trait_groups)) {
    sort(unique(df$group))
  } else {
    trait_groups[trait_groups %in% df$group]
  }
  if (length(grp_levels) < 2L) {
    return(empty_plot(sprintf("Module '%s': fewer than two clinical groups present", module)))
  }
  df <- df[df$group %in% grp_levels, , drop = FALSE]
  df$group <- factor(df$group, levels = grp_levels)

  if (all(is.na(df$cell_type))) df$cell_type <- "all cells"

  # keep cell types where every displayed group has >= min_cells cells, so all
  # pairwise comparisons are computable in every facet
  keep <- vapply(split(df, df$cell_type), function(d) {
    all(table(factor(d$group, levels = grp_levels)) >= min_cells)
  }, logical(1))
  keep_names <- names(keep)[keep]
  if (length(keep_names) == 0L) {
    return(empty_plot(sprintf(
      "Module '%s': no cell type has ≥ %d cells in every clinical group",
      module, min_cells)))
  }
  df <- df[df$cell_type %in% keep_names, , drop = FALSE]
  df$cell_type <- factor(df$cell_type, levels = sort(keep_names))

  comparisons <- utils::combn(grp_levels, 2, simplify = FALSE)

  ggplot2::ggplot(df, ggplot2::aes(x = group, y = score, fill = group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.3) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                          alpha = 0.7, linewidth = 0.3) +
    ggpubr::stat_compare_means(comparisons = comparisons, method = test,
                               size = 2.7, tip.length = 0.01) +
    ggplot2::facet_wrap(~ cell_type, ncol = ncol, scales = scales) +
    ggplot2::labs(x = NULL, y = sprintf("%s module score", module),
                  fill = "Clinical group", title = module) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' avg_log2FC heatmap of differential module eigengenes
#'
#' A tile heatmap of module-eigengene `avg_log2FC` across the pairwise clinical
#' comparisons (columns = modules, rows = comparisons; significance stars from
#' `p_val_adj`). Positive = higher in the test group.
#'
#' @param obj An hdWGCNA Seurat object (for the module ordering).
#' @param dmes The combined `FindDMEs` table (a `wgcna_<group>_dmes` target) with
#'   `module`, `comparison`, `avg_log2FC` and `p_val_adj` columns.
#' @param clip Symmetric clipping applied to `avg_log2FC` for the colour scale
#'   (default 0.5).
#' @returns A `ggplot` object.
#' @export
hdwgcna_dme_heatmap <- function(obj, dmes, clip = 0.5) {
  .hdwgcna_require()
  mods <- levels(hdWGCNA::GetModules(obj)$module)
  mods <- mods[mods != "grey"]
  plot_df <- dmes
  plot_df$module <- factor(as.character(plot_df$module), levels = mods)
  plot_df$avg_log2FC <- pmax(pmin(plot_df$avg_log2FC, clip), -clip)
  ss <- function(p) {
    ifelse(is.na(p), "",
      ifelse(p < 0.001, "***",
        ifelse(p < 0.01, "**",
          ifelse(p < 0.05, "*", ""))))
  }
  plot_df$Significance <- ss(plot_df$p_val_adj)
  plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, "black", "white")
  ggplot2::ggplot(plot_df, ggplot2::aes(y = comparison, x = module, fill = avg_log2FC)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(label = plot_df$Significance, color = plot_df$textcolor) +
    ggplot2::scale_fill_gradient2(low = "purple", mid = "black", high = "yellow") +
    Seurat::RotatedAxis() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    ) +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_equal()
}

#' Build the per-cell-type hdWGCNA report section
#'
#' Generates the Quarto/knitr lines for one cell group's hdWGCNA results --
#' soft-power, dendrogram + modules, module-trait correlation, Enrichr GO
#' enrichment, the module-eigengene dot plot and the differential module
#' eigengenes -- as a single tabset under a `## <group> cells` heading. Mirrors
#' [dea_report_lines()]: emit the returned lines with `knitr::knit_child()` (see
#' that function's usage). The generated chunks read the factory targets
#' (`wgcna_<suffix>`, `_powertest`, `_soft_power`, `_dmes`, `_enrichr`) directly
#' via `targets::tar_read()`, so nothing needs to be pre-bound in the knit
#' environment; the cell-type clustering column / trait must match the
#' [tar_hdwgcna()] call that produced the targets.
#'
#' @param group Cell-group label, exactly as passed to [tar_hdwgcna()] (e.g.
#'   `"MAIT"` or `"CD4 TEM"`). The per-group target names are derived from it.
#' @param clustering_col Cell-type / cluster metadata column used for the
#'   module-eigengene dot plot (default `"predicted.celltype.l2"`).
#' @param trait_col Clinical-trait metadata column (for the prose). When `NULL`,
#'   the module-trait and differential-module-eigengene sections are omitted.
#' @param trait_groups The clinical-group levels (as in [tar_hdwgcna()]). Used to
#'   build the trait labels and the pairwise DME comparisons.
#' @param heading_level Markdown heading level for the `<group> cells` title
#'   (default 2; tabs are one level deeper, inner tabsets two deeper).
#' @param enrichr_dbs Enrichr databases to show, one tab each (default the GO
#'   BP/CC/MF 2023 databases). Use `character()` to omit the enrichment section.
#' @param interactive_enrichr Whether the Enrichr barplots are interactive
#'   (`ggiraph`, default `TRUE`).
#' @returns A character vector of Quarto/knitr lines.
#' @export
hdwgcna_report_lines <- function(group,
                                 clustering_col = "predicted.celltype.l2",
                                 trait_col = NULL,
                                 trait_groups = NULL,
                                 heading_level = 2L,
                                 enrichr_dbs = c(
                                   "GO_Biological_Process_2023",
                                   "GO_Cellular_Component_2023",
                                   "GO_Molecular_Function_2023"
                                 ),
                                 interactive_enrichr = TRUE) {
  stopifnot(is.character(group), length(group) == 1L, nzchar(group))
  sfx <- .hdwgcna_suffix(group)
  net <- paste0("wgcna_", sfx)
  pt <- paste0("wgcna_", sfx, "_powertest")
  sp <- paste0("wgcna_", sfx, "_soft_power")
  dme <- paste0("wgcna_", sfx, "_dmes")
  enr <- paste0("wgcna_", sfx, "_enrichr")
  ms <- paste0("wgcna_", sfx, "_modscore")
  h1 <- .dea_h(heading_level)
  h2 <- .dea_h(heading_level + 1L)
  h3 <- .dea_h(heading_level + 2L)
  has_trait <- !is.null(trait_col) && !is.null(trait_groups) && length(trait_groups) >= 2L
  int_chr <- if (isTRUE(interactive_enrichr)) "TRUE" else "FALSE"

  tl_expr <- if (has_trait) {
    tl <- stats::setNames(trait_groups, paste0("is_", gsub("[^A-Za-z0-9]+", "", trait_groups)))
    paste(deparse(tl), collapse = " ")
  } else {
    "NULL"
  }

  lines <- c(
    sprintf("%s %s cells", h1, group), "",
    "::: {.panel-tabset}", "",
    "```{r}",
    "#| echo: false",
    "#| message: false",
    "#| warning: false",
    sprintf('wgcna_obj <- targets::tar_read("%s")', net),
    sprintf("trait_labels <- %s", tl_expr),
    "```",
    "",
    sprintf("%s Soft-power threshold", h2), "",
    sprintf('A soft-power of `r targets::tar_read("%s")` was chosen based on hdWGCNA computations.', sp), "",
    "```{r}",
    "#| message: false",
    "#| warning: false",
    sprintf("#| label: fig-softpower-%s", sfx),
    sprintf('p <- hdWGCNA::PlotSoftPowers(targets::tar_read("%s"))', pt),
    "patchwork::wrap_plots(p, ncol = 2)",
    "```",
    "",
    sprintf("%s Gene modules", h2), "",
    "```{r}",
    sprintf("#| label: fig-dendro-%s", sfx),
    sprintf('#| fig-cap: "%s hdWGCNA dendrogram"', group),
    'hdWGCNA::PlotDendrogram(wgcna_obj, main = "")',
    "```",
    "",
    "```{r}",
    sprintf("#| label: tbl-modules-%s", sfx),
    sprintf('#| tbl-cap: "Modules identified for %s cells"', group),
    "hdWGCNA::GetModules(wgcna_obj) |>",
    "  dplyr::group_by(module, color) |>",
    '  dplyr::summarise(N = dplyr::n(), .groups = "drop") |>',
    "  DT2::dt2()",
    "```",
    ""
  )

  if (has_trait) {
    lines <- c(lines,
      sprintf("%s Modules associated with the clinical groups", h2), "",
      sprintf(paste0(
        "Each %s gene module is summarised by its harmonised module eigengene ",
        "(hME) and correlated (Spearman) with one-hot indicators of each clinical ",
        "group (`%s`) across all %s cells; p-values are FDR-corrected. Tile colour ",
        "= correlation `r`, tile text = FDR (with significance stars)."
      ), group, trait_col, group), "",
      "```{r}",
      sprintf("#| label: fig-mtcor-%s", sfx),
      "#| fig-width: 10",
      "#| fig-height: 3.5",
      "#| message: false",
      "#| warning: false",
      "scitargets::hdwgcna_module_trait_heatmap(wgcna_obj, trait_labels = trait_labels)",
      "```",
      "",
      "```{r}",
      sprintf("#| label: tbl-mtcor-%s", sfx),
      sprintf('#| tbl-cap: "%s module / clinical-group correlations (r) and FDR, sorted by FDR."', group),
      "DT2::dt2(scitargets::hdwgcna_module_trait_table(wgcna_obj, trait_labels = trait_labels))",
      "```",
      "",
      sprintf("%s Module scores per cell, grouped by clinical condition", h2), "",
      sprintf(paste0(
        "Per-cell module scores (`Seurat::AddModuleScore`) for the %s modules, ",
        "computed across **all cells of the object (every cell type)**. Each column ",
        "is a cell, **grouped (faceted) by clinical condition (`%s`)**; each row is a ",
        "module (scores are row-scaled to a z-score). Up to 300 cells per condition ",
        "are shown."
      ), group, trait_col), "",
      "```{r}",
      sprintf("#| label: fig-modscore-%s", sfx),
      "#| fig-width: 10",
      "#| fig-height: 5",
      "#| message: false",
      "#| warning: false",
      sprintf('scitargets::hdwgcna_module_score_heatmap(targets::tar_read("%s"), trait_col = "condition")', ms),
      "```",
      "",
      sprintf("%s Module-score violins per module (by cell type)", h2), "",
      sprintf(paste0(
        "The same per-cell `AddModuleScore` values, drilled into one module at a ",
        "time (one tab per module). Each panel is a cell type; within a panel the ",
        "scores are shown as violins across the clinical groups (`%s`), filled by ",
        "group. Pairwise group comparisons use an unpaired Wilcoxon rank-sum test ",
        "(`ggpubr::stat_compare_means`). Only cell types with at least 10 cells in ",
        "**every** clinical group are shown."
      ), trait_col), "",
      "```{r}",
      sprintf("#| label: load-modscore-%s", sfx),
      "#| echo: false",
      "#| message: false",
      "#| warning: false",
      sprintf('modscore_df <- targets::tar_read("%s")', ms),
      "```",
      ""
    )
    # One static tab per module (real headings, like the Enrichr tabset below).
    # The module list and a sensible facet height are read from the modscore
    # target now -- targets is available when the report loop runs; a read
    # failure simply omits the violin tabs.
    ms_info <- tryCatch({
      .df <- targets::tar_read_raw(ms)
      .mods <- setdiff(colnames(.df), c("cell", "cell_type", "condition"))
      .sub <- .df[!is.na(.df$cell_type) & .df$condition %in% trait_groups, , drop = FALSE]
      .cnts <- table(.sub$cell_type, factor(.sub$condition, levels = trait_groups))
      .n_kept <- sum(apply(.cnts, 1, function(rw) all(rw >= 10)))
      list(mods = .mods, figh = max(4, 2.6 * ceiling(max(.n_kept, 1) / 2)))
    }, error = function(e) list(mods = character(0), figh = 8))
    if (length(ms_info$mods) > 0L) {
      tg_expr <- paste(deparse(trait_groups), collapse = " ")
      lines <- c(lines, "::: {.panel-tabset}", "")
      for (m in ms_info$mods) {
        lines <- c(lines,
          sprintf("%s %s", h3, m), "",
          "```{r}",
          sprintf("#| label: fig-msviol-%s-%s", sfx, gsub("[^A-Za-z0-9]+", "_", m)),
          "#| echo: false",
          "#| message: false",
          "#| warning: false",
          "#| fig-width: 9",
          sprintf("#| fig-height: %.1f", ms_info$figh),
          sprintf('#| fig-cap: "%s module score across clinical groups, by cell type (pairwise Wilcoxon)."', m),
          sprintf('scitargets::hdwgcna_module_score_violin(modscore_df, module = "%s", trait_groups = %s, min_cells = 10)', m, tg_expr),
          "```",
          ""
        )
      }
      lines <- c(lines, ":::", "")
    }
  }

  lines <- c(lines,
    "```{r}",
    sprintf("#| label: tbl-hub-%s", sfx),
    sprintf('#| tbl-cap: "Top 10 hub genes per %s module (highest kME)."', group),
    "DT2::dt2(hdWGCNA::GetHubGenes(wgcna_obj, n_hubs = 10))",
    "```",
    ""
  )

  if (length(enrichr_dbs) > 0L) {
    db_short <- c(
      GO_Biological_Process_2023 = "bp",
      GO_Cellular_Component_2023 = "cc",
      GO_Molecular_Function_2023 = "mf"
    )
    lines <- c(lines,
      sprintf("%s Module functional enrichment (Enrichr GO)", h2), "",
      "```{r}",
      "#| echo: false",
      "#| message: false",
      "#| warning: false",
      sprintf('enrichr_obj <- targets::tar_read("%s")', enr),
      'go_modules <- if (is.data.frame(enrichr_obj) && nrow(enrichr_obj) > 0L) sort(unique(as.character(enrichr_obj$module[enrichr_obj$module != "grey"]))) else character(0)',
      "```",
      "",
      paste0(
        "GO enrichment was computed with `RunEnrichr` for the ",
        "`r length(go_modules)` module(s) selected by the pipeline ",
        "(set via `enrichr_modules` in the `tar_hdwgcna()` call): ",
        "`r if (length(go_modules)) paste(go_modules, collapse = \", \") else \"none\"`. ",
        "Every module present in the enrichment table is shown below — there is ",
        "no further display filter. Bar length = `-log10(adjusted p-value)`; the ",
        "dashed red line marks the FDR = 0.05 cutoff. Hover a bar for the term, ",
        "overlap, p-values, combined score and genes."
      ), "",
      "::: {.panel-tabset}", "")
    for (db in enrichr_dbs) {
      shortid <- if (db %in% names(db_short)) db_short[[db]] else gsub("[^a-z0-9]+", "", tolower(db))
      pretty <- gsub("_", " ", db)
      lines <- c(lines,
        sprintf("%s %s", h3, pretty), "",
        "```{r}",
        sprintf("#| label: fig-enrichr-%s-%s", shortid, sfx),
        "#| message: false",
        "#| warning: false",
        sprintf('scitargets::hdwgcna_enrichr_barplot(enrichr_obj, "%s", modules = NULL, interactive = %s)', db, int_chr),
        "```",
        ""
      )
    }
    lines <- c(lines, ":::", "")
  }

  lines <- c(lines,
    sprintf("%s Module eigengenes across cell types", h2), "",
    sprintf("Harmonised module eigengenes (hMEs) of the %s-derived modules across cell types.", group), "",
    "```{r}",
    sprintf("#| label: fig-medot-%s", sfx),
    "#| fig-width: 10",
    "#| fig-height: 6",
    "#| message: false",
    "#| warning: false",
    sprintf('scitargets::hdwgcna_me_dotplot(wgcna_obj, group_by = "%s")', clustering_col),
    "```",
    ""
  )

  if (has_trait) {
    comps <- utils::combn(trait_groups, 2, simplify = FALSE)
    lines <- c(lines,
      sprintf("%s Differential module eigengenes between clinical groups", h2), "",
      sprintf(paste0(
        "Modules whose hME differs between two clinical groups (`FindDMEs`, ",
        "**%s cells only**). Positive avg_log2FC = higher in the first (test) ",
        "group; the second group is the reference / control."
      ), group), "",
      "```{r}",
      "#| echo: false",
      "#| message: false",
      "#| warning: false",
      sprintf('dmes_obj <- targets::tar_read("%s")', dme),
      "```",
      "",
      "::: {.panel-tabset}", ""
    )
    for (i in seq_along(comps)) {
      g1 <- comps[[i]][1]
      g2 <- comps[[i]][2]
      cmp <- paste(g1, "vs", g2)
      lines <- c(lines,
        sprintf("%s %s", h3, cmp), "",
        sprintf("Test group: **%s** &nbsp;•&nbsp; reference / control: **%s**.", g1, g2), "",
        "```{r}",
        sprintf("#| label: fig-dmevolc-%d-%s", i, sfx),
        "#| fig-width: 6",
        "#| fig-height: 5",
        "#| message: false",
        "#| warning: false",
        sprintf('hdWGCNA::PlotDMEsVolcano(wgcna_obj, dmes_obj[dmes_obj$comparison == "%s", ])', cmp),
        "```",
        ""
      )
    }
    lines <- c(lines, ":::", "",
      sprintf("%s avg_log2FC heatmap (all comparisons)", h3), "",
      "avg_log2FC of each module across the comparisons (`*` FDR < 0.05, `**` < 0.01, `***` < 0.001).", "",
      "```{r}",
      sprintf("#| label: fig-dmeheat-%s", sfx),
      "#| fig-width: 9",
      "#| fig-height: 4",
      "#| message: false",
      "#| warning: false",
      "scitargets::hdwgcna_dme_heatmap(wgcna_obj, dmes_obj)",
      "```",
      ""
    )
  }

  c(lines, ":::", "")
}

utils::globalVariables(c(
  "module", "clinical_group", "r", "label", "neg_log10", "Term_w",
  "tooltip", "row_id", "comparison", "avg_log2FC", "score", "condition", "cell",
  "cell_type", "group"
))
