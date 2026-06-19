# R/scitargets_dea.R
# S7 result object for differential expression, GO and GSEA enrichment of
# clinical-group comparisons in single-cell RNA-seq data.
#
# The class holds several "levels" of differential expression for the same
# comparison:
#   * "single_cell" : Seurat::FindMarkers on the individual cells
#   * "pseudobulk"   : DESeq2 on per-replicate (e.g. patient) pseudobulk
# several GO ontologies (BP / CC / MF) and several GSEA collections per
# comparison. S7 methods are registered in .onLoad() (see zzz.R).

DataFrame <- S7::new_S3_class("data.frame")

# -----------------------------------------------------------------------------
# Class definition
# -----------------------------------------------------------------------------

#' Differential-expression / GO / GSEA result object
#'
#' An S7 class holding the differential-expression results (single-cell and/or
#' pseudobulk), GO enrichment (one or more ontologies) and GSEA results for a
#' single clinical-group comparison within one cluster. Build instances with
#' [run_dea()].
#'
#' @param comparison_name Length-1 label for the comparison ("group1::group2::cluster").
#' @param group_by Metadata column holding the group labels.
#' @param cluster_by Metadata column holding the clustering (length 0 or 1).
#' @param cluster The cluster value this comparison sits in.
#' @param groups The two group labels being compared.
#' @param group1,group2 The two groups (`group2 = "rest"` means one-vs-rest).
#' @param levels Computed differential-expression levels (subset of
#'   `c("single_cell", "pseudobulk")`).
#' @param pseudobulk_unit Biological-replicate column used for the pseudobulk level.
#' @param n_cells Per-comparison-group cell counts (`data.frame`).
#' @param de Named list (level -> standardized marker `data.frame`).
#' @param go Named list (level -> direction -> ontology -> topGO result).
#' @param gsea Named list (level -> collection -> fgsea result).
#' @param pca Named list (level -> pseudobulk VST PCA + outlier results).
#' @param de_params,go_params,gsea_params Parameter lists recording how each
#'   analysis was run.
#' @param species Organism (`"human"` or `"mouse"`).
#' @param padj_cutoffs Named list of per-level adjusted-p cutoffs.
#' @param status Either `"computed"` or `"not_computed"`.
#' @param reason Length-1 explanation when `status == "not_computed"`.
#' @returns An S7 `scitargets_dea` object generator.
#' @export
scitargets_dea <- S7::new_class(
  "scitargets_dea",
  properties = list(
    comparison_name = S7::class_character,
    group_by = S7::class_character, # metadata column holding the groups
    cluster_by = S7::class_character, # metadata column holding the clustering
    cluster = S7::class_character, # cluster value this comparison sits in
    groups = S7::class_character, # the two group labels being compared
    group1 = S7::class_character,
    group2 = S7::class_character,
    levels = S7::class_character, # subset of c("single_cell","pseudobulk")
    pseudobulk_unit = S7::class_character, # replicate column for pseudobulk (e.g. patient_id)
    n_cells = S7::new_property(DataFrame, default = quote(data.frame())),
    de = S7::class_list, # level -> standardized marker data.frame
    go = S7::class_list, # level -> direction -> ontology -> topGO result
    gsea = S7::class_list, # level -> collection -> fgsea result (stage 3)
    pca = S7::new_property(S7::class_list, default = list()), # level -> VST PCA + outliers (pseudobulk QC)
    de_params = S7::class_list,
    go_params = S7::class_list,
    gsea_params = S7::class_list,
    species = S7::new_property(S7::class_character, default = "human"), # "human" / "mouse"
    padj_cutoffs = S7::new_property( # per-level adjusted-p cutoff:
      S7::class_list, # volcano horizontal line + GO
      default = quote(list(single_cell = 0.05, pseudobulk = 0.05)) # foreground genes
    ),
    status = S7::class_character,
    reason = S7::class_character
  ),
  validator = function(self) {
    problems <- character()

    if (length(self@comparison_name) != 1L) {
      problems <- c(problems, "@comparison_name must be length 1.")
    }
    if (length(self@groups) != 2L) {
      problems <- c(problems, "@groups must contain exactly two group labels.")
    }
    if (length(self@group1) != 1L) {
      problems <- c(problems, "@group1 must be length 1.")
    }
    if (length(self@group2) != 1L) {
      problems <- c(problems, "@group2 must be length 1.")
    }
    if (length(self@cluster) != 1L) {
      problems <- c(problems, "@cluster must be length 1.")
    }
    if (length(self@group_by) != 1L) {
      problems <- c(problems, "@group_by must be length 1.")
    }
    if (length(self@cluster_by) > 1L) {
      problems <- c(problems, "@cluster_by must be length 0 or 1.")
    }
    valid_levels <- c("single_cell", "pseudobulk")
    if (length(self@levels) > 0L && !all(self@levels %in% valid_levels)) {
      problems <- c(
        problems,
        "@levels must only contain 'single_cell' and/or 'pseudobulk'."
      )
    }
    if (length(self@status) != 1L || !self@status %in% c("computed", "not_computed")) {
      problems <- c(problems, "@status must be either 'computed' or 'not_computed'.")
    }
    if (length(self@reason) != 1L) {
      problems <- c(problems, "@reason must be length 1.")
    }
    if (length(self@species) != 1L || !self@species %in% c("human", "mouse")) {
      problems <- c(problems, "@species must be 'human' or 'mouse' (only human and mouse are supported).")
    }
    if (!is.list(self@pca)) {
      problems <- c(problems, "@pca must be a list.")
    }
    if (!is.list(self@padj_cutoffs) ||
      !all(vapply(
        self@padj_cutoffs,
        function(v) is.numeric(v) && length(v) == 1L && v > 0 && v <= 1,
        logical(1L)
      ))) {
      problems <- c(problems, "@padj_cutoffs must be a list of numeric scalars in (0, 1].")
    }

    if (length(problems) > 0L) problems else NULL
  }
)

# -----------------------------------------------------------------------------
# List normalization helper for tar_read() output
# -----------------------------------------------------------------------------

#' Test whether an object is a [scitargets_dea]
#'
#' @param x Any object.
#' @returns A logical scalar.
#' @export
is_scitargets_dea <- function(x) {
  # S7_inherits is the reliable check: inside a package the S7 class name is
  # namespaced ("scitargets::scitargets_dea"), so a bare inherits() on
  # "scitargets_dea" returns FALSE.
  isTRUE(tryCatch(S7::S7_inherits(x, scitargets_dea), error = function(e) FALSE)) ||
    inherits(x, "scitargets_dea") ||
    inherits(x, "scitargets::scitargets_dea")
}

#' Flatten a (possibly nested) list of [scitargets_dea] objects
#'
#' Normalizes the output of `targets::tar_read()` (a list of branches) into a
#' flat named list keyed by comparison name.
#'
#' @param x A [scitargets_dea], or an arbitrarily nested list of them.
#' @returns A named list of [scitargets_dea] objects.
#' @export
normalize_dea_list <- function(x) {
  out <- list()

  add_one <- function(value, name = NULL) {
    if (is_scitargets_dea(value)) {
      nm <- name
      if (is.null(nm) || is.na(nm) || identical(nm, "")) {
        nm <- value@comparison_name
      }
      out[[nm]] <<- value
      return(invisible(NULL))
    }

    if (is.list(value)) {
      value_names <- names(value)
      if (is.null(value_names)) {
        value_names <- rep(NA_character_, length(value))
      }
      for (i in seq_along(value)) {
        add_one(value[[i]], value_names[[i]])
      }
    }

    invisible(NULL)
  }

  add_one(x)
  out
}

# -----------------------------------------------------------------------------
# S7 generics and methods used by Quarto
# -----------------------------------------------------------------------------

#' Differential-expression and enrichment accessors for [scitargets_dea]
#'
#' Generics to extract tables and plots from a [scitargets_dea] object. Most
#' take optional `level` ("single_cell"/"pseudobulk"), `direction` ("up"/"down")
#' and `ontology` ("BP"/"CC"/"MF") arguments, defaulting to the first available.
#'
#' @param x A [scitargets_dea] object.
#' @param ... Passed to methods (e.g. `level`, `direction`, `ontology`).
#' @returns A data.frame, HTML string, or plot object depending on the generic.
#' @name dea_accessors
#' @export
markers_table <- S7::new_generic("markers_table", "x")
#' @rdname dea_accessors
#' @export
go_table <- S7::new_generic("go_table", "x")
#' @rdname dea_accessors
#' @export
go_plot_data <- S7::new_generic("go_plot_data", "x")
#' @rdname dea_accessors
#' @export
go_genes_html <- S7::new_generic("go_genes_html", "x")
#' @rdname dea_accessors
#' @export
volcano_plot <- S7::new_generic("volcano_plot", "x")
#' @rdname dea_accessors
#' @export
pca_plot <- S7::new_generic("pca_plot", "x")
#' @rdname dea_accessors
#' @export
go_barplot <- S7::new_generic("go_barplot", "x")
#' @rdname dea_accessors
#' @export
go_cnetplot <- S7::new_generic("go_cnetplot", "x")
#' Write a [scitargets_dea] object to an xlsx report
#'
#' @param x A [scitargets_dea] object.
#' @param ... Passed to the method, e.g. `out_dir` (output directory) and
#'   `levels` (which differential-expression level(s) to write).
#' @returns The path to the written xlsx file (invisibly via the method).
#' @export
dea_write_xlsx <- S7::new_generic("dea_write_xlsx", "x")

S7::method(markers_table, scitargets_dea) <- function(x, level = NULL, ...) {
  level <- .dea_default_level(x, level)
  tab <- x@de[[level]]
  if (is.null(tab)) data.frame() else tab
}

S7::method(go_table, scitargets_dea) <- function(x,
                                                 direction = c("up", "down", "all"),
                                                 ontology = NULL,
                                                 level = NULL,
                                                 ...) {
  direction <- match.arg(direction)
  res <- .dea_go_result(x, level, direction, ontology)
  if (is.null(res)) data.frame() else res$res.table
}

S7::method(go_plot_data, scitargets_dea) <- function(x,
                                                     direction = c("up", "down", "all"),
                                                     ontology = NULL,
                                                     level = NULL,
                                                     ...) {
  direction <- match.arg(direction)
  res <- .dea_go_result(x, level, direction, ontology)
  if (is.null(res)) data.frame() else res$plot_data
}

S7::method(go_genes_html, scitargets_dea) <- function(x,
                                                      direction = c("up", "down", "all"),
                                                      ontology = NULL,
                                                      level = NULL,
                                                      ...) {
  direction <- match.arg(direction)
  res <- .dea_go_result(x, level, direction, ontology)
  html <- if (is.null(res)) NULL else res$GenesAssociatedToGO$HTML
  if (is.null(html) || is.na(html) || identical(html, "")) {
    html <- "<p>No significant GO terms were detected.</p>"
  }
  html
}

.no_result_plot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 5) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()
}

S7::method(volcano_plot, scitargets_dea) <- function(x,
                                                     level = NULL,
                                                     alpha = NULL,
                                                     interactive = FALSE,
                                                     width_svg = 12,
                                                     height_svg = 7,
                                                     ...) {
  level <- .dea_default_level(x, level)
  # Default the significance cutoff to the per-level adjusted-p cutoff stored on
  # the object (P1.8); an explicit `alpha` still overrides.
  if (is.null(alpha)) {
    alpha <- x@padj_cutoffs[[level]] %||% 0.05
  }
  plot_df <- markers_table(x, level = level)

  if (x@status != "computed" || nrow(plot_df) == 0L) {
    return(.no_result_plot(x@reason %||% "No marker result available."))
  }

  logfc_col <- .logfc_column(plot_df)
  plot_df$significant <- !is.na(plot_df$p_val_adj) & plot_df$p_val_adj < alpha
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$p_val_adj, .Machine$double.xmin))
  plot_df$tooltip <- paste0(
    "Gene: ", plot_df$Gene,
    "\n", logfc_col, ": ", signif(plot_df[[logfc_col]], 3),
    "\nadj. p-value: ", signif(plot_df$p_val_adj, 3)
  )

  base <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[logfc_col]], y = neg_log10_padj)
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = -log10(alpha),
      label = paste0("adj. p = ", alpha),
      hjust = 1.05,
      vjust = -0.5,
      size = 4
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
      name = paste0("p_val_adj < ", alpha)
    ) +
    ggplot2::labs(
      x = "Average log2 fold change",
      y = expression(-log[10](adjusted ~ italic(p))),
      title = "Volcano plot",
      subtitle = paste0(
        "Group 1: ", x@group1,
        " vs Group 2: ", x@group2,
        " | Cluster: ", x@cluster
      )
    ) +
    ggplot2::theme_classic(base_size = 14)

  # ggiraph: gene-name tooltips on hover. ggplot2 (default): static points with
  # ggrepel labels on the SIGNIFICANT genes only.
  if (isTRUE(interactive)) {
    p <- base + ggiraph::geom_point_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = Gene, color = significant),
      alpha = 0.75, size = 2
    )
    return(ggiraph::girafe(
      ggobj = p,
      width_svg = width_svg,
      height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "fill:orange;stroke:black;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  base +
    ggplot2::geom_point(ggplot2::aes(color = significant), alpha = 0.75, size = 2) +
    ggrepel::geom_text_repel(
      data = plot_df[plot_df$significant, , drop = FALSE],
      ggplot2::aes(label = Gene), size = 3, max.overlaps = 20, show.legend = FALSE
    )
}

# Pseudobulk-sample PCA (DESeq2 VST). Colour = comparison role (group1/group2/
# "other"); SHAPE = outlier flag (mt::pca.outlier). `interactive = TRUE` returns a
# ggiraph girafe with per-sample tooltips; otherwise a static ggplot with
# ggrepel-labelled sample ids.
S7::method(pca_plot, scitargets_dea) <- function(x,
                                                 level = "pseudobulk",
                                                 interactive = FALSE,
                                                 width_svg = 8,
                                                 height_svg = 6,
                                                 ...) {
  pca <- x@pca[[level]]
  if (is.null(pca) || !identical(pca$status, "computed") || is.null(pca$coords)) {
    return(.no_result_plot(
      (if (is.null(pca)) NULL else pca$reason) %||% "No pseudobulk PCA available."
    ))
  }
  df <- pca$coords
  df$unit <- sub("^g", "", df$unit)
  ve <- pca$var_explained
  xlab <- sprintf("PC1 (%.1f%%)", 100 * ve[["PC1"]])
  ylab <- if (is.finite(ve[["PC2"]])) sprintf("PC2 (%.1f%%)", 100 * ve[["PC2"]]) else "PC2"
  df$Outlier <- ifelse(df$outlier, "outlier", "kept")
  df$tooltip <- paste0(
    "Sample: ", df$unit, "\nGroup: ", df$group,
    "\nrole: ", df$comparison_role,
    "\noutlier: ", ifelse(df$outlier, "yes", "no")
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::scale_shape_manual(values = c(kept = 16, outlier = 17), name = "Outlier") +
    ggplot2::labs(
      x = xlab, y = ylab, color = "Role", shape = "Outlier",
      title = "Pseudobulk PCA (DESeq2 VST)",
      subtitle = paste0(
        "Group 1: ", x@group1, " vs Group 2: ", x@group2,
        " | Cluster: ", x@cluster
      )
    ) +
    ggplot2::theme_classic(base_size = 13)

  if (isTRUE(interactive)) {
    p <- p + ggiraph::geom_point_interactive(
      ggplot2::aes(
        color = comparison_role, shape = Outlier,
        tooltip = tooltip, data_id = sample
      ),
      size = 3, alpha = 0.85
    )
    return(ggiraph::girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "fill:orange;stroke:black;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  p +
    ggplot2::geom_point(
      ggplot2::aes(color = comparison_role, shape = Outlier),
      size = 3, alpha = 0.85
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = unit, color = comparison_role),
      size = 3, max.overlaps = 20, show.legend = FALSE
    )
}

S7::method(go_barplot, scitargets_dea) <- function(x,
                                                   direction = c("up", "down", "all"),
                                                   ontology = NULL,
                                                   level = NULL,
                                                   interactive = FALSE,
                                                   width_svg = 10,
                                                   height_svg = 7,
                                                   ...) {
  direction <- match.arg(direction)
  res <- .dea_go_result(x, level, direction, ontology)
  gop <- if (is.null(res)) NULL else res$plot_data
  onto_label <- .dea_default_ontology(x, ontology)

  if (is.null(gop) || nrow(gop) == 0L) {
    return(.no_result_plot(
      (if (is.null(res)) NULL else res$reason) %||% "No significant GO terms were detected."
    ))
  }

  cap <- paste0(
    "background: ", length(res$background_genes),
    " ; foreground: ", length(res$foreground_genes)
  )
  title <- if (identical(direction, "all")) {
    paste0("GO:", onto_label, " enrichment in all DE genes")
  } else {
    paste0("GO:", onto_label, " enrichment in ", direction, "-regulated genes")
  }

  # ggiraph: per-bar tooltip with the term, raw + adjusted p-value, enrichment
  # score and the foreground genes annotated to the term.
  if (isTRUE(interactive)) {
    pv <- if ("pval" %in% names(gop)) suppressWarnings(as.numeric(gop$pval)) else NA_real_
    gop$tooltip <- paste0(
      "Term: ", gop$Term,
      "\nraw p-value: ", signif(pv, 3),
      "\nadj. p-value: ", signif(gop$adjpval, 3),
      "\nEnrichment: ", gop$Enrichment,
      "\nForeground genes: ", gop$Genes
    )
    p <- ggplot2::ggplot(gop, ggplot2::aes(x = Term_wrapped, y = neg_log10_adjpval)) +
      ggiraph::geom_col_interactive(
        ggplot2::aes(tooltip = tooltip, data_id = GO.ID),
        fill = "cornflowerblue"
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "", y = "- log10 adjusted P-value", title = title, caption = cap) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 12),
        panel.background = ggplot2::element_blank()
      )
    return(ggiraph::girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "fill:orange;stroke:black;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  p <- ggplot2::ggplot(
    data = gop,
    ggplot2::aes(x = Term_wrapped, y = neg_log10_adjpval)
  ) +
    ggplot2::geom_bar(stat = "identity", fill = "cornflowerblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "", y = "- log10 adjusted P-value") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 6),
      axis.text.y = ggplot2::element_text(size = 16),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Enrichment),
      position = ggplot2::position_dodge(width = 0),
      vjust = 0.5,
      hjust = 1.1
    ) +
    ggplot2::ggtitle(title)

  ggpubr::annotate_figure(p, bottom = cap)
}

S7::method(go_cnetplot, scitargets_dea) <- function(x,
                                                    direction = c("up", "down", "all"),
                                                    ontology = NULL,
                                                    level = NULL,
                                                    interactive = TRUE,
                                                    width_svg = 10,
                                                    height_svg = 7,
                                                    ...) {
  direction <- match.arg(direction)
  res <- .dea_go_result(x, level, direction, ontology)
  cnet <- if (is.null(res)) NULL else res$cnetplot.data

  if (is.null(cnet) || is.null(cnet$edges) || is.null(cnet$nodes)) {
    return(.no_result_plot(
      (if (is.null(res)) NULL else res$reason) %||% "No GO-gene network is available."
    ))
  }

  edge_df <- cnet$edges
  cnet_nodes <- cnet$nodes

  p_cnet <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey75",
      alpha = 0.5,
      linewidth = 0.4
    ) +
    ggiraph::geom_point_interactive(
      data = cnet_nodes,
      ggplot2::aes(
        x = xpos,
        y = ypos,
        color = node_type,
        size = node_size,
        tooltip = tooltip,
        data_id = name
      ),
      alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      data = subset(cnet_nodes, node_type == "GO term"),
      ggplot2::aes(
        x = xpos,
        y = ypos,
        label = short_label,
        color = node_type
      ),
      size = 3,
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("Gene" = "#3B82F6", "GO term" = "#EF4444")
    ) +
    ggplot2::scale_size_continuous(range = c(2.5, 8), guide = "none") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = NULL)

  if (!isTRUE(interactive)) {
    return(p_cnet)
  }

  out <- ggiraph::girafe(
    ggobj = p_cnet,
    width_svg = width_svg,
    height_svg = height_svg
  )

  ggiraph::girafe_options(
    out,
    ggiraph::opts_hover(css = "stroke:black;stroke-width:2px;"),
    ggiraph::opts_toolbar(saveaspng = TRUE),
    ggiraph::opts_selection(type = "none"),
    ggiraph::opts_tooltip(
      css = "background-color:white;color:black;padding:8px;border:1px solid #999;border-radius:4px;"
    )
  )
}

# S7 method: write the report (summary + per-level markers + per-ontology GO +
# GSEA). `levels` restricts which differential-expression level(s) are written
# (NULL = all computed); the summary sheet is always written.
S7::method(dea_write_xlsx, scitargets_dea) <- function(x, out_dir = "./out/clinical_markers",
                                                       levels = NULL, ...) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  path <- file.path(out_dir, .dea_xlsx_filename(x))
  wb <- openxlsx2::wb_workbook()

  onto <- .dea_ontologies(x)
  write_levels <- if (is.null(levels)) .dea_levels(x) else intersect(levels, .dea_levels(x))

  .write_summary_sheet(wb, x, shown_levels = write_levels)

  # Prefix the per-level sheet names only when several levels share one workbook
  # (avoids sc_markers / pb_markers collisions). A per-level file -- the
  # write_dea_xlsx() default -- needs no prefix.
  level_prefix <- c(single_cell = "sc", pseudobulk = "pb")
  multi_level <- length(write_levels) > 1L

  for (level in write_levels) {
    prefix <- if (multi_level) level_prefix[[level]] %||% level else NULL

    .write_table_sheet(
      wb,
      sheet = .sheet_name(prefix, "markers"),
      data = markers_table(x, level = level)
    )

    for (direction in c("up", "down", "all")) {
      for (o in onto) {
        res <- .dea_go_result(x, level = level, direction = direction, ontology = o)
        reason <- if (is.null(res)) NA_character_ else res$reason
        .write_table_sheet(
          wb,
          sheet = .sheet_name(prefix, "GO", o, direction),
          data = go_table(x, direction = direction, ontology = o, level = level),
          empty_data = .go_empty_message(reason, direction, o)
        )
      }
    }

    # GSEA: one sheet per collection for this level.
    for (collection in .dea_gsea_collections(x, level = level)) {
      gres <- .dea_gsea_result(x, level = level, collection = collection)
      greason <- if (is.null(gres)) NA_character_ else gres$reason
      .write_table_sheet(
        wb,
        sheet = .sheet_name(prefix, "GSEA", collection),
        data = gsea_table(x, collection = collection, level = level),
        empty_data = data.frame(
          message = if (is.null(greason) || is.na(greason)) {
            paste0("No GSEA result for ", collection, ".")
          } else {
            greason
          },
          stringsAsFactors = FALSE
        )
      )
    }
  }

  openxlsx2::wb_save(wb, path, overwrite = TRUE)
  path
}

#' Write a [scitargets_dea] branch to per-level xlsx reports
#'
#' Convenience wrapper for the targets pipeline: a mapped branch returns a named
#' list of length one, so the single object is extracted before writing. Writes
#' one xlsx per differential-expression level into a level subfolder of
#' `out_dir` (`out_dir/single_cell/...xlsx` and `out_dir/pseudobulk/...xlsx`),
#' so single-cell and pseudobulk results are kept separate.
#'
#' @param x A [scitargets_dea] object, or a length-one list containing one.
#' @param out_dir Output directory; per-level subfolders are created under it.
#' @param levels Which levels to write (default both single_cell and pseudobulk).
#' @returns The paths to the written xlsx files.
#' @export
write_dea_xlsx <- function(x, out_dir = "./out/clinical_markers",
                           levels = c("single_cell", "pseudobulk")) {
  x <- .extract_single_dea(x)
  vapply(
    levels,
    function(lv) dea_write_xlsx(x, out_dir = file.path(out_dir, lv), levels = lv),
    character(1)
  )
}

# =============================================================================
# Comparison selection + pseudobulk (DESeq2) level
# -----------------------------------------------------------------------------
# These helpers are wired together by the single-entry dispatcher `run_dea`
# (defined at the end of this file), which builds one scitargets_dea object per
# comparison from the requested differential-expression level(s).
# =============================================================================

# Resolve which group pairs to compare, from the `groups` specification.
#   * NA (default)        -> all pairwise comparisons of `group_levels`
#   * two valid groups    -> exactly that pair (group1 vs group2)
#   * one valid group     -> that group vs "rest" (all other cells)
# Invalid group names raise an error.
resolve_group_comparisons <- function(group_levels, groups = NA) {
  group_levels <- unique(as.character(stats::na.omit(group_levels)))

  if (length(groups) == 1L && (is.na(groups) || identical(groups, NA))) {
    if (length(group_levels) < 2L) {
      return(list())
    }
    return(utils::combn(group_levels, 2L, simplify = FALSE))
  }

  groups <- as.character(groups)
  invalid <- setdiff(groups, group_levels)
  if (length(invalid) > 0L) {
    stop(
      "Invalid group name(s): ", paste(invalid, collapse = ", "),
      ". Available groups: ", paste(group_levels, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(groups) == 2L) {
    return(list(as.character(groups)))
  }
  if (length(groups) == 1L) {
    # one group versus all the others
    return(list(c(groups, "rest")))
  }

  stop(
    "`groups` must be NA (all pairwise), a single group (one-vs-rest), ",
    "or exactly two groups.",
    call. = FALSE
  )
}


#' Enumerate clinical-group comparisons to run
#'
#' Returns one row per (cluster x group-pair) comparison that has both groups
#' present (with at least `min_cells_per_group` cells) inside the cluster.
#'
#' @param meta A metadata data.frame (e.g. `seurat_obj@meta.data`).
#' @param group_by Name of the column holding the group labels.
#' @param cluster_by Name of the column holding the clustering (must differ from
#'   `group_by`).
#' @param groups Group selection: `NA` (default) for all pairwise comparisons,
#'   two group names for a single pair, or one group name for one-vs-rest.
#' @param clusters Which clusters to use (`NA` = all). Use this to restrict the
#'   analysis, e.g. keep only CD3+ T-cell Azimuth celltypes and drop Platelet /
#'   Eryth / B cells.
#' @param min_cells_per_group Minimum cells per group within a cluster.
#' @returns A data.frame with columns comparison_name, group1, group2, cluster,
#'   group_by, cluster_by.
#' @export
dea_comparisons <- function(meta,
                            group_by,
                            cluster_by,
                            groups = NA,
                            clusters = NA,
                            min_cells_per_group = 3L) {
  stopifnot(group_by %in% colnames(meta))
  if (!identical(group_by, cluster_by) && !cluster_by %in% colnames(meta)) {
    stop("cluster_by column '", cluster_by, "' not found in metadata.", call. = FALSE)
  }
  if (identical(group_by, cluster_by)) {
    stop("group_by and cluster_by must be different columns.", call. = FALSE)
  }

  meta <- as.data.frame(meta)
  cluster_levels <- .resolve_clusters(meta[[cluster_by]], clusters)
  records <- list()

  for (cl in cluster_levels) {
    in_cluster <- meta[as.character(meta[[cluster_by]]) == cl & !is.na(meta[[cluster_by]]), , drop = FALSE]
    group_counts <- table(as.character(in_cluster[[group_by]]))
    present_groups <- names(group_counts)[group_counts >= min_cells_per_group]

    if (length(groups) == 1L && (is.na(groups) || identical(groups, NA))) {
      pairs <- resolve_group_comparisons(present_groups, NA)
    } else {
      # validate against the full set of group levels, then keep the pairs
      # whose required groups are present in this cluster.
      all_levels <- unique(as.character(stats::na.omit(meta[[group_by]])))
      pairs <- resolve_group_comparisons(all_levels, groups)
      pairs <- Filter(function(p) {
        needed <- p[p != "rest"]
        all(needed %in% present_groups)
      }, pairs)
    }

    for (p in pairs) {
      g1 <- p[[1L]]
      g2 <- p[[2L]]
      records[[length(records) + 1L]] <- data.frame(
        comparison_name = paste(g1, g2, cl, sep = "::"),
        group1 = g1,
        group2 = g2,
        cluster = cl,
        group_by = group_by,
        cluster_by = cluster_by,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(records) == 0L) {
    return(data.frame(
      comparison_name = character(), group1 = character(), group2 = character(),
      cluster = character(), group_by = character(), cluster_by = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, records)
}

# Run DESeq2 on a pseudobulk count matrix.
#   counts   : genes x pseudobulk-samples integer matrix
#   col_data : data.frame with rownames == colnames(counts) and a `group` factor
#   group1   : the level whose up-regulation gives a positive log2 fold change
#   group2   : the reference level (control)
#   test     : DESeq2 test, "Wald" (default) or "LRT".
#   design   : length-1 character design formula passed through as.formula() to
#              DESeq2 (default "~ group"). MUST reference the `group` factor so the
#              group1-vs-group2 contrast can be extracted.
#   reduced  : length-1 character reduced formula (default "~ 1"); only used when
#              test == "LRT" (DESeq2 requires it for the likelihood-ratio test).
#   low_count_filter : length-2 integer c(count, n_samples); keep genes whose count
#              is >= count in at least n_samples pseudobulk samples (default
#              c(10, 3); supersedes the old "any count > 0" filter).
#   alpha    : the downstream adjusted-p (FDR) cutoff, passed to DESeq2::results()
#              so its independent filtering is optimised for that threshold (DESeq2
#              default 0.1; run_dea passes padj_cutoff_pseudobulk).
# The group labels are sanitized for the DESeq2 design only; the returned table
# has no group names (just per-gene stats), so the original labels are unaffected.
# The number of genes before/after the low-count filter is recorded on the result
# as attributes "n_genes_prefilter" / "n_genes_postfilter".
run_deseq2 <- function(counts, col_data, group1, group2,
                       test = c("Wald", "LRT"),
                       design = "~ group",
                       reduced = "~ 1",
                       low_count_filter = c(10L, 3L),
                       shrinkage = TRUE,
                       alpha = 0.1) {
  test <- match.arg(test)
  design_formula <- stats::as.formula(design)

  # Check for reduced formula if LRT is used
  if (identical(test, "LRT")) {
    reduced_formula <- stats::as.formula(reduced)
    full_terms <- attr(stats::terms(design_formula), "term.labels")
    reduced_terms <- attr(stats::terms(reduced_formula), "term.labels")

    if (!"group" %in% full_terms) {
      stop(
        "For pb_test = 'LRT', the full design must contain `group` as a main term.",
        call. = FALSE
      )
    }

    if ("group" %in% reduced_terms) {
      stop(
        "For pb_test = 'LRT', the reduced formula must not contain `group`.",
        call. = FALSE
      )
    }

    non_group_terms <- setdiff(full_terms, "group")
    missing_reduced_terms <- setdiff(non_group_terms, reduced_terms)
    if (length(missing_reduced_terms) > 0L) {
      stop(
        "For pb_test = 'LRT', the reduced formula should retain all non-group terms ",
        "from the full design. Missing from reduced formula: ",
        paste(missing_reduced_terms, collapse = ", "),
        ". For example, if design = '~ batch + sex + group', use reduced = '~ batch + sex'.",
        call. = FALSE
      )
    }

    group_interaction_terms <- grep("(^|:)group(:|$)", full_terms, value = TRUE)
    if (length(group_interaction_terms) > 0L) {
      stop(
        "pb_test = 'LRT' with group interactions is not supported in this pairwise ",
        "GO/GSEA pipeline because the LRT p-value is not a simple signed group1-vs-group2 test. ",
        "Problematic term(s): ",
        paste(group_interaction_terms, collapse = ", "),
        call. = FALSE
      )
    }

    message("LRT test was set for the pseudobulk DE analysis. Scitargets will compute GO enrichment only (GSEA is skipped for LRT, as the LRT statistic is unsigned), but be aware that GO may over-interpret the enrichment! We recommand to use the default 'Wald' test.")
  }

  if (!"group" %in% all.vars(design_formula)) {
    stop("The DESeq2 design formula must reference the `group` term (got '",
      design, "').",
      call. = FALSE
    )
  }
  low_count_filter <- as.integer(low_count_filter)
  if (length(low_count_filter) != 2L || anyNA(low_count_filter) ||
    any(low_count_filter < 0L)) {
    stop("low_count_filter must be two non-negative integers c(count, n_samples).",
      call. = FALSE
    )
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
    alpha <= 0 || alpha > 1) {
    stop("alpha must be a single number in (0, 1].", call. = FALSE)
  }

  col_data <- as.data.frame(col_data)
  s1 <- .dea_safe_level(group1)
  s2 <- .dea_safe_level(group2)
  col_data$group <- factor(.dea_safe_level(col_data$group), levels = c(s2, s1))

  if (!"sample" %in% colnames(col_data)) {
    stop("col_data must contain a 'sample' column matching colnames(counts).", call. = FALSE)
  }
  if (!identical(as.character(col_data$sample), colnames(counts))) {
    stop("col_data$sample must be in the same order as colnames(counts).", call. = FALSE)
  }

  rownames(col_data) <- col_data$sample

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = col_data,
    design = design_formula
  )
  # Low-count gene filter (edgeR::filterByExpr philosophy): keep genes with
  # count >= low_count_filter[1] in at least low_count_filter[2] samples.
  n_prefilter <- nrow(dds)
  keep_gene <- rowSums(DESeq2::counts(dds) >= low_count_filter[1L]) >= low_count_filter[2L]
  dds <- dds[keep_gene, ]
  n_postfilter <- nrow(dds)

  if (identical(test, "LRT")) {
    dds <- DESeq2::DESeq(dds,
      test = "LRT",
      reduced = stats::as.formula(reduced), quiet = TRUE
    )
  } else {
    dds <- DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)
  }
  # `alpha` must match the downstream FDR cutoff so DESeq2's independent filtering
  # is optimised for it. See https://www.biostars.org/p/9551104/
  res <- DESeq2::results(dds, contrast = c("group", s1, s2), alpha = alpha)
  out <- .standardize_deseq2(res)

  # Per-group MEAN DESeq2-normalized count columns (named with the ORIGINAL group
  # labels), aligned to the standardized table by gene name.
  nc <- DESeq2::counts(dds, normalized = TRUE)
  g <- as.character(col_data$group)
  mean1 <- rowMeans(nc[, g == s1, drop = FALSE])
  mean2 <- rowMeans(nc[, g == s2, drop = FALSE])
  out[[paste0("mean_norm_count_", group1)]] <- unname(mean1[out$Gene])
  out[[paste0("mean_norm_count_", group2)]] <- unname(mean2[out$Gene])

  # Shrunken log2 fold change reported ALONGSIDE the MLE. avg_log2FC stays the MLE
  # so the sign-based GO/GSEA logic is unchanged; the shrunken estimate goes in
  # avg_log2FC_shrink (NA when shrinkage is disabled). apeglm preferred, else ashr,
  # else DESeq2's built-in normal prior.
  if (isTRUE(shrinkage)) {
    shrink <- .dea_lfc_shrink(dds, s1, s2)
  } else {
    shrink <- list(
      lfc = stats::setNames(rep(NA_real_, nrow(dds)), rownames(dds)),
      method = "disabled"
    )
  }
  out$avg_log2FC_shrink <- unname(shrink$lfc[out$Gene])

  attr(out, "n_genes_prefilter") <- n_prefilter
  attr(out, "n_genes_postfilter") <- n_postfilter
  attr(out, "deseq2_test") <- test
  attr(out, "lfcShrink_method") <- shrink$method
  attr(out, "n_pb_samples") <- ncol(counts)
  out
}


# -----------------------------------------------------------------------------
# Single-entry dispatcher
# -----------------------------------------------------------------------------

#' Run a single clinical-group differential-expression comparison
#'
#' Builds one [scitargets_dea] object for a comparison (group1 vs group2, within
#' `cluster`), running the requested differential-expression level(s) and the
#' downstream GO and GSEA enrichment.
#'
#' @param seurat_obj A Seurat object.
#' @param group_by Metadata column holding the group labels.
#' @param cluster_by Metadata column holding the clustering (must differ from
#'   `group_by`); pass `character()` / `NA` to compare across all cells.
#' @param cluster The cluster value this comparison sits in.
#' @param group1,group2 The two groups; `group2 = "rest"` means one-vs-rest.
#' @param level One of "both" (default), "single_cell", "pseudobulk".
#' @param run_go,run_gsea Whether to run the enrichment analyses (default TRUE).
#'   GSEA needs a valid signed gene ranking, so it is computed only where that
#'   holds (a warning is emitted, and the level's GSEA tab is omitted from the
#'   report, otherwise):
#'   - **single cell**: only when the `Seurat::FindMarkers` gene universe is
#'     unfiltered, i.e. `sc_logfc_threshold = 0`, `sc_min_pct = 0` and
#'     `sc_only_pos = FALSE`; genes are ranked by `sign(log2FC) * -log10(p_val)`.
#'   - **pseudobulk**: only when `pb_test = "Wald"`; genes are ranked by the signed
#'     DESeq2 Wald statistic (the `stat` column). It is skipped for `pb_test = "LRT"`
#'     (the LRT statistic is unsigned).
#' @param go_params,gsea_params Named lists of GO / GSEA parameters. Pass
#'   `gsea_params$pathways` (from [get_msigdbr_pathways()]) to skip re-querying.
#' @param pseudobulk_unit Biological replicate column for DESeq2 (default
#'   "patient_id").
#' @param min_cells_per_group,min_replicates,min_cells_per_sample Quality
#'   thresholds for the single-cell and pseudobulk levels. `min_replicates`
#'   (pseudobulk only; default 3) is the minimum number of biological replicates
#'   (`pseudobulk_unit`) required per group.
#' @param sc_min_pct Single-cell level only: minimum fraction of cells (in
#'   either group) expressing a gene for it to be tested, passed to
#'   [Seurat::FindMarkers()] as `min.pct` (default 0.1).
#' @param sc_logfc_threshold Single-cell level only: minimum absolute log2 fold
#'   change for a gene to be tested, passed to [Seurat::FindMarkers()] as
#'   `logfc.threshold` (default 0.25).
#' @param sc_test_use Single-cell level only: the differential-expression test
#'   used by [Seurat::FindMarkers()] (its `test.use`, default "wilcox").
#' @param sc_only_pos Single-cell level only: whether to return only
#'   positively-regulated (up in `group1`) genes, passed to
#'   [Seurat::FindMarkers()] as `only.pos` (default FALSE = both directions).
#'   Leave FALSE for the standard two-sided analysis; setting it TRUE drops the
#'   down-regulated genes, so the volcano and the "down" GO/GSEA outputs are
#'   empty for this level.
#' @param pb_test Pseudobulk DESeq2 test: "Wald" (default) or "LRT". Pseudobulk level only.
#'   Use `"Wald"` for the standard pairwise group1-vs-group2 analysis and for downstream signed GO/GSEA.
#'   "LRT"` is intended for testing the contribution of `group` relative to a reduced model; when using `"LRT"`, `pb_reduced` should retain all non-group covariates from `pb_design`.
#' @param pb_design Pseudobulk DESeq2 design as a length-1 character formula
#'   (default "~ group"), passed through [stats::as.formula()] to DESeq2. Must
#'   reference the `group` term so the group1-vs-group2 contrast can be extracted.
#'   Pseudobulk level only.
#' @param pb_reduced Reduced formula (length-1 character, default "~ 1") used only
#'   when `pb_test = "LRT"`. Pseudobulk level only.
#' @param pb_low_count_filter Length-2 integer `c(count, n_samples)`: keep genes
#'   whose count is `>= count` in at least `n_samples` pseudobulk samples (default
#'   `c(10, 3)`). Pseudobulk level only; does not affect the single-cell level.
#'   Covariates referenced by `pb_design` (e.g. `~ batch + group`) are looked up
#'   one value per biological replicate, keyed on `pseudobulk_unit` (and an error is
#'   raised if a covariate is not constant within a unit).
#' @param pb_lfc_shrink Whether to report a shrunken log2 fold change for the
#'   pseudobulk level (default TRUE), added as the `avg_log2FC_shrink` column
#'   alongside the MLE `avg_log2FC`. Method cascade: apeglm (if installed) ->
#'   ashr (if installed) -> DESeq2 built-in normal prior; `run_dea` prints which
#'   method was used. Pseudobulk level only; the single-cell level is unaffected.
#' @param pb_pca Whether to compute a pseudobulk-sample PCA + outlier detection
#'   over ALL samples in the cluster (default TRUE). Uses DESeq2 VST + the
#'   top-variance genes, with outliers flagged by [mt::pca.outlier]. Computed
#'   independently of the DE replicate gate (a QC view) and stored on `@pca`.
#'   Pseudobulk level only.
#' @param pca_n_top_genes Number of top-variance VST genes used for the pseudobulk
#'   PCA (default 500).
#' @param pca_outlier_conf Confidence level passed to [mt::pca.outlier] for the
#'   pseudobulk-sample outlier flag (default 0.975).
#' @param pb_remove_outliers Whether to REMOVE the PCA-flagged outlier replicates
#'   from the pseudobulk DESeq2 computation (default FALSE = keep them). When TRUE,
#'   outlier samples are dropped before the replicate gate, so the
#'   `min_replicates` / `min_cells_per_sample` cutoffs are re-checked on the
#'   post-removal sample set (the comparison becomes not-computed if too few
#'   replicates remain). The `@pca` view always retains every sample with its
#'   outlier flag. Pseudobulk level only.
#' @param assay_sc Assay used for the single-cell FindMarkers test.
#' @param recorrect_umi Passed to [Seurat::FindMarkers()] for the single-cell
#'   level. `run_dea()` does NOT run [Seurat::PrepSCTFindMarkers()] itself and no
#'   longer subsets the object (it selects the comparison via `group.by`), so the
#'   object should be `PrepSCTFindMarkers`-corrected upstream once and this left
#'   at its default `FALSE` (reuse the already-corrected SCT counts).
#' @param species Organism, `"human"` (default) or `"mouse"`. Sets the GO
#'   annotation database (org.Hs.eg.db / org.Mm.eg.db) and the GSEA msigdbr
#'   species; stored on the object.
#' @param padj_cutoff_single_cell,padj_cutoff_pseudobulk Per-level adjusted-p
#'   cutoffs (default 0.05 each). Each drives, for its level, the volcano plot's
#'   horizontal cutoff line and the GO foreground-gene selection. The pseudobulk
#'   cutoff is additionally passed to `DESeq2::results(alpha = )`, so DESeq2's
#'   independent filtering is optimised for that FDR threshold.
#' @returns A [scitargets_dea] object.
#' @export
run_dea <- function(seurat_obj,
                    group_by,
                    cluster_by,
                    cluster,
                    group1,
                    group2,
                    level = c("both", "single_cell", "pseudobulk"),
                    run_go = TRUE,
                    run_gsea = TRUE,
                    go_params = list(),
                    gsea_params = list(),
                    pseudobulk_unit = "patient_id",
                    min_cells_per_group = 3L,
                    min_replicates = 3L,
                    min_cells_per_sample = 10L,
                    sc_min_pct = 0.1,
                    sc_logfc_threshold = 0.25,
                    sc_test_use = "wilcox",
                    sc_only_pos = FALSE,
                    pb_test = c("Wald", "LRT"),
                    pb_design = "~ group",
                    pb_reduced = "~ 1",
                    pb_low_count_filter = c(10L, 3L),
                    pb_lfc_shrink = TRUE,
                    pb_pca = TRUE,
                    pca_n_top_genes = 500L,
                    pca_outlier_conf = 0.975,
                    pb_remove_outliers = FALSE,
                    assay_sc = "SCT",
                    recorrect_umi = FALSE,
                    species = c("human", "mouse"),
                    padj_cutoff_single_cell = 0.05,
                    padj_cutoff_pseudobulk = 0.05) {
  level <- match.arg(level)
  species <- match.arg(species)
  pb_test <- match.arg(pb_test)
  # Per-level adjusted-p cutoffs: drive the volcano horizontal cutoff line and
  # the GO foreground-gene selection for that level (stored on the object).
  padj_cutoffs <- list(
    single_cell = padj_cutoff_single_cell,
    pseudobulk = padj_cutoff_pseudobulk
  )
  levels_req <- if (identical(level, "both")) c("single_cell", "pseudobulk") else level
  params <- .default_go_params(go_params)
  gparams <- .default_gsea_params(gsea_params)
  # Species ("human"/"mouse") is the single source of truth for the GO annotation
  # DB (org.Hs.eg.db / org.Mm.eg.db) and the GSEA msigdbr species; it overrides any
  # organism/species left in go_params/gsea_params.
  sp <- .species_db(species)
  params$organism <- sp$organism
  params$org_db <- sp$org_db
  gparams$species <- sp$msigdbr
  # GSEA plots share the GO `top_nodes` cap (number of terms/sets shown in plots)
  # unless the user set `top_nodes` explicitly in gsea_params -- one knob for both.
  if (is.null(gsea_params$top_nodes)) gparams$top_nodes <- params$top_nodes
  # the (potentially large) pathway list is not stored inside the object
  gparams_store <- gparams
  gparams_store$pathways <- NULL

  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", group_by, "' was not found in seurat_obj@meta.data.", call. = FALSE)
  }
  has_cluster <- length(cluster_by) == 1L && !is.na(cluster_by) && nzchar(cluster_by)
  if (has_cluster) {
    if (identical(group_by, cluster_by)) {
      stop("group_by and cluster_by must be different columns.", call. = FALSE)
    }
    if (!cluster_by %in% colnames(seurat_obj@meta.data)) {
      stop("Column '", cluster_by, "' was not found in seurat_obj@meta.data.", call. = FALSE)
    }
  }

  cluster <- as.character(cluster)
  ref <- if (identical(group2, "rest")) "rest" else group2
  comparison_name <- paste(group1, group2, cluster, sep = "::")
  seed <- tryCatch(targets::tar_seed_get(), error = function(e) 5114L)

  de_params <- list(
    group_by = group_by,
    cluster_by = if (has_cluster) cluster_by else NA_character_,
    cluster = cluster,
    ident.1 = group1,
    ident.2 = ref,
    level = level,
    pseudobulk_unit = pseudobulk_unit,
    min_cells_per_group = min_cells_per_group,
    min_replicates = min_replicates,
    min_cells_per_sample = min_cells_per_sample,
    sc_test = sc_test_use,
    sc_assay = assay_sc,
    sc_min.pct = sc_min_pct,
    sc_logfc.threshold = sc_logfc_threshold,
    sc_only.pos = sc_only_pos,
    sc_recorrect_umi = recorrect_umi,
    pb_method = "DESeq2",
    pb_test = pb_test,
    pb_design = pb_design,
    pb_reduced = pb_reduced,
    pb_low_count_filter = pb_low_count_filter,
    pb_lfc_shrink = pb_lfc_shrink,
    pb_pca = pb_pca,
    pca_n_top_genes = pca_n_top_genes,
    pca_outlier_conf = pca_outlier_conf,
    pb_remove_outliers = pb_remove_outliers
  )

  # Label the two comparison groups on the FULL object (no subsetting -- subset()
  # breaks the shared SCT model state FindMarkers needs). FindMarkers selects via
  # group.by/ident; PseudobulkExpression drops the NA (non-comparison) cells.
  seurat_obj <- SeuratObject::AddMetaData(
    seurat_obj,
    metadata = .dea_group_labels(
      seurat_obj@meta.data, group_by,
      if (has_cluster) cluster_by else character(),
      cluster, group1, group2
    ),
    col.name = ".dea_group"
  )

  grp_n <- table(as.character(seurat_obj$.dea_group))
  n_cells <- data.frame(
    comparison_group = names(grp_n),
    N = as.integer(grp_n),
    stringsAsFactors = FALSE
  )

  # Pseudobulk-sample PCA + outlier detection over the WHOLE cluster (all groups).
  # Computed here -- BEFORE the DE replicate gate -- so it is available even when
  # the DE comparison itself cannot be run. Pseudobulk-only QC; stored on @pca.
  pca <- list()
  if ("pseudobulk" %in% levels_req && isTRUE(pb_pca)) {
    pca$pseudobulk <- tryCatch(
      .de_pseudobulk_pca(
        seurat_obj, group_by, if (has_cluster) cluster_by else character(),
        cluster, group1, group2,
        pseudobulk_unit = pseudobulk_unit,
        min_cells_per_sample = min_cells_per_sample,
        n_top_genes = pca_n_top_genes, outlier_conf = pca_outlier_conf
      ),
      error = function(e) {
        list(
          status = "not_computed",
          reason = paste("PCA error:", conditionMessage(e))
        )
      }
    )
  }

  # Outlier replicates (by unit) to drop from the pseudobulk DE when requested.
  # Identified from the PCA QC; "other"-group outliers are harmless (not in the
  # comparison's aggregation). The @pca view keeps every sample regardless.
  outlier_units <- character()
  if (isTRUE(pb_remove_outliers) &&
    identical(pca$pseudobulk$status, "computed")) {
    co <- pca$pseudobulk$coords
    outlier_units <- unique(co$unit[co$outlier])
  }
  de_params$pb_outlier_units <- outlier_units

  # Guard: both groups must have enough cells overall.
  if (length(grp_n) < 2L ||
    any(is.na(grp_n[c(group1, ref)])) ||
    any(grp_n[c(group1, ref)] < min_cells_per_group)) {
    reason <- paste0(
      "Not computed: at least one group has fewer than ",
      min_cells_per_group, " cells in cluster '", cluster, "'."
    )
    return(scitargets_dea(
      comparison_name = comparison_name,
      group_by = group_by,
      cluster_by = if (has_cluster) cluster_by else character(),
      cluster = cluster,
      groups = c(group1, group2),
      group1 = group1,
      group2 = group2,
      levels = character(),
      pseudobulk_unit = character(),
      n_cells = n_cells,
      de = list(),
      go = stats::setNames(
        lapply(levels_req, function(lv) .empty_go_levels(reason, params)),
        levels_req
      ),
      gsea = list(),
      pca = pca,
      de_params = de_params,
      go_params = params,
      gsea_params = gparams_store,
      species = species,
      padj_cutoffs = padj_cutoffs,
      status = "not_computed",
      reason = reason
    ))
  }

  # GSEA is only computed where the gene ranking is valid:
  #   - single cell: needs the FULL unfiltered gene ranking, so only when
  #     sc_logfc_threshold == 0, sc_min_pct == 0 and sc_only_pos == FALSE.
  #   - pseudobulk: only for the Wald test (ranked by the signed DESeq2 Wald
  #     statistic); the LRT statistic is unsigned, so pseudobulk GSEA is skipped
  #     for pb_test = "LRT".
  sc_gsea_ok <- (sc_logfc_threshold == 0) && (sc_min_pct == 0) && !isTRUE(sc_only_pos)
  pb_gsea_ok <- identical(pb_test, "Wald")
  # Gene-set collections for GSEA. Fetch once, and only if some level will use it;
  # callers can pass a pre-built list via gsea_params$pathways to skip msigdbr.
  gsea_pathways <- NULL
  if (isTRUE(run_gsea) &&
    (("single_cell" %in% levels_req && sc_gsea_ok) ||
      ("pseudobulk" %in% levels_req && pb_gsea_ok))) {
    gsea_pathways <- gsea_params$pathways %||%
      get_msigdbr_pathways(gparams$collections, species = gparams$species)
  }

  de <- list()
  go <- list()
  gsea <- list()
  level_reasons <- list()

  for (lv in levels_req) {
    markers <- tryCatch(
      {
        if (identical(lv, "single_cell")) {
          # No subsetting: FindMarkers selects the two groups via group.by +
          # ident on the FULL (PrepSCTFindMarkers'd, upstream) object, so the SCT
          # model state stays intact. `recorrect_umi` (run_dea arg) is forwarded.
          .de_single_cell(seurat_obj, ".dea_group",
            ident.1 = group1, ident.2 = ref,
            min_pct = sc_min_pct, logfc_threshold = sc_logfc_threshold,
            test_use = sc_test_use, only_pos = sc_only_pos,
            assay = assay_sc, seed = seed,
            recorrect_umi = recorrect_umi
          )
        } else {
          # PseudobulkExpression aggregates by (unit x .dea_group) and drops the
          # NA (non-comparison) cells, so it runs on the full object directly.
          .de_pseudobulk(seurat_obj,
            group1 = group1, group2 = group2,
            pseudobulk_unit = pseudobulk_unit,
            min_replicates = min_replicates,
            min_cells_per_sample = min_cells_per_sample,
            test = pb_test, design = pb_design,
            reduced = pb_reduced,
            low_count_filter = pb_low_count_filter,
            shrinkage = pb_lfc_shrink,
            # DESeq2 independent filtering is optimised at this FDR cutoff; use the
            # object's per-level pseudobulk cutoff (stored on @padj_cutoffs).
            alpha = padj_cutoffs[["pseudobulk"]],
            drop_units = outlier_units
          )
        }
      },
      error = function(e) .dea_failed(conditionMessage(e))
    )

    if (is.null(markers) || .dea_is_failed(markers) || nrow(markers) == 0L) {
      reason_lv <- if (.dea_is_failed(markers)) {
        attr(markers, "reason")
      } else {
        paste0(lv, " level produced no result.")
      }
      level_reasons[[lv]] <- reason_lv
      go[[lv]] <- .empty_go_levels(reason_lv, params)
      next
    }

    de[[lv]] <- markers
    if (identical(lv, "pseudobulk")) {
      .dea_message_shrink(
        comparison_name, attr(markers, "lfcShrink_method"),
        pb_lfc_shrink
      )
    }
    # Task 6: when DE genes exist, launch both topGO GO and GSEA.
    if (isTRUE(run_go)) {
      go[[lv]] <- .go_for_markers(markers, params, fg_cutoff = padj_cutoffs[[lv]])
    } else {
      go[[lv]] <- .empty_go_levels("GO not requested (run_go = FALSE).", params)
    }
    if (isTRUE(run_gsea)) {
      if (identical(lv, "single_cell")) {
        if (sc_gsea_ok) {
          gsea[[lv]] <- compute_gsea(markers, gsea_pathways, gparams)
        } else {
          warning(
            "Skipping single-cell GSEA: it needs the full unfiltered gene ranking. ",
            "Set sc_logfc_threshold = 0, sc_min_pct = 0 and sc_only_pos = FALSE.",
            call. = FALSE
          )
        }
      } else if (identical(lv, "pseudobulk")) {
        if (pb_gsea_ok) {
          gsea[[lv]] <- compute_gsea(markers, gsea_pathways, gparams)
        } else {
          warning(
            "Skipping pseudobulk GSEA for pb_test = 'LRT'; use Wald for signed pairwise GSEA.",
            call. = FALSE
          )
        }
      }
    }
  }

  # Record the LFC-shrinkage method actually used for pseudobulk (apeglm/ashr/
  # normal/none/disabled), so the report + xlsx can document it.
  if (!is.null(de[["pseudobulk"]])) {
    de_params$pb_lfc_shrink_method <- attr(de[["pseudobulk"]], "lfcShrink_method")
  }

  computed_levels <- if (length(de) > 0L) names(de) else character()
  status <- if (length(computed_levels) > 0L) "computed" else "not_computed"
  reason <- if (length(level_reasons) > 0L) {
    paste(
      sprintf("%s: %s", names(level_reasons), unlist(level_reasons)),
      collapse = " | "
    )
  } else {
    NA_character_
  }

  scitargets_dea(
    comparison_name = comparison_name,
    group_by = group_by,
    cluster_by = if (has_cluster) cluster_by else character(),
    cluster = cluster,
    groups = c(group1, group2),
    group1 = group1,
    group2 = group2,
    levels = computed_levels,
    pseudobulk_unit = if ("pseudobulk" %in% computed_levels) pseudobulk_unit else character(),
    n_cells = n_cells,
    de = de,
    go = go,
    gsea = gsea,
    pca = pca,
    de_params = de_params,
    go_params = params,
    gsea_params = gparams_store,
    species = species,
    padj_cutoffs = padj_cutoffs,
    status = status,
    reason = if (is.null(reason)) NA_character_ else reason
  )
}


#' Fetch MSigDB gene-set collections for GSEA
#'
#' @param collections Character vector of collection labels. Supported:
#'   "Hallmark", "GO:BP", "C7:ImmuneSigDB".
#' @param species Species name passed to [msigdbr::msigdbr()].
#' @returns A named list (collection-label -> named list of gene-symbol vectors).
#' @export
get_msigdbr_pathways <- function(collections = names(.MSIGDB_COLLECTION_MAP),
                                 species = "Homo sapiens") {
  out <- list()
  for (label in collections) {
    spec <- .MSIGDB_COLLECTION_MAP[[label]]
    if (is.null(spec)) {
      stop("Unknown GSEA collection label: ", label, call. = FALSE)
    }
    df <- if (is.null(spec$subcollection)) {
      msigdbr::msigdbr(species = species, collection = spec$collection)
    } else {
      msigdbr::msigdbr(
        species = species, collection = spec$collection,
        subcollection = spec$subcollection
      )
    }
    out[[label]] <- split(df$gene_symbol, df$gs_name)
  }
  out
}

# Run fgsea for a single collection.
compute_gsea_collection <- function(ranks, pathways, params, collection) {
  if (length(ranks) < 2L || length(pathways) == 0L) {
    return(.empty_gsea_result(
      "Not enough ranked genes or pathways for GSEA.",
      collection, ranks, params
    ))
  }

  fg <- tryCatch(
    suppressWarnings(fgsea::fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = params$minSize,
      maxSize = params$maxSize,
      eps = params$eps,
      # run serially: the targets pipeline already runs comparisons one at a
      # time, so avoid fgsea spawning a pool of ~16 BiocParallel worker R
      # processes (each holding RAM) for every collection/level.
      nproc = params$nproc %||% 1L
    )),
    error = function(e) NULL
  )

  if (is.null(fg) || nrow(fg) == 0L) {
    return(.empty_gsea_result(
      "fgsea returned no enriched gene sets.",
      collection, ranks, params
    ))
  }

  fg <- as.data.frame(fg)
  if ("leadingEdge" %in% colnames(fg)) {
    fg$leadingEdge <- vapply(
      fg$leadingEdge,
      function(g) paste(g, collapse = ", "),
      character(1)
    )
  }
  fg <- fg[order(fg$padj, fg$pval, na.last = TRUE), , drop = FALSE]
  rownames(fg) <- NULL

  list(
    status = "computed",
    reason = NA_character_,
    collection = collection,
    res.table = fg,
    ranks = ranks,
    params = params
  )
}

# Run GSEA for every requested collection on one level's marker table. The ranking
# metric is chosen by .gsea_ranks(): the signed DESeq2 Wald statistic when a `stat`
# column is present (pseudobulk Wald), otherwise sign(log2FC) * -log10(p).
compute_gsea <- function(markers, pathways_by_collection, gsea_params = list()) {
  params <- .default_gsea_params(gsea_params)
  ranks <- .gsea_ranks(markers)

  if (is.null(pathways_by_collection)) {
    pathways_by_collection <- get_msigdbr_pathways(params$collections, params$species)
  }

  res <- list()
  for (label in names(pathways_by_collection)) {
    res[[label]] <- compute_gsea_collection(
      ranks = ranks,
      pathways = pathways_by_collection[[label]],
      params = params,
      collection = label
    )
  }
  res
}

# Generics + methods -------------------------------------------------------

#' GSEA accessors for [scitargets_dea]
#'
#' @param x A [scitargets_dea] object.
#' @param ... Passed to methods (e.g. `level`, `collection`).
#' @returns `gsea_table` returns a data.frame; `gsea_barplot` returns a ggplot.
#' @name gsea_accessors
#' @export
gsea_table <- S7::new_generic("gsea_table", "x")
#' @rdname gsea_accessors
#' @export
gsea_barplot <- S7::new_generic("gsea_barplot", "x")

S7::method(gsea_table, scitargets_dea) <- function(x, collection = NULL, level = NULL, ...) {
  res <- .dea_gsea_result(x, level = level, collection = collection)
  if (is.null(res) || is.null(res$res.table)) data.frame() else res$res.table
}

S7::method(gsea_barplot, scitargets_dea) <- function(x,
                                                     collection = NULL,
                                                     level = NULL,
                                                     top_n = NULL,
                                                     padj_cutoff = NULL,
                                                     metric = c(
                                                       "signed_nlog10_padj",
                                                       "signed_nlog10_pval",
                                                       "signed_pval"
                                                     ),
                                                     interactive = FALSE,
                                                     width_svg = 10,
                                                     height_svg = 7,
                                                     ...) {
  metric <- match.arg(metric)
  res <- .dea_gsea_result(x, level = level, collection = collection)

  if (is.null(res) || is.null(res$res.table) || nrow(res$res.table) == 0L) {
    return(.no_result_plot(
      (if (is.null(res)) NULL else res$reason) %||% "No GSEA result available."
    ))
  }

  tab <- res$res.table
  params <- res$params
  pc <- padj_cutoff %||% params$padj_cutoff %||% 0.05
  tn <- top_n %||% params$top_nodes %||% 20L

  sig <- tab[!is.na(tab$padj) & tab$padj < pc, , drop = FALSE]
  subtitle <- paste0("gene sets with adjusted p-value < ", pc)
  if (nrow(sig) == 0L) {
    sig <- tab[!is.na(tab$pval), , drop = FALSE]
    subtitle <- paste0(
      "no gene set passed the adjusted p-value < ", pc,
      " cutoff; showing the top gene sets by raw p-value"
    )
  }
  if (nrow(sig) == 0L) {
    return(.no_result_plot(res$reason %||% "No GSEA gene set to display."))
  }

  pv <- pmax(sig$pval, .Machine$double.xmin)
  pa <- pmax(sig$padj, .Machine$double.xmin)
  sig$bar <- switch(metric,
    signed_nlog10_pval = sign(sig$NES) * -log10(pv),
    signed_nlog10_padj = sign(sig$NES) * -log10(pa),
    signed_pval = sign(sig$NES) * sig$pval
  )
  axis_label <- switch(metric,
    signed_nlog10_pval = "sign(NES) x -log10(p-value)",
    signed_nlog10_padj = "sign(NES) x -log10(adj. p-value)",
    signed_pval = "sign(NES) x p-value"
  )

  # Select the top_n gene sets by significance, independent of the display metric.
  # Ranking by abs(bar) would, for the "signed_pval" metric (bar = sign(NES) * pval),
  # pick the LARGEST p-values (the weakest sets) instead of the strongest.
  sig <- sig[order(sig$padj, sig$pval, na.last = TRUE), , drop = FALSE]
  if (nrow(sig) > tn) sig <- sig[seq_len(tn), , drop = FALSE]

  # The bar is colored by the group in which the gene set is enriched:
  # NES > 0 -> enriched in group1, NES < 0 -> enriched in group2 / rest.
  g1 <- x@group1
  g2 <- if (length(x@group2) && nzchar(x@group2)) x@group2 else "rest"
  sig$enriched_group <- ifelse(sig$NES > 0, g1, g2)

  sig <- sig[order(sig$bar), , drop = FALSE]
  sig$pathway_label <- stringr::str_wrap(gsub("_", " ", sig$pathway), width = 45)
  sig$pathway_label <- factor(sig$pathway_label, levels = unique(sig$pathway_label))

  # Build the tooltip column BEFORE constructing the plot, so the (interactive)
  # layer's data carries it.
  le <- if ("leadingEdge" %in% names(sig)) sig$leadingEdge else ""
  sig$tooltip <- paste0(
    "Pathway: ", sig$pathway,
    "\nNES: ", sprintf("%.2f", sig$NES),
    "\np-value: ", signif(sig$pval, 3),
    "\nadj. p-value: ", signif(sig$padj, 3),
    "\nLeading edge: ", le
  )

  base <- ggplot2::ggplot(
    sig,
    ggplot2::aes(x = pathway_label, y = bar, fill = enriched_group)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = axis_label,
      fill = "Enriched in",
      title = paste0("GSEA: ", res$collection),
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

  # ggiraph: per-bar tooltip with the pathway, p-value/padj, NES and leading-edge
  # genes. ggplot2 (default): static bars with the NES printed on each bar.
  if (isTRUE(interactive)) {
    p <- base + ggiraph::geom_col_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = pathway)
    )
    return(ggiraph::girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "fill:orange;stroke:black;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  base +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("NES = %.2f", NES)),
      size = 3,
      hjust = ifelse(sig$bar >= 0, 1.05, -0.05)
    )
}

#' Patient-level cell-state composition barplot
#'
#' Stacked barplot of each patient's cell-state composition (proportion of cells
#' per state), faceted by clinical group. Computed from a metadata data.frame
#' (whole dataset), so it is independent of any single comparison.
#'
#' @param meta A metadata data.frame (e.g. `seurat_obj@meta.data`).
#' @param state_col Metadata column holding the cell-state / cluster annotation
#'   (e.g. "clusters_0.3" or "clinical_azimuth_l2"). Cells with `NA` here are
#'   dropped.
#' @param group_col Clinical-group column (default "cohort_clinic_preterrah").
#' @param patient_col Patient / replicate column (default "patient_id").
#' @param interactive If TRUE, return a ggiraph girafe with per-segment tooltips
#'   (patient, group, state, cells, proportion); otherwise a static ggplot
#'   (default FALSE).
#' @param width_svg,height_svg girafe canvas size (interactive only).
#' @returns A ggplot (static) or girafe (interactive) object.
#' @export
composition_plot <- function(meta,
                             state_col,
                             group_col = "cohort_clinic_preterrah",
                             patient_col = "patient_id",
                             interactive = FALSE,
                             width_svg = NULL,
                             height_svg = 6) {
  comp <- .cell_state_composition(meta, state_col, group_col, patient_col)
  if (nrow(comp) == 0L) {
    return(.no_result_plot("No cells available for the composition plot."))
  }
  # Scale the (interactive) canvas width with the number of samples so the bars
  # stay readable when there are many patients.
  if (is.null(width_svg)) {
    width_svg <- max(10, ceiling(0.45 * length(unique(comp$patient))))
  }
  comp$tooltip <- paste0(
    "Patient: ", comp$patient, "\nGroup: ", comp$group,
    "\nState: ", comp$state,
    "\nCells: ", comp$n_cells, " / ", comp$n_total,
    "\nProportion: ", sprintf("%.1f%%", 100 * comp$proportion)
  )
  # order patients by clinical group, then id
  comp$patient <- factor(
    comp$patient,
    levels = unique(comp$patient[order(comp$group, comp$patient)])
  )

  base <- ggplot2::ggplot(comp, ggplot2::aes(x = patient, y = proportion, fill = state)) +
    ggplot2::facet_grid(~group, scales = "free_x", space = "free_x") +
    ggplot2::labs(
      x = NULL, y = "Proportion of cells", fill = "State",
      title = "Patient cell-state composition",
      subtitle = sprintf("state: %s | grouped by %s", state_col, group_col)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      # No legend: with many cell states it dominates the plot; the interactive
      # tooltips carry the state instead.
      legend.position = "none"
    )

  if (isTRUE(interactive)) {
    p <- base + ggiraph::geom_col_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = interaction(patient, state))
    )
    return(ggiraph::girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "stroke:black;stroke-width:1px;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  base + ggplot2::geom_col()
}

#' Cell-state proportion boxplot by clinical group
#'
#' Per-state boxplot of the patient-level cell-state proportions across clinical
#' groups, with jittered per-patient points. Computed from a metadata data.frame
#' (whole dataset).
#'
#' @param meta A metadata data.frame (e.g. `seurat_obj@meta.data`).
#' @param state_col Metadata column holding the cell-state / cluster annotation.
#'   Cells with `NA` here are dropped.
#' @param group_col Clinical-group column (default "cohort_clinic_preterrah").
#' @param patient_col Patient / replicate column (default "patient_id").
#' @param interactive If TRUE, return a ggiraph girafe whose jittered points carry
#'   per-patient tooltips; otherwise a static ggplot (default FALSE).
#' @param width_svg,height_svg girafe canvas size (interactive only).
#' @returns A ggplot (static) or girafe (interactive) object.
#' @export
composition_boxplot <- function(meta,
                                state_col,
                                group_col = "cohort_clinic_preterrah",
                                patient_col = "patient_id",
                                interactive = FALSE,
                                width_svg = 12,
                                height_svg = NULL) {
  comp <- .cell_state_composition(meta, state_col, group_col, patient_col)
  if (nrow(comp) == 0L) {
    return(.no_result_plot("No cells available for the composition boxplot."))
  }
  # Two facet columns; scale the (interactive) canvas height with the number of
  # cell-state facet rows so each panel stays readable when there are many states.
  if (is.null(height_svg)) {
    height_svg <- max(6, ceiling(length(unique(comp$state)) / 2) * 3L)
  }
  comp$tooltip <- paste0(
    "Patient: ", comp$patient, "\nGroup: ", comp$group,
    "\nState: ", comp$state,
    "\nProportion: ", sprintf("%.1f%%", 100 * comp$proportion)
  )

  base <- ggplot2::ggplot(comp, ggplot2::aes(x = group, y = proportion)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = group), outlier.shape = NA, alpha = 0.5) +
    ggplot2::facet_wrap(~state, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      x = NULL, y = "Proportion of cells", fill = "Group",
      title = "Cell-state proportion by clinical group",
      subtitle = sprintf("state: %s | grouped by %s", state_col, group_col)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  jit <- ggplot2::position_jitter(width = 0.15, height = 0, seed = 1L)
  if (isTRUE(interactive)) {
    p <- base + ggiraph::geom_point_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = patient),
      position = jit, alpha = 0.8, size = 1.8
    )
    return(ggiraph::girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "fill:orange;stroke:black;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  base + ggplot2::geom_jitter(position = jit, alpha = 0.7, size = 1.8)
}

#' UpSet plot of DE-gene overlap across comparisons
#'
#' Builds an UpSet plot (via \pkg{ggupset}) of the differentially-expressed genes
#' shared across comparisons, for one regulation direction.
#'
#' @param x A list of [scitargets_dea] objects (e.g. [normalize_dea_list()]
#'   output), or anything [normalize_dea_list()] accepts.
#' @param direction "both" (any DE gene), "up", or "down".
#' @param level Differential-expression level ("single_cell"/"pseudobulk");
#'   `NULL` uses each object's default (first) level.
#' @param comparisons `NA` (default) = all comparisons that have DE genes, or a
#'   character vector of `comparison_name`s to restrict to.
#' @param padj_cutoff Adjusted-p cutoff for calling a gene DE; `NULL` uses each
#'   object's stored per-level cutoff.
#' @returns A ggplot (UpSet) object, or a placeholder when there are no DE genes.
#' @export
dea_upset_plot <- function(x,
                           direction = c("both", "up", "down"),
                           level = NULL,
                           comparisons = NA,
                           padj_cutoff = NULL) {
  direction <- match.arg(direction)
  if (!requireNamespace("ggupset", quietly = TRUE)) {
    stop("Package 'ggupset' is required for UpSet plots; install it with ",
      "install.packages('ggupset').",
      call. = FALSE
    )
  }
  df <- .dea_upset_data(x,
    direction = direction, level = level,
    comparisons = comparisons, padj_cutoff = padj_cutoff
  )
  if (nrow(df) == 0L) {
    return(.no_result_plot("No differentially expressed genes in the selected comparisons."))
  }
  ttl <- switch(direction,
    up = "up-regulated",
    down = "down-regulated",
    both = "all"
  )
  ggplot2::ggplot(df, ggplot2::aes(x = comparisons)) +
    ggplot2::geom_bar(fill = "steelblue") +
    ggplot2::geom_text(
      stat = "count",
      ggplot2::aes(label = ggplot2::after_stat(count)),
      vjust = -0.3, size = 3
    ) +
    ggupset::scale_x_upset() +
    ggplot2::labs(
      x = NULL, y = "Number of DE genes",
      title = sprintf("DE-gene overlap across comparisons (%s)", ttl)
    ) +
    ggplot2::theme_bw(base_size = 12)
}

#' Build Quarto child lines for the DE-gene UpSet plots
#'
#' Generates Quarto/knitr lines that render one UpSet plot per requested
#' direction (down / up / both). Mirrors [dea_report_lines()]: pass the result to
#' [knitr::knit_child()] with `obj_name` bound (to the comparison list) in the
#' environment. Comparisons with no DE genes for a direction emit a short note
#' instead of an empty plot.
#'
#' @param x A list of [scitargets_dea] objects (used at generation time to detect
#'   whether any DE genes exist for each direction).
#' @param obj_name Name of the variable the generated chunks reference (must be
#'   bound to the same comparison list in the knit environment; default
#'   "clinical_comparisons").
#' @param directions Which UpSet plots to emit: any of "down", "up", "both"
#'   (default all three).
#' @param level Differential-expression level passed to [dea_upset_plot()].
#' @param comparisons `NA` (default, all comparisons with DE genes) or a character
#'   vector of `comparison_name`s to restrict to.
#' @param padj_cutoff Adjusted-p cutoff passed to [dea_upset_plot()].
#' @param fig_dims Optional named list keyed by direction, each a named vector or
#'   list with `width` and/or `height` (Quarto `fig-width`/`fig-height`), e.g.
#'   `list(both = c(width = 10, height = 7))`. Missing values fall back to
#'   8 x 6.
#' @param heading_level Markdown heading level for the section title (default 2).
#' @param heading Section heading text.
#' @param fig_id A unique id used to build the figure chunk labels.
#' @returns A character vector of Quarto/knitr lines.
#' @export
dea_upset_lines <- function(x,
                            obj_name = "clinical_comparisons",
                            directions = c("down", "up", "both"),
                            level = NULL,
                            comparisons = NA,
                            padj_cutoff = NULL,
                            fig_dims = NULL,
                            heading_level = 2L,
                            heading = "DE-gene overlap (UpSet)",
                            fig_id = "upset") {
  directions <- match.arg(directions, several.ok = TRUE)
  h <- paste(rep("#", heading_level), collapse = "")
  h2 <- paste(rep("#", heading_level + 1L), collapse = "")
  dlab <- c(down = "Down-regulated", up = "Up-regulated", both = "All DE genes")
  lit <- function(v) paste(deparse(v), collapse = " ")
  getd <- function(fd, key, default) {
    if (is.null(fd)) {
      return(default)
    }
    v <- tryCatch(fd[[key]], error = function(e) NULL)
    if (is.null(v) || is.na(v)) default else v
  }

  lines <- c(paste(h, heading), "")
  for (dir in directions) {
    lines <- c(lines, paste(h2, dlab[[dir]]), "")
    has_de <- nrow(.dea_upset_data(x,
      direction = dir, level = level,
      comparisons = comparisons,
      padj_cutoff = padj_cutoff
    )) > 0L
    if (!has_de) {
      lines <- c(lines, "No DE genes were found in the comparisons to plot.", "")
      next
    }
    fd <- if (is.null(fig_dims)) NULL else fig_dims[[dir]]
    lines <- c(
      lines,
      "```{r}",
      "#| message: false",
      "#| warning: false",
      sprintf("#| fig-width: %s", getd(fd, "width", 8)),
      sprintf("#| fig-height: %s", getd(fd, "height", 6)),
      sprintf("#| label: fig-%s-%s", .dea_sanitize_id(fig_id, dir), dir),
      sprintf(
        "dea_upset_plot(%s, direction = %s, level = %s, comparisons = %s, padj_cutoff = %s)",
        obj_name, shQuote(dir), lit(level), lit(comparisons), lit(padj_cutoff)
      ),
      "```",
      ""
    )
  }
  lines
}

#' Testability heatmap across comparisons
#'
#' Heatmap of a per-(cluster x comparison) testability metric across a list of
#' comparisons (the supervisor's testability summary).
#'
#' @param x A list of [scitargets_dea] objects (or anything
#'   [normalize_dea_list()] accepts).
#' @param metric "testable" (>= min replicates per group), "patients" (min
#'   replicates per group), "cells" (min cells per group), or "outliers" (count of
#'   flagged outlier replicates removed). All computed from the stored PCA coords.
#' @param level Differential-expression level (default "pseudobulk").
#' @param after_removal For "testable"/"patients"/"cells": `TRUE` (default)
#'   excludes the flagged outlier replicates (the post-removal view DESeq2 used);
#'   `FALSE` keeps them (before removal). Ignored for `metric = "outliers"`.
#' @param interactive If TRUE, return a ggiraph girafe with per-tile tooltips.
#' @param width_svg,height_svg girafe canvas size (interactive only).
#' @returns A ggplot (static) or girafe (interactive) object.
#' @export
dea_testability_heatmap <- function(x,
                                    metric = c("testable", "patients", "cells", "outliers"),
                                    level = "pseudobulk",
                                    after_removal = TRUE,
                                    interactive = FALSE,
                                    width_svg = 10,
                                    height_svg = 6) {
  metric <- match.arg(metric)
  df <- .dea_testability_data(x, level = level, metric = metric, after_removal = after_removal)
  if (is.null(df) || nrow(df) == 0L) {
    return(.no_result_plot("No comparisons available for the testability heatmap."))
  }
  lab <- switch(metric,
    testable = "Testable",
    patients = "Min replicates/group",
    cells = "Min cells/group",
    outliers = "Outliers"
  )
  ttl <- switch(metric,
    testable = "Testable: cell state x comparison",
    patients = "Patients per group: cell state x comparison",
    cells = "Cells per group: cell state x comparison",
    outliers = "Outliers removed: cell state x comparison"
  )
  # patients/cells/testable have a before/after-removal meaning; outliers does not.
  stage <- if (metric == "outliers") {
    "flagged outlier replicates removed from DESeq2"
  } else if (isTRUE(after_removal)) {
    "after outlier removal"
  } else {
    "before outlier removal"
  }
  df$label <- if (metric == "testable") {
    ifelse(df$value > 0, "yes", "no")
  } else {
    ifelse(is.na(df$value), "NA", format(df$value))
  }
  df$tooltip <- paste0(
    "Cell state: ", df$cluster, "\nComparison: ", df$comparison,
    "\n", lab, ": ", df$label
  )

  base <- ggplot2::ggplot(df, ggplot2::aes(x = comparison, y = cluster, fill = value)) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "#2c7fb8", na.value = "white") +
    ggplot2::labs(
      x = NULL, y = NULL, fill = lab, title = ttl,
      subtitle = sprintf("level: %s | %s", level, stage)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  if (isTRUE(interactive)) {
    p <- base +
      ggiraph::geom_tile_interactive(
        ggplot2::aes(tooltip = tooltip, data_id = interaction(cluster, comparison)),
        color = "white"
      ) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 3)
    return(ggiraph::girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "stroke:black;stroke-width:1px;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  base +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3)
}

#' Build Quarto child lines for the testability heatmaps
#'
#' Quarto/knitr lines rendering one testability heatmap per requested metric.
#' Mirrors [dea_upset_lines()]: pass to [knitr::knit_child()] with `obj_name`
#' bound to the comparison list.
#'
#' @param x A list of [scitargets_dea] objects.
#' @param obj_name Variable the generated chunks reference (default
#'   "clinical_comparisons").
#' @param metrics Any of "testable", "patients", "cells" (default all three);
#'   each is shown as a Before/After-outlier-removal tabset when any comparison
#'   has flagged outliers. An "Outliers removed" heatmap is appended in that case.
#' @param level Differential-expression level (default "pseudobulk").
#' @param interactive Render the heatmaps interactively (ggiraph).
#' @param fig_dims Optional named list keyed by metric with `width`/`height`.
#' @param heading_level Section heading level (default 2).
#' @param heading Section heading text.
#' @param fig_id Unique id for the chunk labels.
#' @returns A character vector of Quarto/knitr lines.
#' @export
dea_testability_lines <- function(x,
                                  obj_name = "clinical_comparisons",
                                  metrics = c("testable", "patients", "cells"),
                                  level = "pseudobulk",
                                  interactive = FALSE,
                                  fig_dims = NULL,
                                  heading_level = 2L,
                                  heading = "Testability summary",
                                  fig_id = "testability") {
  metrics <- match.arg(metrics, several.ok = TRUE)
  h <- .dea_h(heading_level)
  h2 <- .dea_h(heading_level + 1L)
  h3 <- .dea_h(heading_level + 2L)
  mlab <- c(
    testable = "Testable state x comparison",
    patients = "Patients per group",
    cells = "Cells per group"
  )
  ibool <- if (isTRUE(interactive)) "TRUE" else "FALSE"
  getd <- function(fd, key, default) {
    if (is.null(fd)) {
      return(default)
    }
    v <- tryCatch(fd[[key]], error = function(e) NULL)
    if (is.null(v) || is.na(v)) default else v
  }
  chunk <- function(m, after_removal, fd, tag) {
    c(
      "```{r}",
      "#| message: false",
      "#| warning: false",
      sprintf("#| fig-width: %s", getd(fd, "width", 8)),
      sprintf("#| fig-height: %s", getd(fd, "height", 6)),
      sprintf("#| label: fig-%s", .dea_sanitize_id(fig_id, m, tag)),
      sprintf(
        "dea_testability_heatmap(%s, metric = %s, level = %s, after_removal = %s, interactive = %s)",
        obj_name, shQuote(m), shQuote(level),
        if (isTRUE(after_removal)) "TRUE" else "FALSE", ibool
      ),
      "```",
      ""
    )
  }

  has_out <- .dea_any_outliers(x, level = level)
  lines <- c(paste(h, heading), "")
  if (has_out) {
    lines <- c(lines, sprintf(
      "_Outliers were detected; %s metrics are shown before and after their removal (the post-removal view matches what DESeq2 used)._",
      paste(intersect(metrics, c("patients", "cells", "testable")), collapse = " / ")
    ), "")
  }
  for (m in metrics) {
    fd <- if (is.null(fig_dims)) NULL else fig_dims[[m]]
    lines <- c(lines, paste(h2, mlab[[m]]), "")
    if (has_out) {
      lines <- c(
        lines, "::: {.panel-tabset}", "",
        paste(h3, "Before outlier removal"), "", chunk(m, FALSE, fd, "before"),
        paste(h3, "After outlier removal"), "", chunk(m, TRUE, fd, "after"),
        ":::", ""
      )
    } else {
      lines <- c(lines, chunk(m, TRUE, fd, "after"))
    }
  }
  # how many outlier replicates were removed per comparison (only if any)
  if (has_out) {
    fd <- if (is.null(fig_dims)) NULL else fig_dims[["outliers"]]
    lines <- c(lines, paste(h2, "Outliers removed"), "", chunk("outliers", TRUE, fd, "count"))
  }
  lines
}

#' GSEA-yield heatmap across comparisons
#'
#' Heatmap of the number of significantly enriched gene sets
#' (`padj < padj_cutoff`) per (cell state x comparison) for one MSigDB
#' collection -- the supervisor's "GSEA-yield" summary, showing where (and how
#' much) pathway-level signal appeared across the comparison grid. A tile is
#' white (`NA`) when the collection was not computed for that comparison, and
#' grey / 0 when it was computed but no gene set passed the cutoff.
#'
#' @param x A list of [scitargets_dea] objects (or anything
#'   [normalize_dea_list()] accepts).
#' @param collection MSigDB collection name (e.g. "Hallmark", "GO:BP", "C7").
#' @param level Differential-expression level (default "pseudobulk").
#' @param padj_cutoff Adjusted-p cutoff for counting a gene set as enriched;
#'   `NULL` (default) uses each comparison's stored GSEA cutoff (else 0.05).
#' @param interactive If TRUE, return a ggiraph girafe with per-tile tooltips.
#' @param width_svg,height_svg girafe canvas size (interactive only).
#' @returns A ggplot (static) or girafe (interactive) object.
#' @export
dea_gsea_yield_heatmap <- function(x,
                                   collection,
                                   level = "pseudobulk",
                                   padj_cutoff = NULL,
                                   interactive = FALSE,
                                   width_svg = 10,
                                   height_svg = 6) {
  df <- .dea_gsea_yield_data(x,
    collection = collection, level = level,
    padj_cutoff = padj_cutoff
  )
  if (is.null(df) || nrow(df) == 0L || !any(!is.na(df$value))) {
    return(.no_result_plot(sprintf(
      "No GSEA results for collection '%s' (%s level).", collection, level
    )))
  }
  cutoff_disp <- if (is.null(padj_cutoff)) "comparison cutoff" else format(padj_cutoff)
  df$label <- ifelse(is.na(df$value), "NA", format(df$value))
  df$tooltip <- paste0(
    "Cell state: ", df$cluster, "\nComparison: ", df$comparison,
    "\nEnriched gene sets: ", df$label
  )

  base <- ggplot2::ggplot(df, ggplot2::aes(x = comparison, y = cluster, fill = value)) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "#238b45", na.value = "white") +
    ggplot2::labs(
      x = NULL, y = NULL, fill = "Gene sets",
      title = sprintf("Significantly enriched gene sets: %s", collection),
      subtitle = sprintf("level: %s | padj < %s", level, cutoff_disp)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  if (isTRUE(interactive)) {
    p <- base +
      ggiraph::geom_tile_interactive(
        ggplot2::aes(tooltip = tooltip, data_id = interaction(cluster, comparison)),
        color = "white"
      ) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 3)
    return(ggiraph::girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_hover(css = "stroke:black;stroke-width:1px;"),
        ggiraph::opts_tooltip(
          css = "background-color:white;color:black;padding:6px;border:1px solid grey;"
        )
      )
    ))
  }

  base +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3)
}

#' Build Quarto child lines for the GSEA-yield heatmaps
#'
#' Quarto/knitr lines rendering one GSEA-yield heatmap per MSigDB collection (the
#' number of significantly enriched gene sets per cell state x comparison).
#' Mirrors [dea_upset_lines()] / [dea_testability_lines()]: pass the result to
#' [knitr::knit_child()] with `obj_name` bound (to the comparison list) in the
#' environment. When no GSEA was computed at `level`, a short note is emitted.
#'
#' @param x A list of [scitargets_dea] objects (used at generation time to find
#'   the collections present at `level`).
#' @param obj_name Name of the variable the generated chunks reference (must be
#'   bound to the same comparison list in the knit environment; default
#'   "clinical_comparisons").
#' @param collections `NULL` (default) = every collection present across the
#'   comparisons at `level`, or a character vector of collection names to
#'   restrict to (in the given order).
#' @param level Differential-expression level (default "pseudobulk").
#' @param padj_cutoff Adjusted-p cutoff passed to [dea_gsea_yield_heatmap()];
#'   `NULL` uses each comparison's stored GSEA cutoff.
#' @param interactive Render the heatmaps interactively (ggiraph).
#' @param fig_dims Optional named list keyed by collection, each a named vector
#'   or list with `width` and/or `height` (Quarto `fig-width`/`fig-height`).
#'   Missing values fall back to 8 x 6.
#' @param heading_level Markdown heading level for the section title (default 2).
#' @param heading Section heading text.
#' @param fig_id A unique id used to build the figure chunk labels.
#' @returns A character vector of Quarto/knitr lines.
#' @export
dea_gsea_yield_lines <- function(x,
                                 obj_name = "clinical_comparisons",
                                 collections = NULL,
                                 level = "pseudobulk",
                                 padj_cutoff = NULL,
                                 interactive = FALSE,
                                 fig_dims = NULL,
                                 heading_level = 2L,
                                 heading = "GSEA yield (enriched gene sets)",
                                 fig_id = "gsea-yield") {
  h <- paste(rep("#", heading_level), collapse = "")
  h2 <- paste(rep("#", heading_level + 1L), collapse = "")
  ibool <- if (isTRUE(interactive)) "TRUE" else "FALSE"
  lit <- function(v) paste(deparse(v), collapse = " ")
  getd <- function(fd, key, default) {
    if (is.null(fd)) {
      return(default)
    }
    v <- tryCatch(fd[[key]], error = function(e) NULL)
    if (is.null(v) || is.na(v)) default else v
  }

  cols <- collections %||% .dea_gsea_yield_collections(x, level = level)
  lines <- c(paste(h, heading), "")
  if (length(cols) == 0L) {
    return(c(lines, "No GSEA results were found in the comparisons to plot.", ""))
  }
  for (col in cols) {
    lines <- c(lines, paste(h2, col), "")
    has_yield <- any(!is.na(.dea_gsea_yield_data(
      x,
      collection = col, level = level, padj_cutoff = padj_cutoff
    )$value))
    if (!has_yield) {
      lines <- c(lines, "No GSEA results were found in the comparisons to plot.", "")
      next
    }
    fd <- if (is.null(fig_dims)) NULL else fig_dims[[col]]
    lines <- c(
      lines,
      "```{r}",
      "#| message: false",
      "#| warning: false",
      "#| column: page",
      sprintf("#| fig-width: %s", getd(fd, "width", 8)),
      sprintf("#| fig-height: %s", getd(fd, "height", 6)),
      sprintf("#| label: fig-%s", .dea_sanitize_id(fig_id, level, col)),
      sprintf(
        "dea_gsea_yield_heatmap(%s, collection = %s, level = %s, padj_cutoff = %s, interactive = %s)",
        obj_name, shQuote(col), shQuote(level), lit(padj_cutoff), ibool
      ),
      "```",
      ""
    )
  }
  lines
}

#' Build Quarto child lines for the cell-state composition plots
#'
#' Quarto/knitr lines for the patient cell-state composition barplot +
#' proportion boxplot, shown ONCE (e.g. before the per-comparison sections). The
#' generated chunks call [composition_plot()] / [composition_boxplot()] on a
#' metadata object bound in the knit environment under `meta_name`.
#'
#' @param meta_name Name of the metadata variable in the knit environment
#'   (default "clinical_meta").
#' @param state_col,group_col,patient_col Columns passed to the composition
#'   functions (state default has none; group/patient default to
#'   "cohort_clinic_preterrah" / "patient_id").
#' @param interactive Render the plots interactively (ggiraph) or static.
#' @param heading_level,heading Section heading level / text (for the TOC).
#' @param fig_id Unique id for the chunk labels.
#' @returns A character vector of Quarto/knitr lines.
#' @export
composition_lines <- function(meta_name = "clinical_meta",
                              state_col,
                              group_col = "cohort_clinic_preterrah",
                              patient_col = "patient_id",
                              interactive = FALSE,
                              heading_level = 2L,
                              heading = "Cell-state composition",
                              fig_id = "composition") {
  h <- paste(rep("#", heading_level), collapse = "")
  q <- function(s) shQuote(s)
  cargs <- sprintf(
    "%s, state_col = %s, group_col = %s, patient_col = %s, interactive = %s",
    meta_name, q(state_col), q(group_col), q(patient_col),
    if (isTRUE(interactive)) "TRUE" else "FALSE"
  )
  c(
    paste(h, heading), "",
    "```{r}", "#| message: false", "#| warning: false", "#| column: screen-outset",
    sprintf("#| label: fig-%s-bar", .dea_sanitize_id(fig_id)),
    "#| fig-cap: \"Patient cell-state composition (stacked proportions).\"",
    sprintf("composition_plot(%s)", cargs),
    "```", "",
    "```{r}", "#| message: false", "#| warning: false", "#| column: screen-outset",
    sprintf("#| label: fig-%s-box", .dea_sanitize_id(fig_id)),
    "#| fig-cap: \"Cell-state proportion by clinical group.\"",
    sprintf("composition_boxplot(%s)", cargs),
    "```", ""
  )
}

#' Heatmap of the top DE genes across all cells (Seurat)
#'
#' Plots a [Seurat::DoHeatmap()] of the given genes across all cells of a Seurat
#' object (cell-level expression). Genes not present in the object are dropped;
#' the requested genes are scaled on the fly.
#'
#' @param seurat_obj A Seurat object (e.g. `seurat_merged_filtered`). If it is not
#'   a valid Seurat object a placeholder is returned.
#' @param features Character vector of genes to show.
#' @param group_by Metadata column to group cells by (default the active idents).
#' @returns A ggplot (DoHeatmap), or a placeholder when nothing can be plotted.
#' @export
top_de_heatmap <- function(seurat_obj, features, group_by = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    return(.no_result_plot("No valid Seurat object for the top-DE heatmap."))
  }
  features <- unique(features[features %in% rownames(seurat_obj)])
  if (length(features) == 0L) {
    return(.no_result_plot("No top DE genes available for the heatmap."))
  }
  # Seurat::DoHeatmap() mismatches its internal label vectors when the grouping
  # column has NA cells: it builds one x position per plotted cell but uses
  # sort(group.use), and sort() silently drops NAs -> "differing number of rows"
  # (off by the NA count). Keep only cells that have a non-missing group label.
  if (!is.null(group_by) && group_by %in% colnames(seurat_obj[[]])) {
    grp <- seurat_obj[[group_by]][, 1]
    keep <- colnames(seurat_obj)[!is.na(grp)]
    if (length(keep) == 0L) {
      return(.no_result_plot("No cells with a non-missing group label for the heatmap."))
    }
    if (length(keep) < ncol(seurat_obj)) {
      seurat_obj <- subset(seurat_obj, cells = keep)
    }
  }
  obj <- Seurat::ScaleData(seurat_obj, features = features, verbose = FALSE)
  Seurat::DoHeatmap(obj, features = features, group.by = group_by %||% "ident") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
}

#' Build Quarto child lines for the top-DE-gene heatmap
#'
#' Quarto/knitr lines for a single [top_de_heatmap()] shown AFTER all comparison
#' tabs. The genes are the union of the top-`n_top` DE genes (by |fold-change|
#' among significant genes) across the comparisons that the report shows (those
#' passing `clusters_to_show` / `groups_to_show`). Emitted only when `seurat_obj`
#' is a valid Seurat object AND at least one shown comparison has DE genes.
#'
#' @param x A list of [scitargets_dea] objects.
#' @param seurat_obj The Seurat object used for the cell-level heatmap, checked at
#'   generation time; if not a Seurat object nothing is emitted.
#' @param seurat_name Name of the Seurat variable the generated chunk references
#'   (must be bound in the knit environment; default "seurat_merged_filtered").
#' @param level Differential-expression level ("single_cell"/"pseudobulk");
#'   `NULL` = each object's default.
#' @param padj_cutoff Adjusted-p cutoff for "DE"; `NULL` = each object's per-level
#'   cutoff.
#' @param clusters_to_show,groups_to_show Same display filters as
#'   [dea_report_lines()]; the heatmap pools genes only from the shown comparisons.
#' @param group_by Metadata column to group the heatmap cells by.
#' @param chunk_opts Named list of Quarto chunk options for the heatmap chunk:
#'   `n_top` (number of top genes per comparison, default 20), plus `column`
#'   (default "screen-outset"), `label`, `caption`, `fig_width`, `fig_height`.
#' @param heading_level,heading Section heading level / text (for the TOC).
#' @param fig_id Unique id for the chunk label.
#' @returns A character vector of Quarto/knitr lines (empty when nothing to show).
#' @export
dea_top_de_lines <- function(x,
                             seurat_obj = NULL,
                             seurat_name = "seurat_merged_filtered",
                             level = NULL,
                             padj_cutoff = NULL,
                             clusters_to_show = NA,
                             groups_to_show = NA,
                             group_by = NULL,
                             chunk_opts = list(),
                             heading_level = 2L,
                             heading = "Top DE-gene heatmap",
                             fig_id = "topde") {
  if (!inherits(seurat_obj, "Seurat")) {
    return(character())
  }
  n_top <- chunk_opts[["n_top"]] %||% 20L
  genes <- .dea_top_de_genes(x,
    level = level, n_top = n_top, padj_cutoff = padj_cutoff,
    clusters_to_show = clusters_to_show, groups_to_show = groups_to_show
  )
  if (length(genes) == 0L) {
    return(character())
  }

  h <- paste(rep("#", heading_level), collapse = "")
  lit <- function(v) paste(deparse(v), collapse = " ")
  opts <- .dea_chunk_opts(
    chunk_opts,
    default_label = sprintf("fig-%s", .dea_sanitize_id(fig_id, level %||% "")),
    default_caption = sprintf(
      "Top DE genes (top %d per shown comparison by |log2FC|, %s level) across all cells.",
      n_top, level %||% "default"
    )
  )
  c(
    paste(h, heading), "",
    "```{r}",
    opts,
    sprintf(
      "top_de_heatmap(%s, features = %s, group_by = %s)",
      seurat_name, lit(genes), lit(group_by)
    ),
    "```",
    ""
  )
}

#' Build Quarto child lines for one [scitargets_dea] comparison
#'
#' Generates the markdown / knitr child text that renders a full comparison
#' block: a heading, then one section per computed level (single-cell, then
#' pseudobulk), each containing a tabset of volcano plot, DE table, GO enrichment
#' (per ontology and direction) and GSEA (per collection) shown side-by-side.
#' Pass the result to [knitr::knit_child()] with an environment in which
#' `obj_name` is bound to the object.
#'
#' @param x A [scitargets_dea] object.
#' @param fig_id A unique id used to build figure/table chunk labels.
#' @param cluster_word How to refer to the cluster in the heading (e.g.
#'   "cluster" or "Azimuth L2 annotation").
#' @param obj_name Name of the variable the generated chunks reference (must be
#'   bound in the knit environment; defaults to "comparison").
#' @param heading_level Markdown heading level for the comparison title.
#' @param levels Which differential-expression level(s) to render: `NULL`
#'   (default) for all computed levels, or a subset of
#'   `c("single_cell", "pseudobulk")`. Use this to split the report into a
#'   single-cell document and a pseudobulk document.
#' @param meta Optional metadata data.frame (e.g. `seurat_obj@meta.data`) used to
#'   build the per-replicate cell counts in the "Comparison summary" tab.
#' @param interactive_plots Which plot types to render interactively (ggiraph).
#'   `NA` (default) renders every plot as a static ggplot2; otherwise a character
#'   vector among `c("composition", "pca", "volcano", "go", "gsea")` naming the
#'   plots to make interactive (others stay static).
#' @param gsea_metric Bar metric for the GSEA barplot. `"signed_nlog10_padj"`
#'   (default) plots `sign(NES) * -log10(adj. p-value)`, i.e. the FDR-adjusted
#'   enrichment significance of each gene set; `"signed_nlog10_pval"` uses the
#'   raw enrichment p-value and `"signed_pval"` the signed raw p-value. Note this
#'   only controls the *display*: the gene-level ranking fed to GSEA always uses
#'   the raw per-gene p-value, independent of this argument.
#' @param composition_state_col Optional metadata column for the patient
#'   cell-state composition barplot shown in the comparison header, before the
#'   analysis tabset (e.g. "clusters_0.3" or "clinical_azimuth_l2"). `NULL`
#'   (default) omits it. The generated chunk calls [composition_plot()] on a
#'   `meta` object that must be bound in the knit environment.
#' @param composition_group_col,composition_patient_col Clinical-group and
#'   patient columns for the composition plot (defaults "cohort_clinic_preterrah"
#'   / "patient_id").
#' @param clusters_to_show Display filter (default `NA` = all clusters). A
#'   character vector of cluster names; the comparison is rendered only if its
#'   cluster is listed (returns `character()` to skip otherwise).
#' @param groups_to_show Display filter (default `NA` = all groups). A character
#'   vector of group names: 1 name -> any comparison involving that group; exactly
#'   2 -> only the comparison between those two; >2 -> any comparison whose BOTH
#'   groups are listed.
#' @param available_clusters,available_groups Optional full sets of valid cluster
#'   / group names. When supplied, a `clusters_to_show` / `groups_to_show` entry
#'   not in these sets produces a warning callout instead of silently filtering.
#' @returns A character vector of Quarto/knitr lines.
#' @export
dea_report_lines <- function(x,
                             fig_id,
                             cluster_word = "cluster",
                             obj_name = "comparison",
                             heading_level = 2L,
                             levels = NULL,
                             meta = NULL,
                             interactive_plots = NA,
                             composition_state_col = NULL,
                             composition_group_col = "cohort_clinic_preterrah",
                             composition_patient_col = "patient_id",
                             clusters_to_show = NA,
                             groups_to_show = NA,
                             available_clusters = NULL,
                             available_groups = NULL,
                             gsea_metric = c(
                               "signed_nlog10_padj",
                               "signed_nlog10_pval",
                               "signed_pval"
                             )) {
  stopifnot(is_scitargets_dea(x))
  gsea_metric <- match.arg(gsea_metric)
  # `interactive_plots`: NA (default) -> every plot is a static ggplot2; otherwise
  # a character vector among c("pca","volcano","go","gsea") naming the plot types
  # to render interactively (ggiraph). Validate up front.
  interactive_plots <- .dea_interactive_set(interactive_plots)
  h <- paste(rep("#", heading_level), collapse = "")
  title <- sprintf(
    "%s %s versus %s in %s %s",
    h, x@group1, x@group2, cluster_word, x@cluster
  )

  # P2.14.2 display filters. If `available_*` are supplied, invalid requested
  # names produce a WARNING callout (listing the valid names). Otherwise the
  # comparison is silently skipped (return character()) when it does not match.
  bad <- .dea_filter_invalid(
    clusters_to_show, groups_to_show,
    available_clusters, available_groups
  )
  if (length(bad) > 0L) {
    return(c(
      title, "", "::: {.callout-warning}",
      "#### Invalid filter argument(s)", "", bad, "", ":::", ""
    ))
  }
  if (!.dea_show_comparison(x, clusters_to_show, groups_to_show)) {
    return(character())
  }

  # Patient cell-state composition barplot, emitted in the comparison header
  # BEFORE the analysis tabset(s). Whole-dataset (from `meta`, which must be bound
  # in the knit environment); independent of the comparison. NULL state col omits it.
  comp_lines <- character()
  if (!is.null(composition_state_col)) {
    comp_int <- if ("composition" %in% interactive_plots) "TRUE" else "FALSE"
    cargs <- sprintf(
      "meta, state_col = %s, group_col = %s, patient_col = %s, interactive = %s",
      shQuote(composition_state_col), shQuote(composition_group_col),
      shQuote(composition_patient_col), comp_int
    )
    comp_lines <- c(
      "```{r}",
      "#| echo: false",
      "#| message: false",
      "#| warning: false",
      "#| column: page",
      sprintf("#| label: fig-composition-bar-%s", .dea_sanitize_id(fig_id)),
      sprintf("composition_plot(%s)", cargs),
      "```",
      "",
      "```{r}",
      "#| echo: false",
      "#| message: false",
      "#| warning: false",
      "#| column: page",
      sprintf("#| label: fig-composition-box-%s", .dea_sanitize_id(fig_id)),
      sprintf("composition_boxplot(%s)", cargs),
      "```",
      ""
    )
  }

  # Not computed at all: show the comparison summary (groups/cluster/columns +
  # cell counts) so the reader can see why it was skipped.
  if (!identical(x@status, "computed")) {
    return(c(title, "", comp_lines, .dea_summary_content(x, meta = meta)))
  }

  show_levels <- if (is.null(levels)) .dea_levels(x) else intersect(levels, .dea_levels(x))

  # The requested level(s) were not computed for this comparison (e.g. pseudobulk
  # with too few biological replicates): still show the comparison summary + a note.
  if (length(show_levels) == 0L) {
    want <- paste(vapply(levels, .dea_level_label, character(1)), collapse = " / ")
    return(c(
      title, "",
      comp_lines,
      "::: {.callout-note}",
      sprintf("#### %s not available", want), "",
      sprintf("No %s result for this comparison (e.g. too few replicates).", want), "",
      ":::", "",
      .dea_summary_content(x, meta = meta)
    ))
  }

  header_lines <- c(
    title, "",
    sprintf("**Group 1:** `%s` (positive log2FC / NES)", x@group1),
    "",
    sprintf("**Group 2:** `%s` (reference)", x@group2),
    "",
    sprintf("**%s:** `%s`", cluster_word, x@cluster),
    "",
    sprintf("**Levels shown:** %s.", paste(vapply(show_levels, .dea_level_label, character(1)),
      collapse = ", "
    )),
    ""
  )

  # With a single level (the usual case for the split single-cell / pseudobulk
  # documents) the level heading is redundant, so the analysis tabset starts one
  # level below the comparison title. With several levels, each gets its own
  # section heading and the tabset starts one level deeper.
  single_level <- length(show_levels) == 1L
  base_level <- if (single_level) heading_level + 1L else heading_level + 2L

  body <- character()
  for (level in show_levels) {
    if (!single_level) {
      body <- c(body, sprintf("%s %s", .dea_h(heading_level + 1L), .dea_level_label(level)), "")
    }
    body <- c(
      body,
      .dea_level_lines(x,
        level = level, fig_id = fig_id,
        obj_name = obj_name, base_level = base_level, meta = meta,
        interactive_plots = interactive_plots,
        gsea_metric = gsea_metric
      )
    )
  }

  c(header_lines, comp_lines, body)
}

#' Build the "DEA Comparisons" report section (all comparisons as child tabs)
#'
#' Emits one report section -- a parent tab holding every comparison's
#' [dea_report_lines()] block as a nested child tab -- analogous to the single
#' tab that [composition_lines()] or [dea_testability_lines()] each produce. Call
#' it once instead of looping [dea_report_lines()] by hand, so the per-comparison
#' tabs are grouped under a single "DEA Comparisons" tab (a sibling of the
#' "Cell-state composition", "Testability summary", ... tabs).
#'
#' The generated per-comparison chunks reference each object as
#' `obj_name[["<comparison name>"]]`, so only the LIST variable (`obj_name`) needs
#' to be bound in the environment passed to [knitr::knit_child()].
#'
#' @param x A named list of [scitargets_dea] objects (or anything
#'   [normalize_dea_list()] accepts that yields names). The names index the
#'   chunks and label the child tabs.
#' @param obj_name Name of the LIST variable bound in the knit environment
#'   (default "clinical_comparisons").
#' @param heading,heading_level Parent-tab heading text / markdown level (the
#'   per-comparison tabs are nested one level below).
#' @param fig_id_prefix Prefix used to build each comparison's chunk-label id.
#' @param ... Passed on to [dea_report_lines()] for every comparison (e.g.
#'   `cluster_word`, `levels`, `meta`, `interactive_plots`, `clusters_to_show`,
#'   `groups_to_show`, `gsea_metric`).
#' @returns A character vector of Quarto/knitr lines, or `character()` when no
#'   comparison is shown (e.g. all filtered out).
#' @export
dea_comparisons_lines <- function(x,
                                  obj_name = "clinical_comparisons",
                                  heading = "DEA Comparisons",
                                  heading_level = 2L,
                                  fig_id_prefix = "cmp",
                                  ...) {
  x_list <- normalize_dea_list(x)
  nms <- names(x_list)
  if (length(x_list) == 0L || is.null(nms) || any(!nzchar(nms))) {
    stop("`x` must be a named list of scitargets_dea objects.", call. = FALSE)
  }

  inner <- character()
  for (nm in nms) {
    inner <- c(inner, dea_report_lines(
      x_list[[nm]],
      fig_id = .dea_sanitize_id(fig_id_prefix, nm),
      obj_name = sprintf("%s[[%s]]", obj_name, shQuote(nm)),
      heading_level = heading_level + 1L,
      ...
    ))
  }

  # Nothing to show (e.g. every comparison filtered out): omit the whole tab.
  if (length(inner) == 0L) {
    return(character())
  }

  h <- paste(rep("#", heading_level), collapse = "")
  c(paste(h, heading), "", "::: {.panel-tabset}", "", inner, ":::", "")
}
