# Helper functions for the scitargets package.
# (New home for helpers; more will be moved here over time.)

# Internal: render a data.frame as an interactive DT2 table with download
# buttons (copy / csv / excel / pdf / print). NOT exported -- for use by
# scitargets functions only. Hand-written Quarto chunks should call DT2
# functions (DT2::dt2()) directly instead.
#
# @param data A data.frame (or object coercible with as.data.frame()).
# @param page_length Rows shown per page (default 10).
# @param buttons A list of DT2 button ids (default copy/csv/excel/pdf/print).
# @param ... Further arguments passed to DT2::dt2().
# @returns A DT2 htmlwidget.
# @keywords internal
# @noRd
.scitargets_dt2_tbl <- function(data,
                                page_length = 10L,
                                buttons = list("copy", "csv", "excel", "pdf", "print"),
                                ...) {
  DT2::dt2(
    as.data.frame(data),
    extensions = "Buttons",
    options = list(
      pageLength = page_length,
      layout = list(topEnd = list(buttons = buttons))
    ),
    ...
  )
}


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.go_ontologies <- c("BP", "CC", "MF")

# Map a species ("human" / "mouse") to its annotation DB + msigdbr species name.
# Only human and mouse are supported.
.species_db <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  switch(species,
    human = list(organism = "Human", org_db = "org.Hs.eg.db", msigdbr = "Homo sapiens"),
    mouse = list(organism = "Mouse", org_db = "org.Mm.eg.db", msigdbr = "Mus musculus")
  )
}

.default_go_params <- function(go_params = list()) {
  defaults <- list(
    organism = "Human",
    org_db = "org.Hs.eg.db",
    ontology = "BP",
    statistic = "fisher",
    algorithm = "weight01",
    top_nodes = 50L,
    num_char = 100L,
    p_adj_cutoff = 0.05,
    min_signif_genes = 5L,
    plot_number = 20L,
    node_size = 10L,
    min_foreground_genes = 1L,
    cnet_layout_seed = 5114L,
    keep_topgo_data = FALSE
  )
  params <- utils::modifyList(defaults, go_params)

  # ontology may be one or several of BP / CC / MF (BP is the default).
  params$ontology <- unique(as.character(params$ontology))
  invalid <- setdiff(params$ontology, .go_ontologies)
  if (length(invalid) > 0L) {
    stop(
      "go_params$ontology must be one or more of 'BP', 'CC', 'MF'. Invalid: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(params$ontology) == 0L) {
    params$ontology <- "BP"
  }
  params
}

.html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

.safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
}

.logfc_column <- function(markers) {
  if ("avg_log2FC" %in% colnames(markers)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers)) {
    "avg_logFC"
  } else {
    stop("The marker table must contain either 'avg_log2FC' or 'avg_logFC'.", call. = FALSE)
  }
}

.standardize_markers <- function(markers) {
  markers <- as.data.frame(markers)

  if (!"Gene" %in% colnames(markers)) {
    markers <- tibble::rownames_to_column(markers, var = "Gene")
  }

  required <- c("Gene", "p_val_adj")
  missing <- setdiff(required, colnames(markers))
  if (length(missing) > 0L) {
    stop(
      "The marker table is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  .logfc_column(markers)
  markers
}

# Default accessors --------------------------------------------------------

.dea_levels <- function(x) {
  lv <- x@levels
  if (length(lv) == 0L) lv <- names(x@de)
  lv
}

.dea_default_level <- function(x, level = NULL) {
  lv <- .dea_levels(x)
  if (is.null(level)) {
    if (length(lv) == 0L) {
      return("single_cell")
    }
    return(lv[[1L]])
  }
  if (!level %in% lv) {
    stop(
      "Level '", level, "' is not available for this comparison. Available: ",
      paste(lv, collapse = ", "),
      call. = FALSE
    )
  }
  level
}

.dea_ontologies <- function(x) {
  onto <- x@go_params$ontology %||% "BP"
  unique(as.character(onto))
}

.dea_default_ontology <- function(x, ontology = NULL) {
  onto <- .dea_ontologies(x)
  if (is.null(ontology)) {
    return(onto[[1L]])
  }
  if (!ontology %in% onto) {
    stop(
      "Ontology '", ontology, "' was not computed for this comparison. Available: ",
      paste(onto, collapse = ", "),
      call. = FALSE
    )
  }
  ontology
}

.dea_go_result <- function(x, level, direction, ontology) {
  level <- .dea_default_level(x, level)
  ontology <- .dea_default_ontology(x, ontology)
  go_level <- x@go[[level]]
  if (is.null(go_level)) {
    return(NULL)
  }
  go_dir <- go_level[[direction]]
  if (is.null(go_dir)) {
    return(NULL)
  }
  go_dir[[ontology]]
}

.empty_go_result <- function(reason = NA_character_, direction = NA_character_, params = list()) {
  list(
    status = "not_computed",
    reason = reason,
    direction = direction,
    foreground_genes = character(),
    background_genes = character(),
    params = params,
    res.table = data.frame(),
    plot_data = data.frame(),
    cnetplot.data = NULL,
    GenesAssociatedToGO = list(
      CHARACTER = list(),
      HTML = paste0("<p>", .html_escape(reason %||% "No GO enrichment result."), "</p>")
    ),
    GOdata = NULL
  )
}

# -----------------------------------------------------------------------------
# GO computation
# -----------------------------------------------------------------------------

.prepare_go_plot_data <- function(res_table, min_signif_genes, plot_number) {
  if (is.null(res_table) || nrow(res_table) == 0L) {
    return(data.frame())
  }

  gop <- as.data.frame(res_table)
  gop$Significant <- suppressWarnings(as.numeric(gop$Significant))
  gop$adjpval <- suppressWarnings(as.numeric(gop$adjpval))

  if (nrow(gop) == 0L) {
    return(data.frame())
  }

  if (max(gop$Significant, na.rm = TRUE) >= min_signif_genes) {
    gop <- gop[gop$Significant >= min_signif_genes, , drop = FALSE]
  }

  gop <- gop[order(gop$adjpval, decreasing = FALSE), , drop = FALSE]
  if (nrow(gop) > plot_number) {
    gop <- gop[seq_len(plot_number), , drop = FALSE]
  }

  gop$Term_wrapped <- stringr::str_wrap(gop$Term, width = 35)
  gop$Term_wrapped <- factor(gop$Term_wrapped, levels = rev(unique(gop$Term_wrapped)))
  gop$neg_log10_adjpval <- -log10(pmax(gop$adjpval, .Machine$double.xmin))
  gop
}

.build_go_genes_html <- function(res_table, plot_data) {
  if (is.null(res_table) || nrow(res_table) == 0L || !"Genes" %in% colnames(res_table)) {
    return("<p>No significant GO terms were detected.</p>")
  }

  if (!is.null(plot_data) && nrow(plot_data) > 0L) {
    res_table <- res_table[res_table$GO.ID %in% plot_data$GO.ID, , drop = FALSE]
  }

  if (nrow(res_table) == 0L) {
    return("<p>No significant GO terms were retained for plotting.</p>")
  }

  lines <- character()
  for (i in seq_len(nrow(res_table))) {
    lines <- c(
      lines,
      paste0(
        "<strong>Term ", .html_escape(res_table$GO.ID[[i]]), " ",
        .html_escape(res_table$Term[[i]]), ", ",
        .html_escape(res_table$Significant[[i]]),
        " genes:&nbsp;</strong>",
        .html_escape(res_table$Genes[[i]])
      )
    )
  }
  paste(lines, collapse = "<br/>")
}

.build_cnet_data <- function(res_table, min_signif_genes, plot_number, layout_seed = 5114L) {
  if (
    is.null(res_table) ||
      nrow(res_table) == 0L ||
      !all(c("GO.ID", "Term", "adjpval", "Significant", "Genes") %in% colnames(res_table))
  ) {
    return(NULL)
  }

  cnet_tab <- as.data.frame(res_table)
  cnet_tab$adjpval <- suppressWarnings(as.numeric(cnet_tab$adjpval))
  cnet_tab$Significant <- suppressWarnings(as.numeric(cnet_tab$Significant))
  cnet_tab <- cnet_tab[!is.na(cnet_tab$Genes) & cnet_tab$Genes != "", , drop = FALSE]

  if (nrow(cnet_tab) == 0L) {
    return(NULL)
  }

  if (max(cnet_tab$Significant, na.rm = TRUE) >= min_signif_genes) {
    cnet_tab <- cnet_tab[cnet_tab$Significant >= min_signif_genes, , drop = FALSE]
  }

  cnet_tab <- cnet_tab[order(cnet_tab$adjpval, decreasing = FALSE), , drop = FALSE]
  if (nrow(cnet_tab) > plot_number) {
    cnet_tab <- cnet_tab[seq_len(plot_number), , drop = FALSE]
  }

  if (nrow(cnet_tab) == 0L) {
    return(NULL)
  }

  cnet_tab$GO_label <- paste0(stringr::str_wrap(cnet_tab$Term, width = 35), "\n", cnet_tab$GO.ID)

  cnet_edges <- cnet_tab |>
    dplyr::select(GO.ID, Term, GO_label, adjpval, Significant, Genes) |>
    tidyr::separate_rows(Genes, sep = ",\\s*") |>
    dplyr::mutate(Genes = stringr::str_trim(Genes)) |>
    dplyr::filter(!is.na(Genes), Genes != "") |>
    dplyr::transmute(
      from = GO_label,
      to = Genes,
      GO.ID = GO.ID,
      Term = Term,
      adjpval = adjpval,
      Significant = Significant
    )

  if (nrow(cnet_edges) == 0L) {
    return(NULL)
  }

  go_nodes <- cnet_edges |>
    dplyr::distinct(name = from, GO.ID, Term, adjpval, Significant) |>
    dplyr::mutate(
      node_type = "GO term",
      short_label = stringr::str_wrap(Term, width = 25),
      node_size = pmax(-log10(pmax(adjpval, .Machine$double.xmin)), 1),
      tooltip = paste0(
        "GO term: ", Term,
        "\nGO ID: ", GO.ID,
        "\nAdj. p-value: ", signif(adjpval, 3),
        "\nSignificant genes: ", Significant
      )
    )

  gene_nodes <- cnet_edges |>
    dplyr::distinct(name = to) |>
    dplyr::mutate(
      GO.ID = NA_character_,
      Term = name,
      adjpval = NA_real_,
      Significant = NA_real_,
      node_type = "Gene",
      short_label = NA_character_,
      node_size = 0.7,
      tooltip = paste0("Gene: ", name)
    )

  cnet_nodes <- dplyr::bind_rows(go_nodes, gene_nodes)

  cnet_graph <- igraph::graph_from_data_frame(
    d = dplyr::select(cnet_edges, from, to),
    vertices = cnet_nodes,
    directed = FALSE
  )

  set.seed(layout_seed)
  coords <- igraph::layout_with_fr(cnet_graph)

  cnet_nodes <- as.data.frame(cnet_nodes, stringsAsFactors = FALSE)
  cnet_nodes$xpos <- as.numeric(coords[, 1L])
  cnet_nodes$ypos <- as.numeric(coords[, 2L])

  from_pos <- cnet_nodes |>
    dplyr::transmute(from = name, x = xpos, y = ypos)

  to_pos <- cnet_nodes |>
    dplyr::transmute(to = name, xend = xpos, yend = ypos)

  edge_df <- cnet_edges |>
    dplyr::left_join(from_pos, by = "from") |>
    dplyr::left_join(to_pos, by = "to")

  list(
    edges = edge_df,
    nodes = cnet_nodes,
    graph = cnet_graph
  )
}

# Single-ontology topGO run (with its own BH adjustment).
.compute_topgo_result <- function(fg.genes,
                                  bg.genes,
                                  direction = NA_character_,
                                  go_params = list()) {
  params <- .default_go_params(go_params)

  # topGO's S4 class "topGOdata" and the org.Hs.eg.db annotation must be
  # registered before methods::new() is called. Load the namespaces explicitly
  # so the function works whether or not these packages are attached.
  if (!params$org_db %in% c("org.Hs.eg.db", "org.Mm.eg.db")) {
    stop("Only human (org.Hs.eg.db) and mouse (org.Mm.eg.db) are supported for GO.", call. = FALSE)
  }
  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("Package 'topGO' is required for GO enrichment.", call. = FALSE)
  }
  # org.Hs.eg.db / org.Mm.eg.db are in Suggests, so ATTACH the needed one with
  # require() (topGO::annFUN.org resolves the annotation package by name and needs
  # it on the search path). Clear error if the Suggests package is not installed.
  if (!require(params$org_db, character.only = TRUE, quietly = TRUE)) {
    stop("Package '", params$org_db, "' (", params$organism,
      ") is required for GO enrichment; it is in Suggests -- please install it.",
      call. = FALSE
    )
  }
  # A single topGO run uses a single ontology; callers loop over ontologies.
  ontology <- params$ontology[[1L]]
  if (!ontology %in% .go_ontologies) {
    stop("GO ontology must be one of 'BP', 'CC', or 'MF'.", call. = FALSE)
  }

  bg.genes <- unique(stats::na.omit(as.character(bg.genes)))
  fg.genes <- unique(stats::na.omit(as.character(fg.genes)))
  fg.genes <- intersect(fg.genes, bg.genes)

  if (length(bg.genes) == 0L) {
    return(.empty_go_result("No background genes were supplied.", direction, params))
  }
  if (length(fg.genes) < params$min_foreground_genes) {
    return(.empty_go_result("No foreground genes passed the differential-expression filters.", direction, params))
  }

  gene_list <- integer(length(bg.genes))
  names(gene_list) <- bg.genes
  gene_list[intersect(names(gene_list), fg.genes)] <- 1L
  gene_list <- factor(gene_list)

  GOdata <- methods::new(
    "topGOdata",
    description = paste0("GO analysis: ", direction, " (", ontology, ")"),
    ontology = ontology,
    allGenes = gene_list,
    annot = topGO::annFUN.org,
    mapping = params$org_db,
    ID = "SYMBOL",
    nodeSize = params$node_size
  )

  res_result <- topGO::runTest(
    GOdata,
    statistic = params$statistic,
    algorithm = params$algorithm
  )

  # All scored GO terms define the multiple-testing universe.
  # Do NOT limit this to params$top_nodes before BH correction.
  all_scores <- topGO::score(res_result)
  all_scores <- all_scores[!is.na(all_scores)]

  n_scored <- length(all_scores)

  if (n_scored < 1L) {
    out <- .empty_go_result("topGO returned no GO terms.", direction, params)
    out$foreground_genes <- fg.genes
    out$background_genes <- bg.genes
    if (isTRUE(params$keep_topgo_data)) out$GOdata <- GOdata
    return(out)
  }

  # Request all scored terms from GenTable so that the adjusted p-values are
  # computed over all tested terms. The final table is capped to params$top_nodes
  # only after p-value adjustment and significance filtering.
  res_table <- topGO::GenTable(
    GOdata,
    pval = res_result,
    topNodes = n_scored,
    numChar = params$num_char
  )

  if (is.null(res_table) || nrow(res_table) == 0L) {
    out <- .empty_go_result("topGO returned no GO terms.", direction, params)
    out$foreground_genes <- fg.genes
    out$background_genes <- bg.genes
    if (isTRUE(params$keep_topgo_data)) out$GOdata <- GOdata
    return(out)
  }

  res_table <- as.data.frame(res_table, stringsAsFactors = FALSE)

  # Use numeric topGO scores, keyed by GO.ID.
  # This avoids relying on the formatted p-value column returned by GenTable().
  res_table$pval <- unname(all_scores[res_table$GO.ID])

  # BH adjustment is done individually within this ontology, but across ALL
  # scored GO terms for that ontology, not only the displayed top_nodes terms.
  res_table$adjpval <- stats::p.adjust(res_table$pval, method = "BH")

  res_table <- res_table[!is.na(res_table$adjpval) & res_table$adjpval <= params$p_adj_cutoff, , drop = FALSE]

  if (nrow(res_table) == 0L) {
    out <- .empty_go_result("No GO terms passed the adjusted p-value threshold.", direction, params)
    out$foreground_genes <- fg.genes
    out$background_genes <- bg.genes
    if (isTRUE(params$keep_topgo_data)) out$GOdata <- GOdata
    return(out)
  }

  # Deterministic ordering before truncating the final output.
  res_table <- res_table[
    order(res_table$adjpval, res_table$pval, res_table$GO.ID), ,
    drop = FALSE
  ]

  # The final GO table should not show more than params$top_nodes rows.
  top_nodes <- min(params$top_nodes, nrow(res_table))
  res_table <- res_table[seq_len(top_nodes), , drop = FALSE]

  res_table$Annotated <- suppressWarnings(as.numeric(res_table$Annotated))
  res_table$Significant <- suppressWarnings(as.numeric(res_table$Significant))
  res_table$Expected <- suppressWarnings(as.numeric(res_table$Expected))
  res_table$Enrichment <- format(
    round(
      .safe_divide(res_table$Significant, length(fg.genes)) /
        .safe_divide(res_table$Annotated, length(bg.genes)),
      2L
    ),
    nsmall = 2L
  )
  res_table$Genes <- NA_character_

  myterms <- res_table$GO.ID
  term_genes <- topGO::genesInTerm(GOdata, myterms)
  sig_genes <- topGO::sigGenes(GOdata)
  genes_associated <- list()

  for (i in seq_along(myterms)) {
    term_id <- myterms[[i]]
    term_name <- res_table$Term[[i]]
    genes_for_term <- unique(intersect(sig_genes, term_genes[[term_id]]))
    genes_associated[[as.character(term_name)]] <- genes_for_term
    res_table$Genes[[i]] <- paste(genes_for_term, collapse = ", ")
  }

  rownames(res_table) <- NULL
  res_table <- dplyr::relocate(res_table, Genes, .after = Enrichment)

  plot_data <- .prepare_go_plot_data(
    res_table = res_table,
    min_signif_genes = params$min_signif_genes,
    plot_number = params$plot_number
  )

  cnet_data <- .build_cnet_data(
    res_table = res_table,
    min_signif_genes = params$min_signif_genes,
    plot_number = params$plot_number,
    layout_seed = params$cnet_layout_seed
  )

  list(
    status = "computed",
    reason = NA_character_,
    direction = direction,
    ontology = ontology,
    foreground_genes = fg.genes,
    background_genes = bg.genes,
    params = params,
    res.table = res_table,
    plot_data = plot_data,
    cnetplot.data = cnet_data,
    GenesAssociatedToGO = list(
      CHARACTER = genes_associated,
      HTML = .build_go_genes_html(res_table, plot_data)
    ),
    GOdata = if (isTRUE(params$keep_topgo_data)) GOdata else NULL
  )
}

# Run topGO for every requested ontology and return a named list (by ontology).
# The BH adjustment happens inside each per-ontology run, so it is done
# individually for each ontology, as required.
.compute_go_results <- function(fg.genes,
                                bg.genes,
                                direction = NA_character_,
                                go_params = list()) {
  params <- .default_go_params(go_params)
  ontologies <- params$ontology

  res <- list()
  for (onto in ontologies) {
    onto_params <- params
    onto_params$ontology <- onto
    res[[onto]] <- suppressMessages(.compute_topgo_result(
      fg.genes = fg.genes,
      bg.genes = bg.genes,
      direction = direction,
      go_params = onto_params
    ))
  }
  res
}

# -----------------------------------------------------------------------------
# Differential expression helpers (used by the run_dea dispatcher below)
# -----------------------------------------------------------------------------

# Single-cell differential expression with Seurat::FindMarkers.
# `group_col` is the metadata column holding the two idents to compare.
.de_single_cell <- function(seurat_obj, group_col, ident.1, ident.2,
                            min_pct = 0.1, logfc_threshold = 0.25, test_use = "wilcox",
                            only_pos = FALSE,
                            assay = "SCT", seed = 5114L, recorrect_umi = TRUE) {
  markers <- Seurat::FindMarkers(
    seurat_obj,
    assay = assay,
    test.use = test_use,
    group.by = group_col,
    ident.1 = ident.1,
    ident.2 = ident.2,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    random.seed = seed,
    recorrect_umi = recorrect_umi
  )
  .standardize_markers(markers)
}

# Run GO (both directions, all ontologies) for one level's marker table.
# `fg_cutoff` is the adjusted-p cutoff used to pick the foreground (DE) genes for
# GO. It defaults to params$p_adj_cutoff but run_dea passes the per-level
# padj cutoff (P1.8). It is distinct from params$p_adj_cutoff, which downstream
# (.compute_go_results) uses to filter the GO TERMS by their adjusted p-value.
.go_for_markers <- function(markers, params, fg_cutoff = params$p_adj_cutoff) {
  logfc_col <- .logfc_column(markers)
  background_genes <- markers$Gene
  up_genes <- markers$Gene[
    !is.na(markers$p_val_adj) &
      markers$p_val_adj <= fg_cutoff &
      !is.na(markers[[logfc_col]]) &
      markers[[logfc_col]] > 0
  ]
  down_genes <- markers$Gene[
    !is.na(markers$p_val_adj) &
      markers$p_val_adj <= fg_cutoff &
      !is.na(markers[[logfc_col]]) &
      markers[[logfc_col]] < 0
  ]

  list(
    up = .compute_go_results(up_genes, background_genes, "up", params),
    down = .compute_go_results(down_genes, background_genes, "down", params)
  )
}

# Empty GO container (used when a comparison is not computed), one entry per
# ontology and direction.
.empty_go_levels <- function(reason, params) {
  onto <- params$ontology
  make_dir <- function(direction) {
    stats::setNames(
      lapply(onto, function(o) {
        p <- params
        p$ontology <- o
        .empty_go_result(reason, direction, p)
      }),
      onto
    )
  }
  list(up = make_dir("up"), down = make_dir("down"))
}


# -----------------------------------------------------------------------------
# xlsx report (merged from the former ExportMarkersTables.R)
# -----------------------------------------------------------------------------

.safe_filename_component <- function(x) {
  x <- as.character(x)

  # Keep group names readable, but avoid path separators and invalid filename chars.
  x <- gsub("/", "-", x, fixed = TRUE)
  x <- gsub("\\\\", "-", x)
  x <- gsub('[<>:"|?*]+', "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("-+", "-", x)
  x <- gsub("^[_-]+|[_-]+$", "", x)

  if (!nzchar(x)) {
    x <- "unknown"
  }

  x
}

.dea_xlsx_filename <- function(x) {
  group1 <- .safe_filename_component(x@group1)
  group2 <- .safe_filename_component(x@group2)
  cluster <- .safe_filename_component(x@cluster)

  paste0(group1, "_vs_", group2, "_cluster_", cluster, ".xlsx")
}

# Excel sheet names are capped at 31 characters and must be unique.
.sheet_name <- function(...) {
  nm <- paste(c(...), collapse = "_")
  # Excel forbids \ / ? * [ ] : in sheet names; replace anything that is not a
  # letter, digit, underscore or hyphen, and cap at 31 characters.
  nm <- gsub("[^A-Za-z0-9_-]+", "_", nm)
  if (nchar(nm) > 31L) nm <- substr(nm, 1L, 31L)
  nm
}

.go_input_summary_df <- function(x, level) {
  go_level <- x@go[[level]]
  onto <- .dea_ontologies(x)

  rows <- list()
  for (direction in c("up", "down")) {
    for (o in onto) {
      res <- if (is.null(go_level)) NULL else go_level[[direction]][[o]]
      rows[[length(rows) + 1L]] <- data.frame(
        direction = direction,
        ontology = o,
        go_status = if (is.null(res)) NA_character_ else res$status %||% NA_character_,
        reason = if (is.null(res)) NA_character_ else res$reason %||% NA_character_,
        foreground_genes = if (is.null(res)) 0L else length(res$foreground_genes),
        background_genes = if (is.null(res)) 0L else length(res$background_genes),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.params_to_df <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(data.frame(parameter = character(), value = character()))
  }

  data.frame(
    parameter = names(x),
    value = vapply(
      x,
      function(value) paste(as.character(value), collapse = ", "),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

.go_empty_message <- function(reason, direction, ontology) {
  if (is.null(reason) || is.na(reason) || !nzchar(reason)) {
    reason <- paste0(
      "No significant GO:", ontology, " terms were detected for ",
      direction, "-regulated genes."
    )
  }
  data.frame(message = reason, stringsAsFactors = FALSE)
}

.write_table_sheet <- function(wb, sheet, data, empty_data = NULL) {
  openxlsx::addWorksheet(wb, sheet)

  data <- as.data.frame(data)
  if (nrow(data) == 0L && !is.null(empty_data)) {
    data <- as.data.frame(empty_data)
  }

  openxlsx::writeData(wb, sheet = sheet, x = data)

  if (nrow(data) > 0L) {
    openxlsx::freezePane(wb, sheet = sheet, firstRow = TRUE)
  }

  if (ncol(data) > 0L) {
    openxlsx::setColWidths(
      wb,
      sheet = sheet,
      cols = seq_len(ncol(data)),
      widths = "auto"
    )
  }

  invisible(NULL)
}

.write_summary_block <- function(wb, sheet, title, data, start_row, title_style) {
  openxlsx::writeData(
    wb,
    sheet = sheet,
    x = title,
    startRow = start_row,
    startCol = 1,
    colNames = FALSE
  )

  openxlsx::addStyle(
    wb,
    sheet = sheet,
    style = title_style,
    rows = start_row,
    cols = 1
  )

  data <- as.data.frame(data)

  openxlsx::writeData(
    wb,
    sheet = sheet,
    x = data,
    startRow = start_row + 1L,
    startCol = 1
  )

  start_row + nrow(data) + 3L
}

# One-block data.frame documenting the pseudobulk DESeq2 statistics for the xlsx
# summary sheet (P2.1g). DESeq2::results adjusts p-values with Benjamini-Hochberg
# (its default pAdjustMethod = "BH").
.pb_stats_df <- function(x) {
  dp <- x@de_params
  pb_test <- dp$pb_test %||% "Wald"
  lcf <- dp$pb_low_count_filter %||% c(10L, 3L)
  shrink_requested <- isTRUE(dp$pb_lfc_shrink)
  shrink_method <- dp$pb_lfc_shrink_method %||%
    attr(x@de[["pseudobulk"]], "lfcShrink_method")
  shrink_val <- if (shrink_requested && !is.null(shrink_method) &&
    !shrink_method %in% c("none", "disabled")) {
    paste0("yes (", shrink_method, ")")
  } else if (shrink_requested) {
    "yes"
  } else {
    "no"
  }
  test_val <- paste0(
    pb_test,
    if (identical(pb_test, "LRT")) " (likelihood-ratio test)" else " (two-sided)"
  )
  data.frame(
    field = c(
      "DE engine",
      "normalization",
      "statistical test",
      "p-value adjustment",
      "design formula",
      if (identical(pb_test, "LRT")) "reduced formula" else NULL,
      "low-count gene filter",
      "min replicates per group",
      "log2FC shrinkage",
      "fold-change columns"
    ),
    value = c(
      "DESeq2",
      "median-of-ratios (DESeq2 size factors)",
      test_val,
      "Benjamini-Hochberg (FDR)",
      dp$pb_design %||% "~ group",
      if (identical(pb_test, "LRT")) (dp$pb_reduced %||% "~ 1") else NULL,
      sprintf("count >= %s in >= %s samples", lcf[1], lcf[2]),
      as.character(dp$min_replicates %||% 3L),
      shrink_val,
      "avg_log2FC = MLE; avg_log2FC_shrink = shrunken"
    ),
    stringsAsFactors = FALSE
  )
}

.write_summary_sheet <- function(wb, x) {
  sheet <- "summary"
  openxlsx::addWorksheet(wb, sheet)

  title_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fontSize = 13
  )

  comparison_df <- data.frame(
    field = c(
      "group1",
      "group2",
      "group_by",
      "cluster_by",
      "cluster",
      "levels",
      "status",
      "reason"
    ),
    value = c(
      x@group1,
      x@group2,
      x@group_by,
      paste(x@cluster_by, collapse = ", "),
      x@cluster,
      paste(x@levels, collapse = ", "),
      x@status,
      x@reason
    ),
    stringsAsFactors = FALSE
  )

  summary_blocks <- list(
    "Comparison summary" = comparison_df,
    "Cell counts" = as.data.frame(x@n_cells),
    "Differential-expression parameters" = .params_to_df(x@de_params),
    "GO parameters" = .params_to_df(x@go_params),
    "GSEA parameters" = .params_to_df(x@gsea_params)
  )

  # Pseudobulk DESeq2 statistics block (P2.1g): document the test + p-value
  # adjustment so the pseudobulk xlsx is self-describing.
  if ("pseudobulk" %in% x@levels) {
    summary_blocks <- append(
      summary_blocks,
      list("Pseudobulk DESeq2 statistics" = .pb_stats_df(x)),
      after = which(names(summary_blocks) == "Differential-expression parameters")
    )
  }

  for (level in .dea_levels(x)) {
    summary_blocks[[paste0("GO input gene sets (", level, ")")]] <-
      .go_input_summary_df(x, level)
  }

  row <- 1L
  for (block_name in names(summary_blocks)) {
    row <- .write_summary_block(
      wb = wb,
      sheet = sheet,
      title = block_name,
      data = summary_blocks[[block_name]],
      start_row = row,
      title_style = title_style
    )
  }

  openxlsx::freezePane(wb, sheet = sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet, cols = 1:10, widths = "auto")

  invisible(NULL)
}

.extract_single_dea <- function(x) {
  x_list <- normalize_dea_list(x)

  if (length(x_list) != 1L) {
    stop(
      "Expected exactly one scitargets_dea object per branch, got ",
      length(x_list),
      ".",
      call. = FALSE
    )
  }

  x <- x_list[[1L]]

  if (!is_scitargets_dea(x)) {
    stop("Input is not a scitargets_dea object.", call. = FALSE)
  }

  x
}


# Restrict a set of available cluster levels to the user-requested ones.
#   * NA (default) -> keep every cluster present in the data
#   * a vector     -> keep only those clusters (e.g. CD3+ T-cell celltypes),
#                     erroring on names that are not present.
.resolve_clusters <- function(cluster_levels, clusters = NA) {
  cluster_levels <- unique(as.character(stats::na.omit(cluster_levels)))
  if (length(clusters) == 1L && (is.na(clusters) || identical(clusters, NA))) {
    return(cluster_levels)
  }
  clusters <- as.character(clusters)
  invalid <- setdiff(clusters, cluster_levels)
  if (length(invalid) > 0L) {
    stop(
      "Invalid cluster name(s): ", paste(invalid, collapse = ", "),
      ". Available clusters: ", paste(cluster_levels, collapse = ", "),
      call. = FALSE
    )
  }
  intersect(cluster_levels, clusters)
}

# Standardize a DESeq2 results object into the common marker schema used by the
# GO / GSEA / plotting code (Gene, avg_log2FC, p_val, p_val_adj, ...).
.standardize_deseq2 <- function(res) {
  df <- as.data.frame(res)
  df <- tibble::rownames_to_column(df, var = "Gene")

  out <- data.frame(
    Gene = df$Gene,
    avg_log2FC = df$log2FoldChange,
    p_val = df$pvalue,
    p_val_adj = df$padj,
    baseMean = df$baseMean,
    lfcSE = df$lfcSE,
    stat = df$stat,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$p_val_adj, out$p_val, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Make a string safe to use as a DESeq2 factor level / column name (letters,
# numbers, '_' and '.' only), avoiding the DESeq2 "characters other than ..."
# note for group names like "CR1M/6M". Used only internally for the design; the
# original group labels are kept everywhere else (object, plots, tables, text).
.dea_safe_level <- function(x) {
  x <- gsub("[^A-Za-z0-9_.]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_+|_+$", "", x)
}

# Label each cell with its comparison group WITHOUT subsetting the object: cells
# in the chosen cluster that belong to group1/group2 get that group label (or
# "rest" when group2 == "rest"); every other cell is NA. Returns a named vector
# aligned to rownames(meta). Replaces the old .subset_cluster_groups (subset()
# broke the shared SCT model state FindMarkers needs); FindMarkers selects the
# two groups via group.by/ident and PseudobulkExpression drops the NA cells.
.dea_group_labels <- function(meta, group_by, cluster_by, cluster, group1, group2) {
  keep_cluster <- if (length(cluster_by) == 1L && nzchar(cluster_by)) {
    as.character(meta[[cluster_by]]) == cluster & !is.na(meta[[cluster_by]])
  } else {
    rep(TRUE, nrow(meta))
  }
  grp <- as.character(meta[[group_by]])
  labels <- if (identical(group2, "rest")) {
    ifelse(keep_cluster & !is.na(grp),
      ifelse(grp == group1, group1, "rest"), NA_character_
    )
  } else {
    ifelse(keep_cluster & grp %in% c(group1, group2), grp, NA_character_)
  }
  stats::setNames(labels, rownames(meta))
}

# Sentinel returned when a differential-expression level cannot be computed.
# (R does not allow setting attributes on NULL, so we use a tagged empty object
# that can carry the reason.)
.dea_failed <- function(reason) {
  structure(list(), class = "dea_failed", reason = reason)
}
.dea_is_failed <- function(x) inherits(x, "dea_failed")

# Shrink the group1-vs-group2 log2 fold change. Cascade: apeglm (preferred; needs
# the contrast to be a model coefficient -- it is, since `group`'s reference is set
# to group2) -> ashr (accepts an arbitrary contrast) -> DESeq2's built-in normal
# prior (always available). Returns list(lfc = named numeric, method = character).
.dea_lfc_shrink <- function(dds, s1, s2) {
  coef_name <- paste0("group_", s1, "_vs_", s2)
  rn <- DESeq2::resultsNames(dds)
  as_named <- function(res) stats::setNames(res$log2FoldChange, rownames(res))

  # apeglm (preferred) -- only when the contrast is a coefficient and pkg present.
  if (coef_name %in% rn && requireNamespace("apeglm", quietly = TRUE)) {
    res <- tryCatch(
      DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      return(list(lfc = as_named(res), method = "apeglm"))
    }
  }
  # ashr (fallback) -- works with an arbitrary contrast.
  if (requireNamespace("ashr", quietly = TRUE)) {
    res <- tryCatch(
      DESeq2::lfcShrink(dds, contrast = c("group", s1, s2), type = "ashr", quiet = TRUE),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      return(list(lfc = as_named(res), method = "ashr"))
    }
  }
  # normal (DESeq2 built-in) -- no extra package required.
  res <- tryCatch(
    if (coef_name %in% rn) {
      DESeq2::lfcShrink(dds, coef = coef_name, type = "normal", quiet = TRUE)
    } else {
      DESeq2::lfcShrink(dds, contrast = c("group", s1, s2), type = "normal", quiet = TRUE)
    },
    error = function(e) NULL
  )
  if (!is.null(res)) {
    return(list(lfc = as_named(res), method = "normal"))
  }

  # all methods failed -> NA shrinkage (MLE remains in avg_log2FC).
  list(lfc = stats::setNames(rep(NA_real_, nrow(dds)), rownames(dds)), method = "none")
}

# Fold-change column to rank DE genes by: for pseudobulk use the shrunken LFC
# (avg_log2FC_shrink) when shrinkage was actually applied, otherwise the MLE
# avg_log2FC (also the single-cell case).
.dea_top_de_fc_col <- function(obj, level, de) {
  if (identical(level, "pseudobulk")) {
    used <- isTRUE(obj@de_params$pb_lfc_shrink)
    m <- obj@de_params$pb_lfc_shrink_method
    if (used && !is.null(m) && !m %in% c("none", "disabled") &&
      "avg_log2FC_shrink" %in% names(de)) {
      return("avg_log2FC_shrink")
    }
  }
  "avg_log2FC"
}

# Union of the top-`n_top` DE genes (by absolute fold-change among adjusted-p
# significant genes) across the comparisons that pass the clusters_to_show /
# groups_to_show display filters (i.e. the comparisons dea_report_lines shows).
.dea_top_de_genes <- function(x, level = NULL, n_top = 20L, padj_cutoff = NULL,
                              clusters_to_show = NA, groups_to_show = NA) {
  x_list <- normalize_dea_list(x)
  genes <- character()
  for (obj in x_list) {
    if (!.dea_show_comparison(obj, clusters_to_show, groups_to_show)) next
    lv <- level %||% .dea_default_level(obj, NULL)
    de <- obj@de[[lv]]
    if (is.null(de) || nrow(de) == 0L || !"p_val_adj" %in% names(de)) next
    cutoff <- padj_cutoff %||% obj@padj_cutoffs[[lv]] %||% 0.05
    sig <- de[!is.na(de$p_val_adj) & de$p_val_adj <= cutoff, , drop = FALSE]
    if (nrow(sig) == 0L) next
    fc <- sig[[.dea_top_de_fc_col(obj, lv, sig)]]
    sig <- sig[order(abs(fc), decreasing = TRUE), , drop = FALSE]
    genes <- c(genes, sig$Gene[seq_len(min(n_top, nrow(sig)))])
  }
  unique(genes[!is.na(genes) & nzchar(genes)])
}

# Turn a user named list of Quarto chunk options into `#|` option lines for a
# generated chunk. `column` defaults to "screen-outset"; `label` and `caption`
# (fig-cap) fall back to the supplied defaults; `fig_width` / `fig_height` are
# emitted only when present in `opts`. Always includes message/warning = false.
.dea_chunk_opts <- function(opts = list(), default_label, default_caption) {
  if (is.null(opts)) opts <- list()
  lines <- c(
    "#| message: false",
    "#| warning: false",
    sprintf("#| column: %s", opts$column %||% "screen-outset"),
    sprintf("#| label: %s", opts$label %||% default_label),
    sprintf("#| fig-cap: \"%s\"", opts$caption %||% default_caption)
  )
  if (!is.null(opts$fig_width)) {
    lines <- c(lines, sprintf("#| fig-width: %s", opts$fig_width))
  }
  if (!is.null(opts$fig_height)) {
    lines <- c(lines, sprintf("#| fig-height: %s", opts$fig_height))
  }
  lines
}

# P2.14.2 display filters for dea_report_lines.
# Is NA / NULL / empty? (the "show everything" sentinel)
.dea_filter_all <- function(v) length(v) == 0L || (length(v) == 1L && is.na(v))

# Should this comparison be rendered given clusters_to_show / groups_to_show?
.dea_show_comparison <- function(x, clusters_to_show = NA, groups_to_show = NA) {
  if (!.dea_filter_all(clusters_to_show)) {
    if (!x@cluster %in% as.character(clusters_to_show)) {
      return(FALSE)
    }
  }
  if (!.dea_filter_all(groups_to_show)) {
    g <- as.character(groups_to_show)
    pair <- c(x@group1, x@group2)
    keep <- if (length(g) == 1L) {
      g %in% pair # any comparison involving that group
    } else if (length(g) == 2L) {
      setequal(g, pair) # only the comparison between the two
    } else {
      all(pair %in% g) # both groups among the listed set
    }
    if (!isTRUE(keep)) {
      return(FALSE)
    }
  }
  TRUE
}

# Invalid-name messages for the display filters, but only when the caller passed
# the universe of valid names (available_clusters / available_groups).
.dea_filter_invalid <- function(clusters_to_show, groups_to_show,
                                available_clusters = NULL, available_groups = NULL) {
  msgs <- character()
  if (!is.null(available_clusters) && !.dea_filter_all(clusters_to_show)) {
    bad <- setdiff(as.character(clusters_to_show), as.character(available_clusters))
    if (length(bad) > 0L) {
      msgs <- c(msgs, sprintf(
        "Unknown cluster(s) in `clusters_to_show`: %s. Available: %s.",
        paste(bad, collapse = ", "), paste(available_clusters, collapse = ", ")
      ))
    }
  }
  if (!is.null(available_groups) && !.dea_filter_all(groups_to_show)) {
    bad <- setdiff(as.character(groups_to_show), as.character(available_groups))
    if (length(bad) > 0L) {
      msgs <- c(msgs, sprintf(
        "Unknown group(s) in `groups_to_show`: %s. Available: %s.",
        paste(bad, collapse = ", "), paste(available_groups, collapse = ", ")
      ))
    }
  }
  msgs
}

# Build the UpSet membership data for DE genes across comparisons. Returns a
# data.frame with one row per DE gene and a list-column `comparisons` naming the
# comparisons in which the gene is DE (in the requested direction), suitable for
# ggupset::scale_x_upset(). `comparisons` arg: NA = all comparisons (that have DE
# genes), or a character vector of comparison_name(s) to restrict to.
.dea_upset_data <- function(x, direction = c("both", "up", "down"), level = NULL,
                            comparisons = NA, padj_cutoff = NULL) {
  direction <- match.arg(direction)
  x_list <- normalize_dea_list(x)
  want_all <- length(comparisons) == 1L && is.na(comparisons)
  sets <- list()
  for (obj in x_list) {
    nm <- obj@comparison_name
    if (!want_all && !(nm %in% comparisons)) next
    lv <- level %||% .dea_default_level(obj, NULL)
    de <- obj@de[[lv]]
    if (is.null(de) || nrow(de) == 0L || !"p_val_adj" %in% names(de)) next
    cutoff <- padj_cutoff %||% obj@padj_cutoffs[[lv]] %||% 0.05
    sig <- !is.na(de$p_val_adj) & de$p_val_adj <= cutoff
    lfc <- de$avg_log2FC
    genes <- switch(direction,
      up   = de$Gene[sig & !is.na(lfc) & lfc > 0],
      down = de$Gene[sig & !is.na(lfc) & lfc < 0],
      both = de$Gene[sig]
    )
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    if (length(genes) > 0L) sets[[nm]] <- genes
  }
  if (length(sets) == 0L) {
    return(data.frame())
  }
  all_genes <- sort(unique(unlist(sets, use.names = FALSE)))
  memb <- lapply(all_genes, function(g) {
    names(sets)[vapply(sets, function(s) g %in% s, logical(1L))]
  })
  df <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  df$comparisons <- memb # list-column for ggupset
  df
}

# Testability grid (state/cluster x comparison) across a list of comparisons.
# metric: "testable" (1 if the level was computed), "patients" (min #replicates per
# group, from the stored PCA coords), or "cells" (min cells per group, from @n_cells).
# `after_removal` toggles the BEFORE / AFTER outlier-removal view for the
# patients/cells/testable metrics (computed from the PCA coords so before vs after
# are directly comparable): TRUE excludes the flagged outlier replicates, FALSE
# keeps them. `outliers` (count of flagged replicates) ignores `after_removal`.
.dea_testability_data <- function(x, level = "pseudobulk",
                                  metric = c("testable", "patients", "cells", "outliers"),
                                  after_removal = TRUE) {
  metric <- match.arg(metric)
  x_list <- normalize_dea_list(x)
  if (length(x_list) == 0L) {
    return(data.frame())
  }
  rows <- lapply(x_list, function(obj) {
    ref <- if (identical(obj@group2, "rest")) "rest" else obj@group2
    co <- obj@pca[[level]]$coords
    min_rep <- obj@de_params$min_replicates %||% 3L
    # per-group patient / cell counts from coords, optionally excluding outliers
    grp_counts <- function() {
      keep <- co$comparison_role %in% c(obj@group1, ref) & (!after_removal | !co$outlier)
      sub <- co[keep, , drop = FALSE]
      list(
        n1 = length(unique(sub$unit[sub$comparison_role == obj@group1])),
        n2 = length(unique(sub$unit[sub$comparison_role == ref])),
        c1 = sum(sub$n_cells[sub$comparison_role == obj@group1], na.rm = TRUE),
        c2 = sum(sub$n_cells[sub$comparison_role == ref], na.rm = TRUE)
      )
    }
    val <- switch(metric,
      testable = {
        if (is.null(co)) {
          as.numeric(identical(obj@status, "computed") && level %in% obj@levels)
        } else {
          g <- grp_counts()
          as.numeric(min(g$n1, g$n2) >= min_rep)
        }
      },
      patients = if (is.null(co)) {
        NA_real_
      } else {
        g <- grp_counts()
        min(g$n1, g$n2)
      },
      cells = {
        if (is.null(co)) {
          nc <- as.data.frame(obj@n_cells) # fallback (pre-removal) when no PCA
          g1 <- nc$N[nc$comparison_group == obj@group1]
          g2 <- nc$N[nc$comparison_group == ref]
          if (length(g1) && length(g2)) min(g1, g2) else NA_real_
        } else {
          g <- grp_counts()
          min(g$c1, g$c2)
        }
      },
      outliers = {
        # flagged outlier replicates among the comparison's two groups (those
        # removed from DESeq2 when pb_remove_outliers = TRUE).
        if (is.null(co)) NA_real_ else sum(co$outlier & co$comparison_role %in% c(obj@group1, ref))
      }
    )
    data.frame(
      cluster = obj@cluster,
      comparison = paste(obj@group1, "vs", obj@group2),
      value = as.numeric(val), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Any flagged outlier replicates across the comparisons (drives the before/after
# tabset in dea_testability_lines).
.dea_any_outliers <- function(x, level = "pseudobulk") {
  d <- .dea_testability_data(x, level = level, metric = "outliers")
  nrow(d) > 0L && any(d$value > 0, na.rm = TRUE)
}

# Patient-level cell-state composition from a metadata data.frame. For each
# (patient x state) returns the cell count, the patient's total cells, the
# proportion of that state among the patient's cells, and the patient's clinical
# group. Cells with an NA state are dropped (e.g. non-CD3 cells in
# clinical_azimuth_l2). Proportions (not raw counts) are the biological display.
.cell_state_composition <- function(meta, state_col, group_col, patient_col) {
  md <- as.data.frame(meta)
  for (cc in c(state_col, group_col, patient_col)) {
    if (!cc %in% colnames(md)) {
      stop("Column '", cc, "' was not found in the metadata.", call. = FALSE)
    }
  }
  md <- md[!is.na(md[[state_col]]), , drop = FALSE]
  if (nrow(md) == 0L) {
    return(data.frame())
  }
  pat <- as.character(md[[patient_col]])
  st <- as.character(md[[state_col]])

  tab <- as.data.frame(table(patient = pat, state = st), stringsAsFactors = FALSE)
  tab$n_cells <- as.integer(tab$Freq)
  tab$Freq <- NULL
  totals <- tapply(tab$n_cells, tab$patient, sum)
  tab$n_total <- as.integer(totals[tab$patient])
  tab$proportion <- ifelse(tab$n_total > 0L, tab$n_cells / tab$n_total, NA_real_)

  # one clinical group per patient (first observed)
  pg <- md[!duplicated(pat), , drop = FALSE]
  tab$group <- as.character(pg[[group_col]])[match(tab$patient, as.character(pg[[patient_col]]))]
  tab <- tab[tab$n_total > 0L, , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

# Emit a user-facing message describing the pseudobulk LFC-shrinkage outcome for
# one comparison (printed by run_dea when the pseudobulk level is computed).
.dea_message_shrink <- function(comparison_name, method, requested) {
  if (!isTRUE(requested)) {
    message(sprintf(
      "[run_dea] %s: pseudobulk LFC shrinkage disabled (avg_log2FC_shrink = NA).",
      comparison_name
    ))
    return(invisible(NULL))
  }
  desc <- switch(method %||% "none",
    apeglm = "apeglm (DESeq2::lfcShrink type='apeglm'; adaptive Bayesian shrinkage on the group coefficient)",
    ashr   = "ashr (DESeq2::lfcShrink type='ashr'; adaptive shrinkage on the group1-vs-group2 contrast)",
    normal = "normal (DESeq2 built-in lfcShrink type='normal'; apeglm/ashr not installed)",
    none   = "FAILED -- no shrinkage method succeeded (avg_log2FC_shrink = NA)",
    method
  )
  message(sprintf(
    "[run_dea] %s: pseudobulk LFC shrinkage = %s. MLE kept in avg_log2FC; shrunken in avg_log2FC_shrink.",
    comparison_name, desc
  ))
  invisible(NULL)
}

# Add design covariates to the pseudobulk colData. When `design` references
# variables beyond `group` (e.g. "~ origin + group"), look up one value per
# replicate from the cell metadata, taking the first row per `covariate_key`
# (dplyr::filter(!duplicated(key))), and warn if a covariate is not constant
# within a key. `col_data$unit` holds the (sanitized) pseudobulk_unit values, so
# the lookup is joined on pseudobulk_unit.
.dea_add_covariates <- function(col_data, md, design, covariate_key, pseudobulk_unit) {
  vars <- setdiff(all.vars(stats::as.formula(design)), "group")

  if (length(vars) == 0L) {
    return(col_data)
  }

  if (!"unit" %in% colnames(col_data)) {
    stop("col_data must contain a 'unit' column.", call. = FALSE)
  }
  if (!covariate_key %in% colnames(md)) {
    stop("covariate_key column '", covariate_key, "' not found in metadata.", call. = FALSE)
  }
  if (!pseudobulk_unit %in% colnames(md)) {
    stop("pseudobulk_unit column '", pseudobulk_unit, "' not found in metadata.", call. = FALSE)
  }
  miss <- setdiff(vars, colnames(md))
  if (length(miss) > 0L) {
    stop(
      "design references covariate(s) not found in metadata: ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  key <- as.character(md[[covariate_key]])
  for (v in vars) {
    val <- md[[v]]
    n_distinct <- tapply(
      val, key,
      function(z) length(unique(z[!is.na(z)]))
    )
    if (any(n_distinct > 1L, na.rm = TRUE)) {
      warning(
        "Covariate '", v, "' is not constant within '", covariate_key,
        "'; using the first value per ", covariate_key, ".",
        call. = FALSE
      )
    }
  }
  keep <- !is.na(key) & !duplicated(key)
  lut <- md[keep, , drop = FALSE]

  .sanitize_unit <- function(u) ifelse(grepl("^[0-9]", u), paste0("g", u), u)

  m <- match(
    as.character(col_data$unit),
    .sanitize_unit(as.character(lut[[pseudobulk_unit]]))
  )

  if (anyNA(m)) {
    missing_units <- unique(as.character(col_data$unit)[is.na(m)])
    stop(
      "Could not match pseudobulk unit(s) in col_data back to metadata column '",
      pseudobulk_unit,
      "': ",
      paste(missing_units, collapse = ", "),
      call. = FALSE
    )
  }
  for (v in vars) {
    covariate_values <- lut[[v]][m]
    if (anyNA(covariate_values)) {
      stop(
        "Covariate '", v, "' contains missing values after matching pseudobulk units. ",
        "DESeq2 cannot use missing values in the design matrix.",
        call. = FALSE
      )
    }
    # Preserve numeric covariates such as age, RIN, percent.mt, etc.
    # Convert character/logical covariates to factors for DESeq2 design use.
    if (is.character(covariate_values) || is.logical(covariate_values)) {
      covariate_values <- factor(covariate_values)
    }
    col_data[[v]] <- covariate_values
  }

  col_data
}

# Aggregate raw RNA counts to pseudobulk samples (replicate-unit x group) using
# Seurat::PseudobulkExpression, and run DESeq2.
# Returns a standardized marker data.frame, or a `.dea_failed()` sentinel
# carrying a reason when the comparison cannot be run (too few replicates).
.de_pseudobulk <- function(seurat_obj, group1, group2,
                           pseudobulk_unit = "patient_id",
                           min_replicates = 3L,
                           min_cells_per_sample = 10L,
                           test = "Wald",
                           design = "~ group",
                           reduced = "~ 1",
                           low_count_filter = c(10L, 3L),
                           covariate_key = pseudobulk_unit,
                           shrinkage = TRUE,
                           drop_units = character()) {
  md <- seurat_obj@meta.data
  if (!pseudobulk_unit %in% colnames(md)) {
    stop("pseudobulk_unit column '", pseudobulk_unit, "' not found in metadata.", call. = FALSE)
  }

  # Aggregate (sum) raw counts per (replicate-unit x comparison group).
  pb <- Seurat::PseudobulkExpression(
    seurat_obj,
    assays = "RNA",
    layer = "counts",
    group.by = c(pseudobulk_unit, ".dea_group"),
    method = "aggregate",
    return.seurat = TRUE,
    verbose = FALSE
  )

  counts <- as.matrix(SeuratObject::LayerData(pb, assay = "RNA", layer = "counts"))
  pb_md <- pb@meta.data

  # Recover the replicate unit and group for each pseudobulk sample. When
  # return.seurat = TRUE, Seurat keeps the splitting variables in the metadata.
  if (all(c(pseudobulk_unit, ".dea_group") %in% colnames(pb_md))) {
    col_data <- data.frame(
      sample = rownames(pb_md),
      unit = as.character(pb_md[[pseudobulk_unit]]),
      group = as.character(pb_md[[".dea_group"]]),
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "Could not recover pseudobulk sample metadata from PseudobulkExpression output.",
      call. = FALSE
    )
  }
  col_data <- col_data[match(colnames(counts), col_data$sample), , drop = FALSE]

  # number of cells contributing to each pseudobulk sample (unit x group), so
  # that pseudobulk samples built from too few cells can be dropped.
  # PseudobulkExpression sanitizes the replicate-unit values (e.g. it prepends a
  # "g" to ids that start with a digit, "01-032" -> "g01-032"), so the unit in
  # `col_data` (from the pseudobulk metadata) will not match the raw metadata.
  # Sanitize the same way before matching, and fall back to keeping all samples
  # if the mapping still cannot be recovered (so pseudobulk is not silently
  # dropped because of an unexpected name change).
  .sanitize_unit <- function(u) ifelse(grepl("^[0-9]", u), paste0("g", u), u)
  cell_counts <- as.data.frame(
    table(
      unit = .sanitize_unit(as.character(md[[pseudobulk_unit]])),
      group = as.character(seurat_obj$.dea_group)
    ),
    stringsAsFactors = FALSE
  )
  col_data$n_cells <- cell_counts$Freq[match(
    paste(col_data$unit, col_data$group, sep = "\r"),
    paste(cell_counts$unit, cell_counts$group, sep = "\r")
  )]
  if (all(is.na(col_data$n_cells))) {
    col_data$n_cells <- min_cells_per_sample # mapping failed: keep all samples
  } else {
    col_data$n_cells[is.na(col_data$n_cells)] <- 0L
  }

  keep_sample <- col_data$n_cells >= min_cells_per_sample
  counts <- counts[, keep_sample, drop = FALSE]
  col_data <- col_data[keep_sample, , drop = FALSE]

  # Optionally drop outlier replicates (by unit) BEFORE the replicate gate, so the
  # min-replicate / min-cell cutoffs are re-checked on the post-removal sample set.
  n_dropped <- 0L
  if (length(drop_units) > 0L) {
    drop <- col_data$unit %in% drop_units
    n_dropped <- sum(drop)
    counts <- counts[, !drop, drop = FALSE]
    col_data <- col_data[!drop, , drop = FALSE]
  }

  ref <- if (identical(group2, "rest")) "rest" else group2
  n_per_group <- table(col_data$group)
  needed <- c(group1, ref)
  if (length(n_per_group) < 2L ||
    any(is.na(n_per_group[needed])) ||
    any(n_per_group[needed] < min_replicates)) {
    return(.dea_failed(paste0(
      "Pseudobulk not computed: fewer than ", min_replicates,
      " biological replicates (", pseudobulk_unit, ") with >= ",
      min_cells_per_sample, " cells in at least one group",
      if (n_dropped > 0L) paste0(" after removing ", n_dropped, " outlier sample(s)") else "",
      "."
    )))
  }

  # Add any design covariates (one row per replicate) before fitting.
  col_data <- .dea_add_covariates(col_data, md, design, covariate_key, pseudobulk_unit)

  run_deseq2(counts, col_data,
    group1 = group1, group2 = ref,
    test = test, design = design, reduced = reduced,
    low_count_filter = low_count_filter, shrinkage = shrinkage
  )
}

# Pseudobulk-sample PCA + outlier detection for one comparison's CLUSTER. Unlike
# .de_pseudobulk (which keeps only group1/group2), this aggregates EVERY patient
# with cells in the cluster (all groups) so the PCA can show the other samples
# too; samples are tagged group1 / group2 / "other". Runs DESeq2 VST
# (getVarianceStabilizedData) on the cluster pseudobulk, PCA on the top-variance
# genes, and mt::pca.outlier on the PC scores. Computed independently of the DE
# replicate gate (it is a QC view). Returns a list with status/coords/var_explained
# /outlier_samples, or status "not_computed" + reason when it cannot be built.
.de_pseudobulk_pca <- function(seurat_obj, group_by, cluster_by, cluster,
                               group1, group2, pseudobulk_unit = "patient_id",
                               min_cells_per_sample = 10L, n_top_genes = 500L,
                               outlier_conf = 0.975) {
  md <- seurat_obj@meta.data
  if (!pseudobulk_unit %in% colnames(md)) {
    return(list(
      status = "not_computed",
      reason = paste0("pseudobulk_unit '", pseudobulk_unit, "' not in metadata.")
    ))
  }
  # label cluster cells with their FULL group (NA elsewhere -> dropped on aggregate)
  keep_cluster <- if (length(cluster_by) == 1L && nzchar(cluster_by)) {
    as.character(md[[cluster_by]]) == cluster & !is.na(md[[cluster_by]])
  } else {
    rep(TRUE, nrow(md))
  }
  grp <- as.character(md[[group_by]])
  pca_group <- ifelse(keep_cluster & !is.na(grp), grp, NA_character_)
  seurat_obj <- SeuratObject::AddMetaData(
    seurat_obj, stats::setNames(pca_group, rownames(md)),
    col.name = ".pca_group"
  )

  pb <- Seurat::PseudobulkExpression(
    seurat_obj,
    assays = "RNA", layer = "counts",
    group.by = c(pseudobulk_unit, ".pca_group"),
    method = "aggregate", return.seurat = TRUE, verbose = FALSE
  )
  counts <- as.matrix(SeuratObject::LayerData(pb, assay = "RNA", layer = "counts"))
  pb_md <- pb@meta.data
  if (!all(c(pseudobulk_unit, ".pca_group") %in% colnames(pb_md))) {
    return(list(
      status = "not_computed",
      reason = "Could not recover pseudobulk sample metadata for PCA."
    ))
  }
  col_data <- data.frame(
    sample = rownames(pb_md), unit = as.character(pb_md[[pseudobulk_unit]]),
    group = as.character(pb_md[[".pca_group"]]), stringsAsFactors = FALSE
  )
  col_data <- col_data[match(colnames(counts), col_data$sample), , drop = FALSE]

  # drop pseudobulk samples built from too few cells (sanitize unit like Seurat does)
  .san <- function(u) ifelse(grepl("^[0-9]", u), paste0("g", u), u)
  cc <- as.data.frame(
    table(
      unit = .san(as.character(md[[pseudobulk_unit]])),
      group = as.character(seurat_obj$.pca_group)
    ),
    stringsAsFactors = FALSE
  )
  col_data$n_cells <- cc$Freq[match(
    paste(col_data$unit, col_data$group, sep = "\r"),
    paste(cc$unit, cc$group, sep = "\r")
  )]
  if (all(is.na(col_data$n_cells))) {
    col_data$n_cells <- min_cells_per_sample
  } else {
    col_data$n_cells[is.na(col_data$n_cells)] <- 0L
  }
  keep <- col_data$n_cells >= min_cells_per_sample
  counts <- counts[, keep, drop = FALSE]
  col_data <- col_data[keep, , drop = FALSE]

  if (ncol(counts) < 3L) {
    return(list(
      status = "not_computed",
      reason = sprintf(
        "PCA needs >= 3 pseudobulk samples in cluster '%s' (got %d).",
        cluster, ncol(counts)
      )
    ))
  }

  # DESeq2 VST on the cluster pseudobulk -> top-variance genes -> PCA.
  vmat <- tryCatch(
    {
      dds <- DESeq2::DESeqDataSetFromMatrix(
        round(counts), data.frame(row.names = col_data$sample),
        design = ~1
      )
      dds <- dds[rowSums(DESeq2::counts(dds)) > 0L, ]
      dds <- DESeq2::estimateSizeFactors(dds)
      suppressMessages(DESeq2::getVarianceStabilizedData(
        DESeq2::estimateDispersions(dds, fitType = "parametric")
      ))
    },
    error = function(e) NULL
  )
  if (is.null(vmat) || nrow(vmat) < 2L) {
    return(list(
      status = "not_computed",
      reason = "VST failed (too few genes/samples) for PCA."
    ))
  }
  rv <- apply(vmat, 1L, stats::var)
  top <- utils::head(order(rv, decreasing = TRUE), min(n_top_genes, nrow(vmat)))
  X <- t(vmat[top, , drop = FALSE]) # samples x genes
  pca <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  pct <- pca$sdev^2 / sum(pca$sdev^2)

  out_names <- tryCatch(
    names(mt::pca.outlier(X, conf.level = outlier_conf, plot = FALSE)$outlier),
    error = function(e) character()
  )

  ref <- if (identical(group2, "rest")) "rest" else group2
  role <- ifelse(col_data$group == group1, group1,
    ifelse(!identical(group2, "rest") & col_data$group == ref, ref, "other")
  )
  coords <- data.frame(
    sample = col_data$sample, unit = col_data$unit, group = col_data$group,
    comparison_role = role, n_cells = col_data$n_cells,
    PC1 = pca$x[, 1L], PC2 = if (ncol(pca$x) >= 2L) pca$x[, 2L] else 0,
    outlier = col_data$sample %in% out_names, stringsAsFactors = FALSE
  )

  list(
    status = "computed", coords = coords,
    var_explained = c(PC1 = pct[1L], PC2 = if (length(pct) >= 2L) pct[2L] else NA_real_),
    n_samples = nrow(coords), n_top_genes = length(top),
    outlier_samples = out_names, outlier_conf = outlier_conf, reason = NA_character_
  )
}

# =============================================================================
# GSEA (fgsea + MSigDB collections via msigdbr)
# =============================================================================

# Friendly collection label -> msigdbr (collection, subcollection).
.MSIGDB_COLLECTION_MAP <- list(
  "Hallmark"       = list(collection = "H", subcollection = NULL),
  "GO:BP"          = list(collection = "C5", subcollection = "GO:BP"),
  "C7:ImmuneSigDB" = list(collection = "C7", subcollection = "IMMUNESIGDB")
)

.default_gsea_params <- function(gsea_params = list()) {
  defaults <- list(
    collections = names(.MSIGDB_COLLECTION_MAP), # Hallmark, GO:BP, C7 ImmuneSigDB
    species = "Homo sapiens",
    minSize = 10L,
    maxSize = 500L,
    eps = 0,
    nproc = 1L, # serial fgsea (avoid spawning BiocParallel worker processes)
    padj_cutoff = 0.05,
    plot_number = 20L
  )
  params <- utils::modifyList(defaults, gsea_params)
  # `pathways` may be supplied to skip the msigdbr query; keep it out of the
  # validated parameter list but pass it through untouched.
  unknown <- setdiff(params$collections, names(.MSIGDB_COLLECTION_MAP))
  if (length(unknown) > 0L) {
    stop(
      "Unknown GSEA collection(s): ", paste(unknown, collapse = ", "),
      ". Available: ", paste(names(.MSIGDB_COLLECTION_MAP), collapse = ", "),
      call. = FALSE
    )
  }
  params
}

# Ranking metric for GSEA: sign(log2FC) * -log10(raw p-value), sorted decreasing.
# Positive ranks => up in group1; negative => up in group2 / rest.
.gsea_ranks <- function(markers) {
  m <- as.data.frame(markers)
  logfc_col <- .logfc_column(m)
  keep <- !is.na(m$p_val) & !is.na(m[[logfc_col]]) & !is.na(m$Gene) & nzchar(m$Gene)
  m <- m[keep, , drop = FALSE]
  if (nrow(m) == 0L) {
    return(numeric())
  }
  stat <- sign(m[[logfc_col]]) * -log10(pmax(m$p_val, .Machine$double.xmin))
  names(stat) <- m$Gene
  stat <- stat[!duplicated(names(stat))]
  sort(stat, decreasing = TRUE)
}

.empty_gsea_result <- function(reason, collection, ranks = numeric(), params = list()) {
  list(
    status = "not_computed",
    reason = reason,
    collection = collection,
    res.table = data.frame(),
    ranks = ranks,
    params = params
  )
}

# Accessors ----------------------------------------------------------------

.dea_gsea_collections <- function(x, level = NULL) {
  level <- .dea_default_level(x, level)
  g <- x@gsea[[level]]
  if (is.null(g)) character() else names(g)
}

.dea_default_collection <- function(x, collection = NULL, level = NULL) {
  cols <- .dea_gsea_collections(x, level)
  if (length(cols) == 0L) {
    return(NA_character_)
  }
  if (is.null(collection)) {
    return(cols[[1L]])
  }
  if (!collection %in% cols) {
    stop(
      "GSEA collection '", collection, "' was not computed. Available: ",
      paste(cols, collapse = ", "),
      call. = FALSE
    )
  }
  collection
}

.dea_gsea_result <- function(x, level = NULL, collection = NULL) {
  level <- .dea_default_level(x, level)
  g <- x@gsea[[level]]
  if (is.null(g) || length(g) == 0L) {
    return(NULL)
  }
  collection <- .dea_default_collection(x, collection, level)
  if (is.na(collection)) {
    return(NULL)
  }
  g[[collection]]
}

# =============================================================================
# Quarto report generation
# =============================================================================

.dea_sanitize_id <- function(...) {
  id <- paste(c(...), collapse = "-")
  id <- gsub("[^A-Za-z0-9-]+", "-", id)
  id <- gsub("-+", "-", id)
  gsub("^-+|-+$", "", id)
}

# Normalize/validate dea_report_lines' `interactive_plots`: NA/NULL/empty -> none
# (all static ggplot2); otherwise a subset of c("pca","volcano","go","gsea").
.dea_interactive_set <- function(interactive_plots) {
  if (length(interactive_plots) == 0L || all(is.na(interactive_plots))) {
    return(character())
  }
  ip <- as.character(interactive_plots)
  valid <- c("composition", "pca", "volcano", "go", "gsea")
  bad <- setdiff(ip, valid)
  if (length(bad) > 0L) {
    stop("interactive_plots must be NA or a subset of c(",
      paste(sprintf("'%s'", valid), collapse = ","), "); got: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  unique(ip)
}

.dea_level_label <- function(level) {
  switch(level,
    single_cell = "Single-cell analysis (FindMarkers)",
    pseudobulk = "Pseudobulk analysis (DESeq2)",
    level
  )
}

# A row of '#' of length n (markdown heading prefix).
.dea_h <- function(n) paste(rep("#", n), collapse = "")

# Quarto fig-height (inches) for a horizontal barplot with `n` bars: enough room
# per term, clamped to a sensible range.
.dea_fig_height <- function(n, per_term = 0.33, base = 1.6, min_h = 2.5, max_h = 16) {
  n <- suppressWarnings(as.numeric(n))
  if (length(n) == 0L || is.na(n) || n <= 0) {
    return(min_h)
  }
  round(max(min_h, min(max_h, base + per_term * n)), 1)
}

# Number of GO terms a go_barplot will display (the plot_data is already
# truncated to plot_number in .compute_topgo_result).
.dea_go_n_terms <- function(x, level, direction, ontology) {
  res <- .dea_go_result(x, level = level, direction = direction, ontology = ontology)
  if (is.null(res) || is.null(res$plot_data)) 0L else nrow(res$plot_data)
}

# Number of pathways a gsea_barplot will display (mirrors its selection: the
# significant sets, or all if none pass, capped at plot_number).
.dea_gsea_n_terms <- function(x, level, collection) {
  res <- .dea_gsea_result(x, level = level, collection = collection)
  if (is.null(res) || is.null(res$res.table) || nrow(res$res.table) == 0L) {
    return(0L)
  }
  tab <- res$res.table
  pc <- res$params$padj_cutoff %||% 0.05
  tn <- res$params$plot_number %||% 20L
  sig <- tab[!is.na(tab$padj) & tab$padj < pc, , drop = FALSE]
  if (nrow(sig) == 0L) sig <- tab[!is.na(tab$pval), , drop = FALSE]
  min(tn, nrow(sig))
}

# GSEA-yield grid (cell state/cluster x comparison) for one collection: the
# number of significantly enriched gene sets (padj < cutoff) per comparison.
# value is NA when the collection was NOT computed for that comparison/level (so
# white in the heatmap), and 0 when it was computed but no gene set passed the
# cutoff. `padj_cutoff` NULL uses each comparison's stored GSEA cutoff (else 0.05).
.dea_gsea_yield_data <- function(x, collection, level = "pseudobulk",
                                 padj_cutoff = NULL) {
  x_list <- normalize_dea_list(x)
  if (length(x_list) == 0L) {
    return(data.frame())
  }
  rows <- lapply(x_list, function(obj) {
    g <- obj@gsea[[level]]
    val <- NA_real_
    if (!is.null(g) && collection %in% names(g)) {
      res <- g[[collection]]
      rt <- res$res.table
      if (!is.null(rt) && nrow(rt) > 0L && "padj" %in% names(rt)) {
        cutoff <- padj_cutoff %||% res$params$padj_cutoff %||% 0.05
        val <- sum(!is.na(rt$padj) & rt$padj < cutoff)
      } else {
        val <- 0 # collection computed but produced no testable gene set
      }
    }
    data.frame(
      cluster = obj@cluster,
      comparison = paste(obj@group1, "vs", obj@group2),
      value = as.numeric(val), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Collections present across a list of comparisons at `level` (union, in order of
# first appearance) -> one GSEA-yield heatmap per collection.
.dea_gsea_yield_collections <- function(x, level = "pseudobulk") {
  x_list <- normalize_dea_list(x)
  cols <- character()
  for (obj in x_list) {
    g <- obj@gsea[[level]]
    if (!is.null(g)) cols <- c(cols, names(g))
  }
  unique(cols)
}

# Per-group and per-replicate cell counts for the comparison, computed from a
# metadata data.frame at render time (so no pipeline rerun is needed). Returns a
# list(per_group, per_unit, unit) or NULL if `meta` is unusable.
.dea_compare_counts <- function(x, meta) {
  if (is.null(meta)) {
    return(NULL)
  }
  m <- as.data.frame(meta)
  gb <- x@group_by
  if (length(gb) != 1L || !gb %in% colnames(m)) {
    return(NULL)
  }
  cb <- if (length(x@cluster_by) == 1L && !is.na(x@cluster_by) && nzchar(x@cluster_by)) {
    x@cluster_by
  } else {
    NA_character_
  }
  unit <- x@de_params$pseudobulk_unit %||% "patient_id"

  if (!is.na(cb) && cb %in% colnames(m)) {
    m <- m[as.character(m[[cb]]) == as.character(x@cluster) & !is.na(m[[cb]]), , drop = FALSE]
  }
  grp <- as.character(m[[gb]])
  g1 <- x@group1
  g2 <- if (length(x@group2) == 1L && nzchar(x@group2)) x@group2 else "rest"
  if (identical(g2, "rest")) {
    keep <- !is.na(grp)
    comp <- ifelse(grp == g1, g1, "rest")
  } else {
    keep <- grp %in% c(g1, g2)
    comp <- grp
  }
  m <- m[keep, , drop = FALSE]
  comp <- factor(comp[keep], levels = c(g1, g2))

  per_group <- data.frame(
    group = levels(comp),
    n_cells = as.integer(table(comp)),
    stringsAsFactors = FALSE
  )

  per_unit <- NULL
  if (unit %in% colnames(m)) {
    units <- as.character(m[[unit]])
    per_group$n_patients <- vapply(
      levels(comp), function(g) length(unique(units[as.character(comp) == g])), integer(1)
    )
    tab <- table(unit = units, group = comp)
    per_unit <- as.data.frame.matrix(tab[rowSums(tab) > 0, , drop = FALSE])
    per_unit <- cbind(stats::setNames(data.frame(rownames(per_unit)), unit), per_unit)
    rownames(per_unit) <- NULL
  }
  list(per_group = per_group, per_unit = per_unit, unit = unit)
}

# Markdown lines for the "Comparison summary": which groups/cluster (with the
# metadata column names), status/reason, and cell counts per group and per
# replicate. `meta` is an optional metadata data.frame for the per-replicate
# counts; without it, the stored per-group counts are used.
# Per-level summary lines for the "Comparison summary" tab: the adjusted-p cutoff
# (both levels), the DE-gene counts at that cutoff (both levels), and -- for
# pseudobulk only -- whether LFC shrinkage was used and by which method.
.dea_summary_level_lines <- function(x, level) {
  if (is.null(level) || !level %in% x@levels) {
    return(character())
  }
  out <- character()

  cutoff <- x@padj_cutoffs[[level]] %||% 0.05
  out <- c(out, sprintf("- **p-value cutoff (adjusted):** %s", format(cutoff)))

  de <- x@de[[level]]
  if (!is.null(de) && nrow(de) > 0L && "p_val_adj" %in% names(de)) {
    sig <- !is.na(de$p_val_adj) & de$p_val_adj <= cutoff
    n_up <- sum(sig & de$avg_log2FC > 0, na.rm = TRUE)
    n_dn <- sum(sig & de$avg_log2FC < 0, na.rm = TRUE)
    out <- c(out, sprintf(
      "- **Differentially Expressed Genes:** %d up-regulated and %d down-regulated (%d total)",
      n_up, n_dn, n_up + n_dn
    ))
  }

  if (identical(level, "pseudobulk")) {
    requested <- isTRUE(x@de_params$pb_lfc_shrink)
    method <- x@de_params$pb_lfc_shrink_method %||%
      attr(x@de[["pseudobulk"]], "lfcShrink_method")
    shr <- if (requested && !is.null(method) &&
      !method %in% c("none", "disabled")) {
      sprintf("Yes (%s)", method)
    } else if (requested) {
      "Yes"
    } else {
      "No"
    }
    out <- c(out, sprintf("- **Shrinkage:** %s", shr))

    # Whether PCA-flagged outliers were used or removed for the DESeq2 comparison.
    removed <- isTRUE(x@de_params$pb_remove_outliers)
    out <- c(out, sprintf(
      "- **Outliers:** %s for the DESeq2 comparison.",
      if (removed) "removed" else "kept"
    ))
  }

  c(out, "")
}

.dea_summary_content <- function(x, meta = NULL, level = NULL) {
  cluster_txt <- if (length(x@cluster_by) == 1L && nzchar(x@cluster_by)) {
    sprintf("`%s` (metadata column `%s`)", x@cluster, x@cluster_by)
  } else {
    sprintf("`%s`", x@cluster)
  }
  lines <- c(
    sprintf(
      "- **Groups compared:** `%s` vs `%s`  (metadata column `%s`)",
      x@group1, x@group2, x@group_by
    ),
    sprintf("- **Cluster:** %s", cluster_txt),
    sprintf("- **Status:** %s", x@status),
    ""
  )
  if (identical(x@status, "computed")) {
    lines <- c(lines, sprintf(
      "- **Levels computed:** %s",
      paste(x@levels, collapse = ", ")
    ), "")
    lines <- c(lines, .dea_summary_level_lines(x, level))
  } else {
    lines <- c(
      lines,
      "::: {.callout-warning}",
      "#### Comparison not computed", "",
      x@reason %||% "Not computed.", "",
      ":::", ""
    )
  }

  cc <- .dea_compare_counts(x, meta)
  if (!is.null(cc)) {
    lines <- c(
      lines, "**Cells per group**", "",
      as.character(knitr::kable(cc$per_group, format = "pipe")), ""
    )
    if (!is.null(cc$per_unit)) {
      lines <- c(
        lines, sprintf("**Cells per %s and group**", cc$unit), "",
        as.character(knitr::kable(cc$per_unit, format = "pipe")), ""
      )
    }
  } else if (nrow(as.data.frame(x@n_cells)) > 0L) {
    lines <- c(
      lines, "**Cells per group**", "",
      as.character(knitr::kable(as.data.frame(x@n_cells), format = "pipe")), ""
    )
  }
  lines
}

# Build the Quarto/knitr child lines that render one analysis level
# (comparison summary + volcano + DE table + GO + GSEA) for a comparison. The
# generated R chunks reference the object via the variable name `obj_name`
# (which must be bound in the environment passed to knitr::knit_child).
#   base_level : heading level of the analysis tabset tabs.
#                GO plot-type tabs are base_level+1; ontology x direction tabs
#                are base_level+2.
.dea_level_lines <- function(x, level, fig_id, obj_name = "comparison",
                             base_level = 4L, meta = NULL,
                             interactive_plots = character(),
                             gsea_metric = "signed_nlog10_padj") {
  onto <- .dea_ontologies(x)
  collections <- .dea_gsea_collections(x, level = level)
  q <- function(s) shQuote(s)
  bool <- function(b) if (isTRUE(b)) "TRUE" else "FALSE"
  h_tab <- .dea_h(base_level) # analysis tabs
  h_kind <- .dea_h(base_level + 1L) # GO plot-type / GSEA collection
  h_dir <- .dea_h(base_level + 2L) # ontology x direction

  # Pseudobulk PCA: a sentence listing the outlier samples BEFORE the tabset, and
  # a "PCA" tab (with the VST PCA plot) between Comparison summary and Volcano.
  pca <- x@pca[[level]]
  has_pca <- identical(level, "pseudobulk") &&
    !is.null(pca) && identical(pca$status, "computed") && !is.null(pca$coords)
  pre_lines <- character()
  pca_tab <- character()
  if (has_pca) {
    outl <- unique(pca$coords$unit[pca$coords$outlier])
    outl <- sub("^g", "", outl)
    pre_lines <- if (length(outl) > 0L) {
      c(sprintf(
        "**The following samples were identified as outliers:** %s.",
        paste(outl, collapse = ", ")
      ), "")
    } else {
      c("*No samples were identified as outliers.*", "")
    }
    pca_tab <- c(
      paste(h_tab, "PCA"), "",
      "```{r}", "#| column: page", "#| message: false", "#| warning: false",
      sprintf("#| label: fig-pca-%s", .dea_sanitize_id(fig_id, level)),
      sprintf(
        "pca_plot(%s, level = %s, interactive = %s)",
        obj_name, q(level), bool("pca" %in% interactive_plots)
      ),
      "```", ""
    )
  }

  lines <- c(
    pre_lines,
    "::: {.panel-tabset}",
    "",
    paste(h_tab, "Comparison summary"),
    "",
    .dea_summary_content(x, meta = meta, level = level),
    "",
    pca_tab,
    paste(h_tab, "Volcano plot"),
    "",
    "```{r}",
    "#| column: page",
    "#| message: false",
    "#| warning: false",
    sprintf("#| label: fig-volcano-%s", .dea_sanitize_id(fig_id, level)),
    sprintf(
      "volcano_plot(%s, level = %s, interactive = %s)",
      obj_name, q(level), bool("volcano" %in% interactive_plots)
    ),
    "```",
    "",
    paste(h_tab, "DE table"),
    "",
    "```{r}",
    "#| column: page",
    "#| message: false",
    "#| warning: false",
    sprintf("#| label: tbl-de-%s", .dea_sanitize_id(fig_id, level)),
    sprintf(
      "#| tbl-cap: \"Differential-expression results (%s level): %s vs %s in cluster %s.\"",
      level, x@group1, x@group2, x@cluster
    ),
    sprintf(
      "scitargets:::.scitargets_dt2_tbl(markers_table(%s, level = %s), page_length = 20)",
      obj_name, q(level)
    ),
    "```",
    ""
  )

  # GO: separate Barplot / Network / Genes-Terms tabs, each split by
  # ontology x direction. Barplot height scales with the number of terms.
  emit_go_onto_dir <- function(kind_lines_fun) {
    out <- c("::: {.panel-tabset}", "")
    for (o in onto) {
      for (dir in c("up", "down")) {
        out <- c(
          out, paste(h_dir, sprintf("%s %s", o, toupper(dir))), "",
          kind_lines_fun(o, dir), ""
        )
      }
    }
    c(out, ":::", "")
  }

  # GO tab: only emit when there are GO terms to show (skip the whole "Gene
  # Ontology" tab when GO could not be computed, e.g. no significant DE genes
  # -> no GO terms). P2.14.1.
  go_total <- sum(vapply(onto, function(o) {
    sum(vapply(
      c("up", "down"),
      function(d) .dea_go_n_terms(x, level, d, o), numeric(1L)
    ))
  }, numeric(1L)))
  if (go_total > 0L) {
    lines <- c(lines, paste(h_tab, "Gene Ontology"), "", "::: {.panel-tabset}", "")

    # Barplot
    lines <- c(lines, paste(h_kind, "Barplot"), "")
    lines <- c(lines, emit_go_onto_dir(function(o, dir) {
      tag <- .dea_sanitize_id(fig_id, level, "go-bar", o, dir)
      c(
        "```{r}", "#| column: page", "#| message: false", "#| warning: false",
        sprintf("#| fig-height: %s", .dea_fig_height(.dea_go_n_terms(x, level, dir, o))),
        sprintf("#| label: fig-go-%s", tag),
        sprintf(
          "go_barplot(%s, direction = %s, ontology = %s, level = %s, interactive = %s)",
          obj_name, q(dir), q(o), q(level), bool("go" %in% interactive_plots)
        ),
        "```"
      )
    }))

    # Network plot
    lines <- c(lines, paste(h_kind, "Network plot"), "")
    lines <- c(lines, emit_go_onto_dir(function(o, dir) {
      tag <- .dea_sanitize_id(fig_id, level, "go-net", o, dir)
      c(
        "```{r}", "#| column: page", "#| message: false", "#| warning: false",
        sprintf("#| label: fig-cnet-%s", tag),
        sprintf(
          "go_cnetplot(%s, direction = %s, ontology = %s, level = %s, interactive = %s)",
          obj_name, q(dir), q(o), q(level), bool("go" %in% interactive_plots)
        ),
        "```"
      )
    }))

    # Genes-Terms list
    lines <- c(lines, paste(h_kind, "Genes-Terms list"), "")
    lines <- c(lines, emit_go_onto_dir(function(o, dir) {
      c(
        "```{r}", "#| output: asis", "#| message: false", "#| warning: false",
        sprintf(
          "cat(go_genes_html(%s, direction = %s, ontology = %s, level = %s))",
          obj_name, q(dir), q(o), q(level)
        ),
        "```"
      )
    }))

    lines <- c(lines, ":::", "") # close GO (Barplot/Network/Genes) tabset
  } # end if(go_total > 0L): omit the entire "Gene Ontology" tab when no GO terms

  # GSEA: one tab per collection; barplot height scales with the number shown.
  lines <- c(lines, paste(h_tab, "GSEA"), "")
  if (length(collections) == 0L) {
    lines <- c(lines, "_No GSEA results for this level._", "")
  } else {
    lines <- c(lines, "::: {.panel-tabset}", "")
    for (col in collections) {
      tag <- .dea_sanitize_id(fig_id, level, col)
      lines <- c(
        lines,
        paste(h_kind, col),
        "",
        "```{r}",
        "#| column: page",
        "#| message: false",
        "#| warning: false",
        sprintf("#| fig-height: %s", .dea_fig_height(.dea_gsea_n_terms(x, level, col))),
        sprintf("#| label: fig-gsea-%s", tag),
        sprintf(
          "gsea_barplot(%s, collection = %s, level = %s, metric = %s, interactive = %s)",
          obj_name, q(col), q(level), q(gsea_metric),
          bool("gsea" %in% interactive_plots)
        ),
        "```",
        "",
        "```{r}",
        "#| column: page",
        "#| message: false",
        "#| warning: false",
        sprintf("#| label: tbl-gsea-%s", tag),
        sprintf(
          "#| tbl-cap: \"GSEA results: %s collection (%s level), %s vs %s in cluster %s.\"",
          col, level, x@group1, x@group2, x@cluster
        ),
        sprintf(
          "scitargets:::.scitargets_dt2_tbl(gsea_table(%s, collection = %s, level = %s), page_length = 15)",
          obj_name, q(col), q(level)
        ),
        "```",
        ""
      )
    }
    lines <- c(lines, ":::", "")
  }

  lines <- c(lines, ":::", "") # close the analysis (volcano/DE/GO/GSEA) tabset
  lines
}
