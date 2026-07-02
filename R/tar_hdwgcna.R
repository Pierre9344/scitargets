# Local variables that only ever appear inside the bquote()-generated target
# commands (so R CMD check does not flag them as undefined globals).
utils::globalVariables(c(
  "obj", "power_table", "pass", "tg", "cc", "md", "cells", "comparisons",
  "dme_list", "g1", "g2", "b1", "b2", "d"
))

# Lower-cased, syntactically-valid target-name suffix for a group label
# (e.g. "MAIT" -> "mait", "CD8 TEM" -> "cd8_tem").
.hdwgcna_suffix <- function(group) {
  s <- tolower(base::gsub("[^A-Za-z0-9]+", "_", group))
  base::gsub("^_+|_+$", "", s)
}

#' Targets factory for an hdWGCNA co-expression analysis
#'
#' Builds the `targets` pipeline steps for a high-dimensional WGCNA (hdWGCNA)
#' co-expression network analysis, following the hdWGCNA basic + differential
#' module-eigengene tutorials. A single shared metacell-preparation step
#' (`wgcna_prep`) is created once, and the per-group steps below are created for
#' every element of `group` (target names carry the group, e.g. `wgcna_mait`):
#'
#' - `wgcna_prep` — RNA assay prep (join layers, normalise, variable features,
#'   scale) then `SetupForWGCNA` + `MetacellsByGroups` and the metacell embedding.
#' - `wgcna_<group>_powertest` — `SetDatExpr` for the group + `TestSoftPowers`.
#' - `wgcna_<group>_soft_power` — the soft power chosen exactly as
#'   `hdWGCNA::PlotSoftPowers` highlights it (lowest tested power with
#'   `SFT.R.sq >= sft_rsquared`).
#' - `wgcna_<group>` — `ConstructNetwork` (at the selected soft power) +
#'   `ModuleEigengenes` + `ModuleConnectivity` + `ResetModuleNames`, and, when
#'   `trait_col` is given, `ModuleTraitCorrelation` of the modules against one-hot
#'   indicators of each `trait_groups` level. The correlation is computed once
#'   across all of the group's cells (a single constant `group.by`), so the stored
#'   result has only an `"all_cells"` entry and is not split by the object's
#'   cluster identities.
#' - `wgcna_<group>_dmes` — differential module eigengenes (`FindDMEs`) for every
#'   pairwise comparison of `trait_groups` within the group (only created when
#'   `trait_col` and `trait_groups` are supplied).
#' - `wgcna_<group>_enrichr` — module gene-set enrichment via `RunEnrichr`, run
#'   once per `enrichr_dbs` entry and row-bound into a single data frame (each row
#'   tagged with its `db`); created when `run_enrichr = TRUE`.
#'
#' All values passed to the hdWGCNA / WGCNA functions are exposed as arguments.
#' Every step that uses WGCNA multithreading or loads the whole shared Seurat
#' object is **pinned to `main`** and is *not* affected by `deployment` /
#' `controller`: `wgcna_prep`, the power test, the soft-power pick, the network
#' and the Enrichr step. WGCNA's threading is incompatible with `crew` workers,
#' and the power test / network hold a full copy of the object. Only the
#' differential-module-eigengene step (`wgcna_<group>_dmes`, a Wilcoxon test)
#' honours `deployment` and `controller`, so it can be sent to a worker (and to a
#' named controller when a controller group is used). `tom_name` is set to the
#' group label.
#'
#' @param group Character vector of cell-group labels to analyse (the values found
#'   in `clustering_col`, e.g. `"MAIT"`). Per-group target names use a lower-cased,
#'   sanitised version of the label.
#' @param create_prep Whether this call creates the shared `wgcna_prep` target.
#'   `TRUE` (default) builds it; set `FALSE` on the additional `tar_hdwgcna()`
#'   calls that reuse the same `wgcna_prep` (exactly one call must use `TRUE`).
#'   When `FALSE`, `input_obj` is not required.
#' @param input_obj Name of the upstream Seurat-object target to start from
#'   (required when `create_prep = TRUE`).
#' @param wgcna_name Name of the hdWGCNA experiment (`SetupForWGCNA`).
#' @param assay Assay used throughout (default `"RNA"`, log-normalised).
#' @param clustering_col Metadata column holding the cell groups; used as
#'   `ident.group` for the metacells and as `group.by`/`subset_by` downstream.
#' @param patient_col Biological-replicate column; used as the metacell harmony
#'   variable and as `group.by.vars` for `ModuleEigengenes`.
#' @param metacell_group_by `group.by` for `MetacellsByGroups`. Defaults to
#'   `c(clustering_col, patient_col)`.
#' @param reduction Reduction used for the metacell KNN and metacell UMAP.
#' @param gene_select,fraction `SetupForWGCNA` gene selection method and fraction.
#' @param metacell_k,metacell_max_shared `MetacellsByGroups` `k` and `max_shared`.
#' @param metacell_dims Dims used for `RunUMAPMetacells`.
#' @param network_type `networkType` for `TestSoftPowers` / `ConstructNetwork`.
#' @param sft_rsquared Scale-free-topology R^2 threshold for the soft-power pick.
#' @param tom_outdir Directory where `ConstructNetwork` writes the TOM.
#' @param overwrite_tom Whether `ConstructNetwork` overwrites an existing TOM.
#' @param n_threads Threads for `WGCNA::enableWGCNAThreads`.
#' @param deployment Where the *relocatable* steps run: `"main"` (default, the
#'   whole pipeline on the main process) or `"worker"` (send them to a `crew`
#'   worker). Only the `dmes` step is relocatable; the multithreaded /
#'   whole-object steps are always on `"main"` (see Details).
#' @param controller Optional name of a `crew` controller to route the
#'   relocatable steps to (via `targets::tar_resources()`), for use with a
#'   `crew::crew_controller_group()`. `NULL` (default) uses the pipeline's default
#'   controller when `deployment = "worker"`, and is ignored on `"main"`.
#' @param trait_col Metadata column with the (categorical) trait to associate the
#'   modules with, e.g. a clinical group. When `NULL`, the module-trait
#'   correlation and the DME steps are skipped.
#' @param trait_groups Character vector of `trait_col` levels to use. One-hot
#'   encoded for `ModuleTraitCorrelation` and compared pairwise (all combinations,
#'   first element = test group) for `FindDMEs`.
#' @param mtc_cor_method Correlation method for `ModuleTraitCorrelation`
#'   (`"spearman"` or `"pearson"`). Default `"spearman"`: rank-based, so it is
#'   robust to outliers and to the non-normal / non-linear module-eigengene
#'   distributions that Pearson assumes away.
#' @param dme_test Test passed to `FindDMEs` (`test.use`).
#' @param dme_harmonized Whether `FindDMEs` uses the harmonized module eigengenes.
#' @param run_enrichr Whether to create the `wgcna_<group>_enrichr` enrichment
#'   target (default `TRUE`). Set `FALSE` to skip it (e.g. when offline, since
#'   `RunEnrichr` queries the Enrichr web API).
#' @param enrichr_dbs Character vector of Enrichr database names to test, each run
#'   in its own `RunEnrichr` call. Defaults to the 2023 GO BP/CC/MF databases.
#' @param enrichr_max_genes `max_genes` passed to `RunEnrichr` (use `Inf` for all
#'   genes in a module).
#' @param enrichr_modules Which modules to send to the GO/Enrichr analysis. One of:
#'   \describe{
#'     \item{`"all"`}{(default) every non-grey module is tested.}
#'     \item{`"correlated"`}{only modules correlated with at least one
#'       `trait_groups` level in the stored `ModuleTraitCorrelation`
#'       (FDR < `enrichr_corr_cutoff`), via [hdwgcna_associated_modules()].}
#'     \item{`"dmes"`}{only modules with at least one significant differential
#'       module eigengene in `wgcna_<group>_dmes` (`p_val_adj < enrichr_dme_cutoff`
#'       for any pairwise comparison).}
#'     \item{`"both"`}{the union of the `"correlated"` and `"dmes"` module sets.}
#'   }
#'   For every option other than `"all"` the non-selected modules are skipped
#'   entirely, so no Enrichr API calls are made for them. Any option other than
#'   `"all"` requires `trait_col` to be set (the selection reads the module-trait
#'   correlation and/or the DME target); if no module passes the relevant cutoff
#'   the enrichment target is an empty data frame. This only affects the module
#'   selection for the GO analysis -- the `ModuleTraitCorrelation` and `FindDMEs`
#'   targets keep their full, unfiltered results.
#' @param enrichr_corr_cutoff FDR cutoff on the module-trait correlation used to
#'   pick the `"correlated"` modules (default 0.05). Only used when
#'   `enrichr_modules` is `"correlated"` or `"both"`.
#' @param enrichr_dme_cutoff Adjusted p-value cutoff (`p_val_adj`) on the
#'   differential module eigengenes used to pick the `"dmes"` modules
#'   (default 0.05). Only used when `enrichr_modules` is `"dmes"` or `"both"`.
#' @param module_table_file Either `NULL` (default, writes nothing) or a valid
#'   path to a single `.xlsx` document to write, describing the composition of
#'   every module (one worksheet per cell group, sheet named after the group).
#'   The filename is taken from this string and **must end in `.xlsx`** (an error
#'   is raised otherwise); any missing parent directories are created. Each sheet
#'   is the long-format
#'   [hdwgcna_module_composition()] table (one row per gene: `module`, `color`,
#'   `gene`, `kME`). `NULL` (default) writes nothing. When set, a light per-group
#'   `wgcna_<group>_modcomp` target is created plus a single
#'   `wgcna_module_composition` file target that assembles the workbook.
#' @param module_gene_list Whether to create the single nested module gene-list
#'   target `wgcna_module_gene_list` (default `TRUE`), a list with one element per
#'   cell group, each itself a named list of the genes in each module
#'   (`ls[[<cell type>]][[<module>]]`). It is assembled from the light per-group
#'   `wgcna_<group>_modgenes` targets (always created). The aggregate has a fixed
#'   name, so when several `tar_hdwgcna()` calls are spliced into one pipeline set
#'   `module_gene_list = FALSE` on all but one to avoid a duplicate-target error.
#'
#' @returns A list of `targets` objects to splice into a `targets` pipeline. With
#'   `<g>` the lower-cased, sanitised label of each `group` element (e.g. `"MAIT"`
#'   becomes `mait`, `"CD8 TEM"` becomes `cd8_tem`), the following steps are created:
#'   - `wgcna_prep`: the shared metacell preparation (`deployment = "main"`).
#'     Created only when `create_prep = TRUE` (the default); omit it (set
#'     `create_prep = FALSE`) on the other `tar_hdwgcna()` calls that reuse the
#'     same `wgcna_prep`.
#'   - `wgcna_<g>_powertest`: `SetDatExpr` for the group + `TestSoftPowers`
#'     (one per group; `deployment = "main"`, since it loads the whole Seurat
#'     object).
#'   - `wgcna_<g>_soft_power`: the selected soft power, a single number
#'     (one per group; `deployment = "main"`).
#'   - `wgcna_<g>`: `ConstructNetwork` + `ModuleEigengenes` + `ModuleConnectivity`
#'     + `ResetModuleNames`, plus `ModuleTraitCorrelation` when `trait_col` is set
#'     (one per group; `deployment = "main"`). This is the main per-group object.
#'   - `wgcna_<g>_dmes`: differential module eigengenes (`FindDMEs`) for every
#'     pairwise `trait_groups` comparison; created only when `trait_col` is
#'     supplied (one per group; runs on a worker).
#'   - `wgcna_<g>_enrichr`: a data frame of `RunEnrichr` module enrichment results,
#'     one `RunEnrichr` call per `enrichr_dbs` entry, row-bound (each row tagged
#'     with its `db`); created only when `run_enrichr = TRUE` (one per group;
#'     `deployment = "main"`).
#'   - `wgcna_<g>_modcomp` and `wgcna_module_composition`: created only when
#'     `module_table_file` is set. The first is the per-group long-format module
#'     composition data frame; the second (one per call) is a `format = "file"`
#'     target that writes the `.xlsx` workbook named by `module_table_file` (one
#'     sheet per group) and returns its path.
#'   - `wgcna_<g>_modgenes`: the group's module-to-gene list -- a named `list` with
#'     one element per (non-grey) module holding the character vector of that
#'     module's genes (ordered by decreasing kME), via [hdwgcna_module_gene_list()].
#'     Always created (one per group; `deployment = "main"`).
#'   - `wgcna_module_gene_list`: a single nested `list` assembling the per-group
#'     `wgcna_<g>_modgenes` into `ls[[<cell type>]][[<module>]] = genes` across all
#'     of this call's groups. Created when `module_gene_list = TRUE` (the default);
#'     one per call, with a fixed name, so set `module_gene_list = FALSE` on all but
#'     one of several `tar_hdwgcna()` calls (`deployment = "main"`).
#'   - `wgcna_<g>_modscore`: a per-cell module-score `data.frame` with one row per
#'     cell **of the whole object (all cell types)** -- columns `cell`, `cell_type`
#'     (each cell's own `clustering_col` label), `condition` (its `trait_col` value,
#'     or `NA` when `trait_col` is unset) and one numeric column per module --
#'     computed with [hdwgcna_module_scores()] / [Seurat::AddModuleScore()] using the
#'     group's modules. Always created (one per group; `deployment = "main"`).
#' @export
#'
#' @examples
#' \dontrun{
#' list(
#'   scitargets::tar_hdwgcna(
#'     group = "MAIT",
#'     input_obj = "seurat_merged_filtered",
#'     clustering_col = "predicted.celltype.l2",
#'     patient_col = "patient_id",
#'     trait_col = "cohort_clinic_preterrah",
#'     trait_groups = c("CR1M/6M", "IR6M", "NR1M")
#'   )
#' )
#' }
tar_hdwgcna <- function(
  group,
  create_prep = TRUE,
  input_obj = NULL,
  wgcna_name = "hdwgcna",
  assay = "RNA",
  clustering_col = NULL,
  patient_col = NULL,
  metacell_group_by = NULL,
  reduction = "harmony",
  gene_select = "fraction",
  fraction = 0.05,
  metacell_k = 25,
  metacell_max_shared = 10,
  metacell_dims = 1:15,
  network_type = "signed",
  sft_rsquared = 0.8,
  tom_outdir = "./out/hdwgcna/",
  overwrite_tom = TRUE,
  n_threads = 8,
  deployment = "main",
  controller = NULL,
  trait_col = NULL,
  trait_groups = NULL,
  mtc_cor_method = "spearman",
  dme_test = "wilcox",
  dme_harmonized = TRUE,
  run_enrichr = TRUE,
  enrichr_dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023"),
  enrichr_max_genes = 100,
  enrichr_modules = c("all", "correlated", "dmes", "both"),
  enrichr_corr_cutoff = 0.05,
  enrichr_dme_cutoff = 0.05,
  module_table_file = NULL,
  module_gene_list = TRUE
) {
  if (!is.character(group) || length(group) == 0L) {
    stop("`group` must be a non-empty character vector of cell-group labels.")
  }
  if (!is.logical(create_prep) || length(create_prep) != 1L || is.na(create_prep)) {
    stop("`create_prep` must be a single TRUE/FALSE.")
  }
  if (!is.logical(run_enrichr) || length(run_enrichr) != 1L || is.na(run_enrichr)) {
    stop("`run_enrichr` must be a single TRUE/FALSE.")
  }
  if (!is.logical(module_gene_list) || length(module_gene_list) != 1L || is.na(module_gene_list)) {
    stop("`module_gene_list` must be a single TRUE/FALSE.")
  }
  if (isTRUE(run_enrichr) && (!is.character(enrichr_dbs) || length(enrichr_dbs) == 0L)) {
    stop("`enrichr_dbs` must be a non-empty character vector of Enrichr database names.")
  }
  # `enrichr_modules` and the two cutoffs only govern which modules are sent to
  # the GO/Enrichr analysis; they do not change the ModuleTraitCorrelation or the
  # FindDMEs targets (those keep their full, unfiltered results).
  enrichr_modules <- match.arg(enrichr_modules, c("all", "correlated", "dmes", "both"))
  enrichr_restrict <- !identical(enrichr_modules, "all")
  enrichr_need_corr <- enrichr_modules %in% c("correlated", "both")
  enrichr_need_dme <- enrichr_modules %in% c("dmes", "both")
  if (isTRUE(run_enrichr) && enrichr_restrict && is.null(trait_col)) {
    stop("`enrichr_modules = \"", enrichr_modules, "\"` selects modules from the ",
      "module-trait correlation and/or the differential module eigengenes, so ",
      "`trait_col` (and `trait_groups`) must be set.")
  }
  validate_cutoff <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x > 1) {
      stop("`", nm, "` must be a single number in (0, 1].")
    }
  }
  validate_cutoff(enrichr_corr_cutoff, "enrichr_corr_cutoff")
  validate_cutoff(enrichr_dme_cutoff, "enrichr_dme_cutoff")
  if (!is.null(module_table_file)) {
    if (!is.character(module_table_file) || length(module_table_file) != 1L || !nzchar(module_table_file)) {
      stop("`module_table_file` must be NULL or a valid path to a .xlsx document (a single non-empty string).")
    }
    if (!grepl("\\.xlsx$", module_table_file, ignore.case = TRUE)) {
      stop("`module_table_file` must be NULL or a valid path to a .xlsx document; '",
        module_table_file, "' does not end in '.xlsx'.")
    }
  }
  # input_obj only feeds wgcna_prep, so it is required only when this call builds it.
  if (isTRUE(create_prep) &&
    (is.null(input_obj) || !is.character(input_obj) || length(input_obj) != 1L || !nzchar(input_obj))) {
    stop("`input_obj` must be set to the name of the upstream Seurat-object target (a single non-empty string) when `create_prep = TRUE`.")
  }
  if (is.null(clustering_col) || !is.character(clustering_col) || length(clustering_col) != 1L || !nzchar(clustering_col)) {
    stop("`clustering_col` must be set to the metadata column holding the cell groups / clusters (a single non-empty string).")
  }
  if (is.null(patient_col) || !is.character(patient_col) || length(patient_col) != 1L || !nzchar(patient_col)) {
    stop("`patient_col` must be set to the biological-replicate metadata column (a single non-empty string).")
  }
  if (!is.null(trait_col) && (!is.character(trait_col) || length(trait_col) != 1L)) {
    stop("`trait_col` must be NULL or a single metadata column name.")
  }
  if (!is.null(trait_col) && (is.null(trait_groups) || length(trait_groups) < 2L)) {
    stop("When `trait_col` is set, `trait_groups` must list at least two levels.")
  }
  if (is.null(metacell_group_by)) {
    metacell_group_by <- c(clustering_col, patient_col)
  }

  # Deployment of the *relocatable* steps. Every step that uses WGCNA
  # multithreading (or loads the whole shared Seurat object) stays pinned to
  # "main" regardless of these arguments -- WGCNA's threading is incompatible
  # with crew workers. Only the differential-module-eigengene step (`dmes`,
  # a Wilcoxon test) honours `deployment`/`controller`.
  deployment <- match.arg(deployment, c("main", "worker"))
  movable_resources <- if (is.null(controller)) {
    targets::tar_option_get("resources")
  } else {
    targets::tar_resources(crew = targets::tar_resources_crew(controller = controller))
  }

  # ---- shared metacell preparation (created once across all tar_hdwgcna calls) ----
  prep <- if (isTRUE(create_prep)) list(targets::tar_target_raw(
    name = "wgcna_prep",
    command = bquote({
      library(WGCNA)
      library(hdWGCNA)
      WGCNA::enableWGCNAThreads(nThreads = .(n_threads))
      obj <- .(as.name(input_obj))
      SeuratObject::DefaultAssay(obj) <- .(assay)
      obj[[.(assay)]] <- SeuratObject::JoinLayers(obj[[.(assay)]])
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
      obj <- Seurat::ScaleData(obj, features = SeuratObject::VariableFeatures(obj), verbose = FALSE)
      obj <- hdWGCNA::SetupForWGCNA(obj, gene_select = .(gene_select), fraction = .(fraction), wgcna_name = .(wgcna_name))
      obj <- hdWGCNA::MetacellsByGroups(obj, group.by = .(metacell_group_by), reduction = .(reduction), k = .(metacell_k), max_shared = .(metacell_max_shared), ident.group = .(clustering_col))
      obj <- hdWGCNA::NormalizeMetacells(obj)
      obj <- hdWGCNA::ScaleMetacells(obj, features = SeuratObject::VariableFeatures(obj))
      obj <- hdWGCNA::RunPCAMetacells(obj, features = SeuratObject::VariableFeatures(obj))
      obj <- hdWGCNA::RunHarmonyMetacells(obj, group.by.vars = .(patient_col))
      obj <- hdWGCNA::RunUMAPMetacells(obj, reduction = .(reduction), dims = .(metacell_dims))
      obj
    }, where = environment()),
    deployment = "main",
    description = "hdWGCNA: shared metacell preparation (RNA: join layers, normalize, scale, SetupForWGCNA, metacells)"
  )) else list()

  per_group <- lapply(group, function(g) {
    suffix <- .hdwgcna_suffix(g)
    pt_name <- base::paste0("wgcna_", suffix, "_powertest")
    sp_name <- base::paste0("wgcna_", suffix, "_soft_power")
    net_name <- base::paste0("wgcna_", suffix)
    dme_name <- base::paste0("wgcna_", suffix, "_dmes")
    enr_name <- base::paste0("wgcna_", suffix, "_enrichr")

    powertest <- targets::tar_target_raw(
      name = pt_name,
      command = bquote({
        library(WGCNA)
        library(hdWGCNA)
        # Runs on main (deployment = "main"): SetDatExpr loads the whole shared
        # Seurat object (wgcna_prep), so dispatching several groups to crew workers
        # at once would hold multiple full copies in memory. On main the groups run
        # one at a time. Kept single-threaded (no enableWGCNAThreads(), whose socket
        # cluster also clashes under concurrency); the multithreaded step is the
        # network target (ConstructNetwork), also on main.
        obj <- hdWGCNA::SetDatExpr(.(as.name("wgcna_prep")), group_name = .(g), group.by = .(clustering_col), assay = .(assay), layer = "data")
        hdWGCNA::TestSoftPowers(obj, networkType = .(network_type))
      }, where = environment()),
      deployment = "main",
      description = base::paste0("hdWGCNA ", g, ": set datExpr + test soft powers")
    )

    soft_power <- targets::tar_target_raw(
      name = sp_name,
      command = bquote({
        library(hdWGCNA)
        power_table <- hdWGCNA::GetPowerTable(.(as.name(pt_name)))
        pass <- power_table[["SFT.R.sq"]] >= .(sft_rsquared)
        if (any(pass, na.rm = TRUE)) {
          min(power_table[["Power"]][which(pass)])
        } else {
          max(power_table[["Power"]], na.rm = TRUE)
        }
      }, where = environment()),
      deployment = "main",
      description = base::paste0("hdWGCNA ", g, ": soft power (lowest tested power with SFT.R.sq >= ", sft_rsquared, ")")
    )

    network <- targets::tar_target_raw(
      name = net_name,
      command = bquote({
        library(WGCNA)
        library(hdWGCNA)
        WGCNA::enableWGCNAThreads(nThreads = .(n_threads))
        obj <- hdWGCNA::ConstructNetwork(.(as.name(pt_name)), soft_power = .(as.name(sp_name)), tom_name = .(g), tom_outdir = .(tom_outdir), overwrite_tom = .(overwrite_tom))
        obj <- hdWGCNA::ModuleEigengenes(obj, group.by.vars = .(patient_col), assay = .(assay))
        obj <- hdWGCNA::ModuleConnectivity(obj, group.by = .(clustering_col), group_name = .(g), assay = .(assay))
        obj <- hdWGCNA::ResetModuleNames(obj, new_name = .(paste0(g, "-M")))
        if (!is.null(.(trait_col)) && length(.(trait_groups)) > 0L) {
          for (tg in .(trait_groups)) {
            obj@meta.data[[paste0("is_", gsub("[^A-Za-z0-9]+", "", tg))]] <- as.numeric(obj@meta.data[[.(trait_col)]] == tg)
          }
          # single constant grouping so the correlation is computed once across all
          # of the group's cells ("all_cells") instead of falling back to Idents()
          # (which would split the heatmap by the object's cluster identities).
          obj@meta.data[[".mtc_group"]] <- .(g)
          obj <- hdWGCNA::ModuleTraitCorrelation(
            obj,
            traits = paste0("is_", gsub("[^A-Za-z0-9]+", "", .(trait_groups))),
            features = "hMEs", cor_method = .(mtc_cor_method),
            group.by = ".mtc_group",
            subset_by = .(clustering_col), subset_groups = .(g)
          )
        }
        obj
      }, where = environment()),
      deployment = "main",
      description = base::paste0("hdWGCNA ", g, ": network, module eigengenes, connectivity, module-trait correlation")
    )

    out <- list(powertest, soft_power, network)

    if (!is.null(trait_col)) {
      dmes <- targets::tar_target_raw(
        name = dme_name,
        command = bquote({
          library(WGCNA)
          library(hdWGCNA)
          obj <- .(as.name(net_name))
          md <- obj@meta.data
          cells <- md[[.(clustering_col)]] == .(g) & !is.na(md[[.(trait_col)]])
          comparisons <- utils::combn(.(trait_groups), 2, simplify = FALSE)
          dme_list <- lapply(comparisons, function(cc) {
            g1 <- cc[[1]]
            g2 <- cc[[2]]
            b1 <- rownames(md)[cells & md[[.(trait_col)]] == g1]
            b2 <- rownames(md)[cells & md[[.(trait_col)]] == g2]
            d <- tryCatch(
              hdWGCNA::FindDMEs(obj, barcodes1 = b1, barcodes2 = b2, test.use = .(dme_test), harmonized = .(dme_harmonized)),
              error = function(e) {
                warning("FindDMEs failed for ", g1, " vs ", g2, ": ", conditionMessage(e))
                NULL
              }
            )
            if (is.null(d) || nrow(d) == 0L) {
              return(NULL)
            }
            if (!"module" %in% colnames(d)) d$module <- rownames(d)
            rownames(d) <- NULL
            d$group1 <- g1
            d$group2 <- g2
            d$comparison <- paste(g1, "vs", g2)
            d
          })
          do.call(rbind, Filter(Negate(is.null), dme_list))
        }, where = environment()),
        deployment = deployment,
        resources = movable_resources,
        description = base::paste0("hdWGCNA ", g, ": differential module eigengenes (FindDMEs) for all pairwise trait comparisons")
      )
      out <- c(out, list(dmes))
    }

    if (isTRUE(run_enrichr)) {
      enrichr <- targets::tar_target_raw(
        name = enr_name,
        command = bquote({
          library(WGCNA)
          library(hdWGCNA)
          obj <- .(as.name(net_name))
          # Optionally restrict the GO/Enrichr analysis to a subset of modules
          # (`enrichr_modules`): the modules correlated with a clinical group
          # ("correlated"), the modules with a significant differential module
          # eigengene ("dmes"), or their union ("both"). Selected modules are kept;
          # every other non-grey module is relabelled "grey" and its factor level
          # dropped, so RunEnrichr (which iterates `levels(modules$module)` minus
          # grey) only queries the kept modules. This edits this target's local
          # copy of the object only, not the stored network or DME targets.
          run_dbs <- TRUE
          if (.(enrichr_restrict)) {
            keep_mods <- base::character(0)
            if (.(enrichr_need_corr)) {
              keep_mods <- base::union(
                keep_mods,
                scitargets::hdwgcna_associated_modules(obj, p_cutoff = .(enrichr_corr_cutoff))
              )
            }
            dmes <- .(if (enrichr_need_dme) as.name(dme_name) else NULL)
            if (!base::is.null(dmes) && base::nrow(dmes) > 0L) {
              sig <- dmes[!base::is.na(dmes$p_val_adj) & dmes$p_val_adj < .(enrichr_dme_cutoff), , drop = FALSE]
              keep_mods <- base::union(keep_mods, base::unique(base::as.character(sig$module)))
            }
            keep_mods <- base::setdiff(keep_mods, "grey")
            if (base::length(keep_mods) == 0L) {
              warning("No modules selected for the GO analysis (enrichr_modules = \"",
                .(enrichr_modules), "\"); the Enrichr table is empty.")
              run_dbs <- FALSE
            } else {
              mods_df <- hdWGCNA::GetModules(obj)
              mod_chr <- base::as.character(mods_df$module)
              mod_chr[!(mod_chr %in% keep_mods)] <- "grey"
              mods_df$module <- base::factor(mod_chr, levels = c(keep_mods, "grey"))
              obj <- hdWGCNA::SetModules(obj, mods_df)
            }
          }
          # Run each Enrichr database separately and stack the per-database tables
          # (GetEnrichrTable tags each row with its `db`). A failed database (e.g.
          # the Enrichr API being unreachable) is dropped with a warning.
          if (!run_dbs) {
            base::data.frame()
          } else {
            do.call(rbind, Filter(Negate(is.null), lapply(.(enrichr_dbs), function(db) {
              tryCatch({
                o <- hdWGCNA::RunEnrichr(obj, dbs = db, max_genes = .(enrichr_max_genes))
                hdWGCNA::GetEnrichrTable(o)
              }, error = function(e) {
                warning("RunEnrichr failed for ", db, ": ", conditionMessage(e))
                NULL
              })
            })))
          }
        }, where = environment()),
        deployment = "main",
        description = base::paste0("hdWGCNA ", g, ": Enrichr enrichment (", enrichr_modules,
          " modules), one call per database (", base::paste(enrichr_dbs, collapse = ", "), ")")
      )
      out <- c(out, list(enrichr))
    }

    # Per-group long-format module composition (small data frame). Kept as a
    # separate light target so the aggregate workbook can be assembled without
    # holding several full network objects in memory at once.
    if (!is.null(module_table_file)) {
      modcomp <- targets::tar_target_raw(
        name = base::paste0("wgcna_", suffix, "_modcomp"),
        command = bquote({
          library(WGCNA)
          library(hdWGCNA)
          scitargets::hdwgcna_module_composition(.(as.name(net_name)))
        }, where = environment()),
        deployment = "main",
        description = base::paste0("hdWGCNA ", g, ": long-format module composition (GetModules)")
      )
      out <- c(out, list(modcomp))
    }

    # Per-group module -> gene-vector list (the building block of the nested
    # `wgcna_module_gene_list`). Light target (just the GetModules() assignment).
    modgenes <- targets::tar_target_raw(
      name = base::paste0("wgcna_", suffix, "_modgenes"),
      command = bquote({
        library(hdWGCNA)
        scitargets::hdwgcna_module_gene_list(.(as.name(net_name)))
      }, where = environment()),
      deployment = "main",
      description = base::paste0("hdWGCNA ", g, ": module -> gene-vector list")
    )
    out <- c(out, list(modgenes))

    # Per-cell module scores (Seurat::AddModuleScore) for this group's modules,
    # computed across ALL cells of the object (every cell type), carrying each
    # cell's cell type and clinical condition so the report can group by them.
    # Loads the full network object -> "main".
    modscore <- targets::tar_target_raw(
      name = base::paste0("wgcna_", suffix, "_modscore"),
      command = bquote({
        library(Seurat)
        library(hdWGCNA)
        scitargets::hdwgcna_module_scores(
          .(as.name(net_name)),
          trait_col = .(trait_col),
          clustering_col = .(clustering_col),
          assay = .(assay)
        )
      }, where = environment()),
      deployment = "main",
      description = base::paste0("hdWGCNA ", g, ": per-cell module scores across all cells (Seurat::AddModuleScore)")
    )
    out <- c(out, list(modscore))
    out
  })

  # One workbook per tar_hdwgcna() call, with a worksheet per cell group.
  module_tab <- if (!is.null(module_table_file)) {
    modcomp_names <- base::paste0("wgcna_", vapply(group, .hdwgcna_suffix, character(1)), "_modcomp")
    tabs_expr <- base::as.call(c(base::quote(base::list), lapply(modcomp_names, as.name)))
    list(targets::tar_target_raw(
      name = "wgcna_module_composition",
      command = bquote({
        tabs <- .(tabs_expr)
        base::names(tabs) <- .(group)
        out_file <- .(module_table_file)
        base::dir.create(base::dirname(out_file), recursive = TRUE, showWarnings = FALSE)
        wb <- openxlsx2::wb_workbook()
        for (nm in base::names(tabs)) {
          # Excel sheet names: <= 31 chars, no : \\ / ? * [ ]
          sheet <- base::substr(base::gsub("[:\\\\/?*\\[\\]]", "_", nm), 1L, 31L)
          wb <- openxlsx2::wb_add_worksheet(wb, sheet)
          wb <- openxlsx2::wb_add_data(wb, sheet = sheet, x = tabs[[nm]])
        }
        openxlsx2::wb_save(wb, out_file, overwrite = TRUE)
        out_file
      }, where = environment()),
      format = "file",
      deployment = "main",
      description = base::paste0("hdWGCNA: module composition workbook (one sheet per cell group) -> ", module_table_file)
    ))
  } else {
    list()
  }

  # Single nested module gene list across this call's groups:
  #   ls[[<cell type>]][[<module>]] = character vector of genes.
  # Assembled from the light per-group `wgcna_<group>_modgenes` targets. One per
  # tar_hdwgcna() call (fixed name `wgcna_module_gene_list`); if you make several
  # calls, set `module_gene_list = FALSE` on all but one to avoid a name clash.
  gene_list_tab <- if (isTRUE(module_gene_list)) {
    modgenes_names <- base::paste0("wgcna_", vapply(group, .hdwgcna_suffix, character(1)), "_modgenes")
    genes_expr <- base::as.call(c(base::quote(base::list), lapply(modgenes_names, as.name)))
    list(targets::tar_target_raw(
      name = "wgcna_module_gene_list",
      command = bquote({
        ls <- .(genes_expr)
        base::names(ls) <- .(group)
        ls
      }, where = environment()),
      deployment = "main",
      description = "hdWGCNA: nested module gene list ls[[cell_type]][[module]] = genes"
    ))
  } else {
    list()
  }

  c(prep, unlist(per_group, recursive = FALSE), module_tab, gene_list_tab)
}
