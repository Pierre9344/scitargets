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

# Insert a scope suffix before a path's extension, so several clustering
# columns in one tar_hdwgcna() call get one prep file each rather than
# overwriting one another: "out/prep.qs2" + "_0.3" -> "out/prep_0.3.qs2".
.hdwgcna_suffix_path <- function(path, suffix) {
  if (base::is.null(path) || !nzchar(suffix)) {
    return(path)
  }
  ext <- base::regmatches(path, base::regexpr("\\.[^.\\\\/]+$", path))
  base::paste0(base::sub("\\.[^.\\\\/]+$", "", path), suffix, ext)
}

#' Is a file a usable `wgcna_prep` object?
#'
#' [tar_hdwgcna()] can be pointed at a pre-built metacell object with its
#' `prep_path` argument, so that the expensive preparation runs once on a
#' machine with enough memory and every other machine reads the result. This is
#' the check it applies before deciding to reuse such a file.
#'
#' A file is usable when all of the following hold. Each failure is reported as
#' a message and makes the function return `FALSE`, so a broken or half-copied
#' file leads to a rebuild rather than to an obscure error several steps later.
#'
#' 1. The path exists and is not empty.
#' 2. [as_seurat()] can read it, which also means the extension is one of the
#'    supported ones and the contents are a Seurat object.
#' 3. The object carries the hdWGCNA experiment named `wgcna_name`, i.e.
#'    `SetupForWGCNA()` was run on it with that name.
#' 4. That experiment has a metacell object attached, i.e.
#'    `MetacellsByGroups()` was run.
#'
#' The file is read in full, because a truncated copy is exactly the failure
#' this is meant to catch and only a real read finds it. That costs one extra
#' read of the object per pipeline run.
#'
#' @param path Path to the candidate file.
#' @param wgcna_name Name of the hdWGCNA experiment to look for.
#' @returns `TRUE` or `FALSE`.
#' @export
hdwgcna_prep_is_valid <- function(path, wgcna_name = "hdwgcna") {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(FALSE)
  }
  if (!base::file.exists(path)) {
    base::message("hdwgcna_prep_is_valid(): no file at '", path, "'.")
    return(FALSE)
  }
  if (base::isTRUE(base::file.size(path) == 0)) {
    base::message("hdwgcna_prep_is_valid(): '", path, "' is empty.")
    return(FALSE)
  }
  obj <- base::tryCatch(as_seurat(path), error = function(e) {
    base::message("hdwgcna_prep_is_valid(): could not read '", path, "': ",
      base::conditionMessage(e))
    NULL
  })
  if (base::is.null(obj)) {
    return(FALSE)
  }
  if (!wgcna_name %in% base::names(obj@misc)) {
    base::message("hdwgcna_prep_is_valid(): '", path, "' holds a Seurat object ",
      "with no '", wgcna_name, "' hdWGCNA experiment (SetupForWGCNA was not run ",
      "with that wgcna_name).")
    return(FALSE)
  }
  has_metacells <- base::tryCatch(
    !base::is.null(hdWGCNA::GetMetacellObject(obj, wgcna_name = wgcna_name)),
    error = function(e) FALSE
  )
  if (!base::isTRUE(has_metacells)) {
    base::message("hdwgcna_prep_is_valid(): '", path, "' has no metacell object ",
      "for '", wgcna_name, "' (MetacellsByGroups was not run).")
    return(FALSE)
  }
  TRUE
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
#' - `wgcna_<group>_powertest` — `SetDatExpr` for the group + `TestSoftPowers`,
#'   returning the scale-free fit table and the chosen power rather than the
#'   Seurat object.
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
#'
#'   When `clustering_col` holds several columns, `group` must be a **list** of
#'   the same length, one character vector per column, since a clustering at one
#'   resolution has different labels from another. A single character vector is
#'   recycled to every column, which is only sensible when the columns really do
#'   share their labels.
#' @param create_prep Whether this call creates the shared `wgcna_prep` target.
#'   `TRUE` (default) builds it; set `FALSE` on the additional `tar_hdwgcna()`
#'   calls that reuse the same `wgcna_prep` (exactly one call must use `TRUE`).
#'   When `FALSE`, `input_obj` is not required.
#' @param prep_path Where the prepared metacell object is cached, or `NULL` to
#'   keep it in the `targets` store as an ordinary object target. When set
#'   (the default, `"./out/seurat/wgcna_prep.qs2"`), `wgcna_prep` becomes a
#'   `format = "file"` target whose value is this path: it reuses the file when
#'   [hdwgcna_prep_is_valid()] accepts it, and otherwise builds the object,
#'   writes it with [save_seurat()] and returns the path. The extension picks
#'   the format (`.qs2`/`.qs`, `.rds`, `.RData`/`.rda`) and is checked when the
#'   pipeline is defined; note that `.qs` is written by `qs2` here too, so it is
#'   a qs2-format file rather than one the old `qs` package can read. See the
#'   "Caching the prepared object" section for why this exists and how to build
#'   the file by hand. With several `clustering_col` values the scope suffix is
#'   inserted before the extension, so `"out/prep.qs2"` becomes
#'   `out/prep_0.3.qs2` and `out/prep_0.6.qs2`.
#' @param input_obj Name of the upstream Seurat-object target to start from
#'   (required when `create_prep = TRUE`). That target may hold the Seurat
#'   object itself, or the PATH to a file holding it, as a `format = "file"`
#'   target does; [as_seurat()] accepts either and says which file types it
#'   reads.
#' @param wgcna_name Name of the hdWGCNA experiment (`SetupForWGCNA`).
#' @param assay Assay used throughout (default `"RNA"`, log-normalised).
#' @param clustering_col Metadata column holding the cell groups; used as
#'   `ident.group` for the metacells and as `group.by`/`subset_by` downstream.
#'
#'   May be a character **vector**, to run the whole analysis over several
#'   groupings of the same object in one call, for example two clustering
#'   resolutions and an Azimuth annotation. Each column gets its own
#'   `wgcna_prep`, since the metacells depend on the grouping, and its own set of
#'   per-group targets. Every generated name is then suffixed to keep the scopes
#'   apart; the suffix comes from the vector's names when it has any, otherwise
#'   from the column name itself:
#'
#'   ```r
#'   clustering_col = c(`0.3` = "clusters_0.3", `0.6` = "clusters_0.6")
#'   # -> wgcna_prep_0.3, wgcna_1_0.3, ..., wgcna_prep_0.6, wgcna_1_0.6, ...
#'   ```
#'
#'   A single unnamed column produces unsuffixed names, exactly as before, so
#'   existing pipelines are unaffected.
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
#' @param delete_tom Whether to delete the topological overlap matrix as soon as
#'   `ConstructNetwork()` has finished with it (default `TRUE`). The TOM is a
#'   dense gene-by-gene matrix, roughly 3 GB per cell group here, and
#'   `ConstructNetwork()` is the only thing that reads it: the modules,
#'   eigengenes, connectivity and every downstream target are derived inside
#'   that call, so the file is write-only for this pipeline.
#'
#'   Set `FALSE` if you intend to call the hdWGCNA functions that do re-read it,
#'   namely `HubGeneNetworkPlot()`, `ModuleNetworkPlot()` and
#'   `RunModuleUMAP()` / `ModuleUMAPPlot()`. Nothing in this package uses them.
#'   Deleting only removes `<tom_outdir>/<tom_name>_TOM.rda` and any
#'   `<tom_name>_block.N.rda` left behind by a multi-block run; the fitted
#'   network object is untouched.
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
#' @section Caching the prepared object:
#'
#' `wgcna_prep` is the memory peak of this factory, and by a wide margin. It
#' holds the whole Seurat object, adds a normalised assay layer and a dense
#' `scale.data` over the variable features, and only then reduces to metacells.
#' On a dataset of ~50k cells that is comfortably into double-digit GB, which is
#' more than a small server or a memory-capped container has, while the steps
#' below it work on metacells and are far cheaper.
#'
#' `prep_path` exists so that peak has to be paid once, on a machine that can
#' afford it, rather than on every machine that wants the results. When it is
#' set, `wgcna_prep` becomes a `format = "file"` target that:
#'
#' 1. checks the file with [hdwgcna_prep_is_valid()] and, if it passes, returns
#'    the path without building anything;
#' 2. otherwise runs the preparation, writes it with [save_seurat()] and returns
#'    the path.
#'
#' So the usual workflow is to run the pipeline once where there is memory, copy
#' the resulting file across with the rest of `out/`, and run everything else
#' anywhere. **The file wins over the pipeline**: an upstream change does not
#' refresh it, exactly so that a machine which cannot rebuild it does not try.
#' Delete the file to force a rebuild.
#'
#' The object is a plain Seurat object, so it can also be built outside the
#' pipeline. This is what the target does, and a hand-made file only has to
#' match it:
#'
#' ```r
#' library(Seurat); library(hdWGCNA); library(WGCNA)
#' # load the seurat object and basic setup
#' obj <- scitargets::as_seurat("out/seurat/seurat_obj.qs2")
#' DefaultAssay(obj) <- "RNA"
#' obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
#' obj <- NormalizeData(obj, verbose = FALSE)
#' obj <- FindVariableFeatures(obj, verbose = FALSE)
#' obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
#'
#' # WGCNA setup (see https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html#set-up-seurat-object-for-wgcna)
#' obj <- SetupForWGCNA(obj,
#'  gene_select = "fraction",     # gene_select
#'   fraction = 0.05,             # fraction
#'   wgcna_name = "hdwgcna"       # wgcna_name
#'   )
#' obj <- MetacellsByGroups(obj,
#'   group.by    = c("predicted.celltype.l2", "patient_id"), # clustering_col and patient_col
#'   ident.group = "predicted.celltype.l2",                  # metacell_group_by
#'   reduction   = "harmony",                                # reduction
#'    k = 25,                                                # metacell_k
#'    max_shared = 10                                        # metacell_max_shared
#'    )
#' obj <- NormalizeMetacells(obj)
#' obj <- ScaleMetacells(obj, features = VariableFeatures(obj))
#' obj <- RunPCAMetacells(obj, features = VariableFeatures(obj))
#' obj <- RunHarmonyMetacells(
#'    obj,
#'    group.by.vars = "patient_id"      # patient_cols
#'    )
#' obj <- RunUMAPMetacells(obj,
#'    reduction = "harmony",            # reduction
#'    dims = 1:15                       # metacells_dims
#'    )
#' scitargets::save_seurat(obj, "./out/seurat/wgcna_prep.qs2")
#' ```
#'
#' Every argument above must match the `tar_hdwgcna()` call that will read the
#' file (`assay`, `gene_select`, `fraction`, `wgcna_name`, `metacell_group_by`,
#' `clustering_col`, `reduction`, `metacell_k`, `metacell_max_shared`,
#' `metacell_dims`, `patient_col`), because none of them is re-checked: only the
#' structural conditions in [hdwgcna_prep_is_valid()] are. A file prepared with
#' a different `fraction` is accepted and silently changes the gene set the
#' networks are built on.
#'
#' Set `prep_path = NULL` to restore the previous behaviour, where the prepared
#' object is an ordinary target kept in the `targets` store.
#'
#' @returns A list of `targets` objects to splice into a `targets` pipeline. With
#'   `<g>` the lower-cased, sanitised label of each `group` element (e.g. `"MAIT"`
#'   becomes `mait`, `"CD8 TEM"` becomes `cd8_tem`), the following steps are created:
#'   - `wgcna_prep`: the shared metacell preparation (`deployment = "main"`).
#'     Created only when `create_prep = TRUE` (the default); omit it (set
#'     `create_prep = FALSE`) on the other `tar_hdwgcna()` calls that reuse the
#'     same `wgcna_prep`. With `prep_path` set (the default) this is a
#'     `format = "file"` target whose value is that path rather than the object;
#'     see the "Caching the prepared object" section.
#'   - `wgcna_<g>_powertest`: `SetDatExpr` for the group + `TestSoftPowers`
#'     (one per group; `deployment = "main"`, since it loads the whole Seurat
#'     object). Returns a small `list`, not the Seurat object: `power_table`
#'     (the scale-free fit table `hdWGCNA::PlotSoftPowers()` draws),
#'     `soft_power`, `reached_threshold`, `sft_rsquared` and `group`. Keeping
#'     the object here cost several GB of store per cell group to carry a few
#'     hundred numbers; the network step rebuilds `datExpr` from `wgcna_prep`
#'     instead, which is cheap next to `ConstructNetwork`.
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
  prep_path = "./out/seurat/wgcna_prep.qs2",
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
  delete_tom = TRUE,
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
  module_gene_list = TRUE,
  name_suffix = ""
) {
  # ---- several groupings in one call -------------------------------------
  # Each column is its own scope: the metacells depend on the grouping, so
  # nothing is shared between them and the work is simply done once per column.
  # Dispatching here, before any validation, keeps the single-scope path below
  # byte-for-byte what it was, so an existing single-column call still produces
  # unsuffixed target names and nothing already computed is invalidated.
  # Several columns, or a single NAMED one: a name is how the caller asks for a
  # suffix, so one scope out of a factory still gets `wgcna_prep_0.3` rather
  # than the bare legacy name. An unnamed single column stays on the legacy path
  # and keeps its unsuffixed names.
  if (length(clustering_col) > 1L ||
    (length(clustering_col) == 1L && !base::is.null(base::names(clustering_col)))) {
    if (nzchar(name_suffix)) {
      stop("`name_suffix` cannot be combined with several `clustering_col` values; ",
        "the suffixes are derived from `clustering_col`, or from its names when it has any.")
    }
    # Names when supplied, otherwise the column itself, so an unnamed
    # c("clusters_0.3", "clusters_0.6") still yields readable, distinct names.
    keys <- base::names(clustering_col)
    if (base::is.null(keys) || base::any(!nzchar(keys))) keys <- clustering_col
    # Dots are legal in target names and the surrounding pipeline already uses
    # them for resolutions (clusters_0.6, cluster_percent_0.6), so keep them
    # rather than running the group sanitiser, which would give wgcna_prep_0_3.
    suffixes <- base::tolower(base::gsub("[^A-Za-z0-9._]+", "_", keys))
    suffixes <- base::paste0("_", base::gsub("^[._]+|[._]+$", "", suffixes))
    if (base::anyDuplicated(suffixes)) {
      stop("`clustering_col` must give distinct target-name suffixes; got: ",
        base::paste(suffixes, collapse = ", "))
    }
    # One clustering has different labels from another, so `group` must say
    # which labels belong to which column.
    if (!base::is.list(group)) group <- base::rep(base::list(group), length(clustering_col))
    if (length(group) != length(clustering_col)) {
      stop("`group` must be a list with one element per `clustering_col` (",
        length(clustering_col), "), or a single character vector to reuse for all; got ",
        length(group), ".")
    }
    out <- base::lapply(base::seq_along(clustering_col), function(i) {
      tar_hdwgcna(
        group = group[[i]],
        create_prep = create_prep,
        # One prep file per scope, for the same reason as tom_outdir below: the
        # metacells depend on the clustering column, so two resolutions sharing
        # one path would each read the other's object.
        prep_path = .hdwgcna_suffix_path(prep_path, suffixes[i]),
        input_obj = input_obj,
        wgcna_name = wgcna_name,
        assay = assay,
        clustering_col = base::unname(clustering_col[[i]]),
        patient_col = patient_col,
        # Left NULL so it re-defaults to c(<this scope's column>, patient_col)
        # rather than being resolved once from the whole vector.
        metacell_group_by = metacell_group_by,
        reduction = reduction,
        gene_select = gene_select,
        fraction = fraction,
        metacell_k = metacell_k,
        metacell_max_shared = metacell_max_shared,
        metacell_dims = metacell_dims,
        network_type = network_type,
        sft_rsquared = sft_rsquared,
        # ConstructNetwork writes the TOM as <tom_outdir>/<group>, and two
        # resolutions can both have a cluster called "1", so each scope needs its
        # own directory or the second silently overwrites the first's TOM.
        # sub() on the trailing slash first: tom_outdir usually ends in "/", and
        # file.path() would then produce "./data/hdwgcna//0.3". Harmless, but it
        # shows up in every path the step prints.
        tom_outdir = base::file.path(
          base::sub("/+$", "", tom_outdir), base::sub("^_", "", suffixes[i])
        ),
        overwrite_tom = overwrite_tom,
        delete_tom = delete_tom,
        n_threads = n_threads,
        deployment = deployment,
        controller = controller,
        trait_col = trait_col,
        trait_groups = trait_groups,
        mtc_cor_method = mtc_cor_method,
        dme_test = dme_test,
        dme_harmonized = dme_harmonized,
        run_enrichr = run_enrichr,
        enrichr_dbs = enrichr_dbs,
        enrichr_max_genes = enrichr_max_genes,
        enrichr_modules = enrichr_modules,
        enrichr_corr_cutoff = enrichr_corr_cutoff,
        enrichr_dme_cutoff = enrichr_dme_cutoff,
        # One workbook per scope, or each would overwrite the last.
        module_table_file = if (base::is.null(module_table_file)) NULL else {
          base::sub("(\\.xlsx)$", base::paste0(suffixes[i], "\\1"), module_table_file)
        },
        module_gene_list = module_gene_list,
        name_suffix = suffixes[i]
      )
    })
    return(base::unlist(out, recursive = FALSE))
  }
  if (!is.character(name_suffix) || length(name_suffix) != 1L || is.na(name_suffix)) {
    stop("`name_suffix` must be a single string.")
  }
  if (!is.logical(delete_tom) || length(delete_tom) != 1L || is.na(delete_tom)) {
    stop("`delete_tom` must be a single TRUE/FALSE.")
  }
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
  # Checked HERE rather than at build time: a path the writer cannot handle
  # should fail when the pipeline is defined, not after the half-hour of work
  # that was supposed to fill it.
  if (!is.null(prep_path)) {
    if (!is.character(prep_path) || length(prep_path) != 1L || !nzchar(prep_path)) {
      stop("`prep_path` must be NULL or a single non-empty file path.")
    }
    .ext <- base::tolower(base::sub(".*\\.", "", base::basename(prep_path)))
    if (!.ext %in% seurat_file_formats()) {
      stop("`prep_path` must end in one of .qs2, .qs, .rds, .RData or .rda; got '",
        base::basename(prep_path), "'.")
    }
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
  prep_name <- base::paste0("wgcna_prep", name_suffix)

  # The preparation itself, identical in both modes. Kept as one expression so
  # the in-store and the on-disk variants cannot drift apart.
  prep_body <- bquote({
    WGCNA::enableWGCNAThreads(nThreads = .(n_threads))
    # as_seurat(): the upstream target may hold the object itself OR a path to
    # it, as a `format = "file"` target does. See ?as_seurat.
    obj <- scitargets::as_seurat(.(as.name(input_obj)))
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
  }, where = environment())

  prep <- if (isTRUE(create_prep)) list(targets::tar_target_raw(
    name = prep_name,
    command = if (base::is.null(prep_path)) {
      bquote({
        library(WGCNA)
        library(hdWGCNA)
        .(prep_body)
      }, where = environment())
    } else {
      bquote({
        library(WGCNA)
        library(hdWGCNA)
        .prep_path <- .(prep_path)
        # A usable file WINS over rebuilding, which is the whole point: the
        # preparation holds the full Seurat object plus a dense scale.data and
        # peaks well above what a small machine has, so it is run once where
        # there is room and the file is copied. It also means a change upstream
        # does NOT refresh this object: delete the file to force a rebuild.
        if (scitargets::hdwgcna_prep_is_valid(.prep_path, .(wgcna_name))) {
          base::message("hdWGCNA: reusing the prepared object at '", .prep_path,
            "' (delete it to rebuild).")
          return(.prep_path)
        }
        base::message("hdWGCNA: building the prepared object into '", .prep_path, "'.")
        obj <- .(prep_body)
        scitargets::save_seurat(obj, .prep_path)
        .prep_path
      }, where = environment())
    },
    # format = "file" in the on-disk mode: the target's VALUE is the path, so
    # targets hashes the object file itself and a file copied in from another
    # machine is noticed. as_seurat() downstream accepts either form, so the
    # rest of the factory is the same in both modes.
    format = if (base::is.null(prep_path)) {
      targets::tar_option_get("format")
    } else {
      "file"
    },
    deployment = "main",
    description = if (base::is.null(prep_path)) {
      "hdWGCNA: shared metacell preparation (RNA: join layers, normalize, scale, SetupForWGCNA, metacells)"
    } else {
      base::paste0(
        "hdWGCNA: shared metacell preparation, cached at ", prep_path,
        " (RNA: join layers, normalize, scale, SetupForWGCNA, metacells). ",
        "Reuses the file when it is valid, otherwise builds and writes it."
      )
    }
  )) else list()

  per_group <- lapply(group, function(g) {
    suffix <- .hdwgcna_suffix(g)
    pt_name <- base::paste0("wgcna_", suffix, "_powertest", name_suffix)
    sp_name <- base::paste0("wgcna_", suffix, "_soft_power", name_suffix)
    net_name <- base::paste0("wgcna_", suffix, name_suffix)
    dme_name <- base::paste0("wgcna_", suffix, "_dmes", name_suffix)
    enr_name <- base::paste0("wgcna_", suffix, "_enrichr", name_suffix)

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
        # as_seurat(): wgcna_prep holds the object itself, or the PATH to it when
        # tar_hdwgcna() was given a `prep_path`. Either form works here.
        obj <- hdWGCNA::SetDatExpr(scitargets::as_seurat(.(as.name(prep_name))), group_name = .(g), group.by = .(clustering_col), assay = .(assay), layer = "data")
        obj <- hdWGCNA::TestSoftPowers(obj, networkType = .(network_type))
        # Return the SCALE-FREE FIT TABLE, not the Seurat object it came in.
        # TestSoftPowers only stashes a small table on a copy of the whole
        # object, so storing the object costs several GB per cell group (it was
        # ~3.5 GB each, ~70 GB across the groups) to keep a few hundred numbers.
        # Nothing downstream needs the object: the network step rebuilds datExpr
        # from wgcna_prep, which is cheap next to ConstructNetwork.
        power_table <- hdWGCNA::GetPowerTable(obj)
        pass <- power_table[["SFT.R.sq"]] >= .(sft_rsquared)
        base::list(
          # The table hdWGCNA::PlotSoftPowers() draws, so a report can redraw it
          # from this target alone.
          power_table = power_table,
          # Lowest power reaching the scale-free fit threshold, else the highest
          # tested; `reached_threshold` says which of the two happened, because
          # the fallback is a very different claim about the network.
          soft_power = if (base::any(pass, na.rm = TRUE)) {
            base::min(power_table[["Power"]][base::which(pass)])
          } else {
            base::max(power_table[["Power"]], na.rm = TRUE)
          },
          reached_threshold = base::any(pass, na.rm = TRUE),
          sft_rsquared = .(sft_rsquared),
          group = .(g)
        )
      }, where = environment()),
      deployment = "main",
      description = base::paste0("hdWGCNA ", g, ": set datExpr + test soft powers; returns the scale-free fit table and the chosen power, not the Seurat object")
    )

    soft_power <- targets::tar_target_raw(
      name = sp_name,
      command = bquote({
        .(as.name(pt_name))[["soft_power"]]
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
        # datExpr is rebuilt here rather than inherited from the power-test
        # target, which now returns only its fit table. SetDatExpr is a subset
        # and is cheap next to ConstructNetwork; carrying the object between the
        # two steps instead cost several GB of store per cell group.
        # as_seurat(): as in the power-test target, wgcna_prep may hold the
        # object or the path to it.
        obj <- hdWGCNA::SetDatExpr(scitargets::as_seurat(.(as.name(prep_name))), group_name = .(g), group.by = .(clustering_col), assay = .(assay), layer = "data")
        obj <- hdWGCNA::ConstructNetwork(obj, soft_power = .(as.name(sp_name)), tom_name = .(g), tom_outdir = .(tom_outdir), overwrite_tom = .(overwrite_tom))
        if (.(delete_tom)) {
          # ConstructNetwork() has already derived the modules from the TOM, and
          # nothing downstream re-reads it, so the file is dead weight from here
          # on: a dense gene-by-gene matrix, gigabytes per cell group. Deleted
          # immediately rather than at the end of the target, so the disk is
          # freed even if a later step of this same target fails.
          # Exact paths rather than a regex, because a group name may contain
          # spaces or regex metacharacters ("CD4 TCM").
          .tom_files <- base::c(
            base::file.path(.(tom_outdir), base::paste0(.(g), "_TOM.rda")),
            base::Sys.glob(base::file.path(.(tom_outdir), base::paste0(.(g), "_block.*.rda")))
          )
          .tom_files <- base::unique(.tom_files[base::file.exists(.tom_files)])
          if (base::length(.tom_files)) {
            .freed <- base::sum(base::file.size(.tom_files), na.rm = TRUE)
            base::file.remove(.tom_files)
            base::message(
              "hdWGCNA ", .(g), ": removed ", base::length(.tom_files),
              " TOM file(s), freeing ",
              base::round(.freed / 1024^3, 2), " GB (delete_tom = TRUE)"
            )
          }
        }
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
        name = base::paste0("wgcna_", suffix, "_modcomp", name_suffix),
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
      name = base::paste0("wgcna_", suffix, "_modgenes", name_suffix),
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
      name = base::paste0("wgcna_", suffix, "_modscore", name_suffix),
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
    modcomp_names <- base::paste0(
      "wgcna_", vapply(group, .hdwgcna_suffix, character(1)), "_modcomp", name_suffix
    )
    tabs_expr <- base::as.call(c(base::quote(base::list), lapply(modcomp_names, as.name)))
    list(targets::tar_target_raw(
      name = base::paste0("wgcna_module_composition", name_suffix),
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
    modgenes_names <- base::paste0(
      "wgcna_", vapply(group, .hdwgcna_suffix, character(1)), "_modgenes", name_suffix
    )
    genes_expr <- base::as.call(c(base::quote(base::list), lapply(modgenes_names, as.name)))
    list(targets::tar_target_raw(
      name = base::paste0("wgcna_module_gene_list", name_suffix),
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
