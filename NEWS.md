# scitargets 1.8.0

- New `as_seurat()`: takes a Seurat object **or a path to a file holding one**, and returns the object. It reads `.qs2`/`.qs` with `qs2::qs_read()`, `.rds` with `readRDS()` and `.RData`/`.rda` with `load()`, and errors with a clear message on anything else, on a file whose contents are not a Seurat object, or on an `.RData` holding more than one object (`save()`/`load()` store the variable's NAME, so a multi-object file has no unambiguous answer).

- `tar_hdwgcna()` and `run_dea()` accept either form for their Seurat input. A pipeline whose merged object is too expensive to rebuild on every machine can write it to disk and declare its target `format = "file"`, so the normalisation runs once where there is enough memory and the other machines read the result; the target's value is then a path rather than an object, and `wgcna_prep` used to fail on it with `no applicable method for 'DefaultAssay<-' applied to an object of class "character"`. A target that still holds the object passes straight through, so nothing changes for existing pipelines.

- `as_seurat()` deliberately does **not** cache. A pipeline whose branches each need the object should memoise on its own side, so that the choice between paying one read per branch and holding the object in memory stays with the pipeline.

- The `cluster_to_use` argument of the `azimuth_annot_pbmc` function is now deprecated.

- `tar_demultiplex_hto` was modified to:
  - not create a new steps with the azimuth annotation. If needed, the annotations will be added directly to the `_singlets` or `_feat_removed_singlets` steps depending on whether features need to be removed or not.

# scitargets 1.7.0

- `azimuth_annot_pbmc()` now uses the **installed** PBMC reference instead of downloading it. `Azimuth:::LoadReference()` treats its `path` as a URL unless the string is an existing local directory, so the previous `reference = "pbmcref"` always fetched from seurat.nygenome.org even when `pbmcref.SeuratData` was installed and loading fine. On a host that cannot reach that server the annotation failed outright: on an HPC compute node behind a proxy answering `CONNECT tunnel failed, response 403`, a pipeline died at the annotation step after 38 minutes of upstream work, with the reference sitting installed in the library the whole time. The default `reference = NULL` now resolves `system.file("azimuth", package = "pbmcref.SeuratData")` and passes that directory, which takes `LoadReference()`'s local branch; it falls back to the name `"pbmcref"` when the package is absent, so behaviour is unchanged where the package was never installed. Anywhere else this is simply faster and reproducible, since the annotation no longer depends on a remote host being up. Pass `reference` explicitly to use a reference from another directory.

- `azimuth_annot_pbmc()` can also take the **homolog table** from a local file, via the new `homolog_table` argument. Fixing the reference alone is not enough on a host without access to seurat.nygenome.org: `RunAzimuth.Seurat()` then calls `ConvertGeneNames()` with the table's URL **hardcoded**, and exposes no argument for it, so the annotation fails a second time a few seconds later. `ConvertGeneNames()` itself does accept a local file, so when a table is configured this function swaps that one function out of Azimuth's namespace for the duration of the call and restores it on exit, including on error. This is a workaround for a hardcoded URL rather than a supported extension point, and it is skipped entirely when no table is configured, leaving the default behaviour untouched. The path is resolved from the argument, then `getOption("scitargets.azimuth_homologs")`, then the `AZIMUTH_HOMOLOGS` environment variable; a configured path that does not exist is ignored rather than raising, so a stale setting degrades to the old download instead of breaking the run. Fetch the file from `https://seurat.nygenome.org/azimuth/references/homologs.rds` (3.4 MB) on a machine with access.

- New exported `local_homolog_table()`, returning the configured homolog table or `NULL`, for the same preflight purpose as `local_pbmcref()`.

- `hdwgcna_report_lines()` gains `name_suffix`, matching the argument of the same name on `tar_hdwgcna()`. Since 1.6.0 a multi-scope `tar_hdwgcna()` call suffixes every target it generates (`wgcna_1_powertest_0.6`), but the report helper still built the unsuffixed names, so a report of such a pipeline could not read its own targets. The suffix is appended last, exactly as the factory appends it, so `name_suffix = "_0.6"` reads `wgcna_1_0.6`, `wgcna_1_powertest_0.6` and so on. The default is `""`, which reproduces the previous names, so single-scope pipelines are unaffected. A wrong suffix fails on a missing target rather than quietly reporting a different grouping's network.

- New exported `local_pbmcref()`, returning the installed reference directory or `NULL`. It checks for **both** `ref.Rds` and `idx.annoy` rather than only the directory, because `LoadReference()` needs exactly those two and a partial install otherwise fails later with a more obscure message. Useful as a preflight check: a pipeline can confirm the reference is usable before spending an hour reaching the step that needs it.

# scitargets 1.6.0

- `tar_hdwgcna()` now accepts a **vector** of `clustering_col`, running the whole analysis over several groupings of the same object in one call (for example two clustering resolutions and an Azimuth annotation). `group` becomes a list with one character vector per column, since a clustering at one resolution has different labels from another; a single character vector is recycled to every column. Each column is a separate scope with its own `wgcna_prep`, because the metacells depend on the grouping, and its own set of per-group targets. Every generated name is suffixed to keep the scopes apart, the suffix taken from the vector's names when it has any and from the column name otherwise, so `clustering_col = c(\`0.3\` = "clusters_0.3", \`0.6\` = "clusters_0.6")` yields `wgcna_prep_0.3`, `wgcna_1_powertest_0.3`, `wgcna_1_0.3`, `wgcna_module_composition_0.3` and the matching set for `0.6`. Dots are preserved in the suffix rather than sanitised to underscores, matching the `clusters_0.6` convention such metadata usually follows. Each scope also gets its own `tom_outdir` subdirectory and its own `module_table_file`, since two resolutions can both contain a cluster called `"1"` and would otherwise overwrite each other's TOM and workbook. A single **unnamed** `clustering_col` behaves exactly as before and produces unsuffixed names, so existing pipelines are unaffected; a single **named** one takes the suffixed path, which is what a caller generating one scope at a time wants.

- **Breaking:** `wgcna_<group>_powertest` now returns a small `list` instead of the Seurat object. `TestSoftPowers()` only stashes a scale-free fit table on a copy of the whole object, so storing the object cost several GB per cell group (~3.5 GB each, ~70 GB across the 20 groups of the pipeline this was found in) to keep a few hundred numbers. The list holds `power_table` (the table `hdWGCNA::PlotSoftPowers()` draws, so a report can redraw it from this target alone), `soft_power`, `reached_threshold`, `sft_rsquared` and `group`. `reached_threshold` is new information: the power pick silently falls back to the highest tested power when nothing reaches the `SFT.R.sq` threshold, which is a materially different claim about the network, and the table alone did not say which had happened. `wgcna_<group>_soft_power` now reads `[["soft_power"]]` from that list, and the network step rebuilds `datExpr` from `wgcna_prep` with `SetDatExpr()` rather than inheriting it, which is cheap next to `ConstructNetwork()`. Code that read the power-test target as a Seurat object must be updated; nothing inside this package did.

- New `delete_tom` argument to `tar_hdwgcna()` (default `TRUE`), deleting the topological overlap matrix as soon as `ConstructNetwork()` has finished with it. The TOM is a dense gene-by-gene matrix, around 3 GB per cell group, and `ConstructNetwork()` is the only thing that reads it: the modules, eigengenes, connectivity and every downstream target are derived inside that call, so the file is write-only for a pipeline that does not plot networks. Both the final `<tom_name>_TOM.rda` and any `<tom_name>_block.N.rda` left behind by a multi-block run are removed, by exact path rather than by pattern so that group names containing spaces or regex metacharacters are handled, and the step reports how much it freed. Deletion happens immediately after `ConstructNetwork()` rather than at the end of the target, so the space is reclaimed even if a later step of the same target fails. Set `delete_tom = FALSE` when you intend to call the hdWGCNA functions that do re-read the TOM, namely `HubGeneNetworkPlot()`, `ModuleNetworkPlot()` and `RunModuleUMAP()` / `ModuleUMAPPlot()`.

- `tar_hdwgcna()` strips a trailing slash from `tom_outdir` before appending a scope subdirectory, so the per-scope paths read `./data/hdwgcna/0.3` rather than `./data/hdwgcna//0.3`.

# scitargets 1.5.0

- New `tar_hdwgcna()` targets factory for an hdWGCNA co-expression analysis (following the hdWGCNA basic + differential module-eigengene tutorials). A single shared metacell-preparation step (`wgcna_prep`) is created once for all groups, and the per-group steps `wgcna_<group>_powertest`, `wgcna_<group>_soft_power`, `wgcna_<group>` and `wgcna_<group>_dmes` are created for each `group` (target names carry a lower-cased, sanitised group label). All values passed to the hdWGCNA / WGCNA functions are exposed as arguments; `input_obj`, `clustering_col` and `patient_col` are required and validated. The soft power is selected as `hdWGCNA::PlotSoftPowers()` highlights it (lowest tested power with `SFT.R.sq >= sft_rsquared`); modules are associated with a clinical trait via `ModuleTraitCorrelation()` (Spearman by default, on one-hot trait indicators) and pairwise `FindDMEs()`. `deployment = "main"` for the prep, power test, soft-power, network and Enrichr steps (the power test and network load the whole shared Seurat object, so running several groups concurrently on workers would hold multiple full copies in memory); only the DME step runs on a worker. `WGCNA` and `hdWGCNA` were added to Suggests.

- New hdWGCNA report helpers, mirroring `dea_report_lines()`. `hdwgcna_report_lines(group, ...)` builds a full Quarto tabset for one cell group's hdWGCNA results (soft-power, dendrogram + modules, module-trait correlation, Enrichr GO enrichment, module-eigengene dot plot, differential module eigengenes); emit it with `knitr::knit_child()`. Its chunks read the `tar_hdwgcna()` targets directly via `tar_read()`. Supporting plot/table helpers are also exported: `hdwgcna_module_trait_heatmap()`, `hdwgcna_module_trait_table()`, `hdwgcna_associated_modules()`, `hdwgcna_enrichr_barplot()`, `hdwgcna_me_dotplot()` and `hdwgcna_dme_heatmap()`.

- `tar_demultiplex_hto()` and `tar_hdwgcna()` gain `deployment` (`"main"`/`"worker"`) and `controller` arguments to choose where their steps run, including routing to a named controller of a `crew::crew_controller_group()` (via `targets::tar_resources()`). Steps that must stay on the main process are **not** affected: the raw-load step of `tar_demultiplex_hto()` (Seurat metadata creation bugs on crew workers) and every WGCNA-multithreaded / whole-object step of `tar_hdwgcna()` (`wgcna_prep`, power test, soft-power, network, Enrichr) — only `tar_hdwgcna()`'s `dmes` step is relocatable. Both default to `deployment = "main"`, preserving previous behaviour.

# scitargets 1.4.1

- Better handling of numeric covariates for pseudobulk DESeq2 test (previously forcibly converted to character).

- Added parameters to run_dea to allow the user to set `Seurat::FindMarkers` test.use, group.by, only.pos, min.pct, and logfc.threshold inputs.

- Added a message call warning agains interpretating pseudobulk GO when using LRT test. GO is still computed but the default Wald test should be used for the enrichment computations.

- Better handling of GO p-values enrichment. Previously, the adjustment was only realized in the "top N" terms returned by topGO. Now the adjustment is done on all tested terms.

- GSEA is now only computed where the signed gene ranking is valid. For pseudobulk it runs only with `pb_test = "Wald"` and ranks genes by the signed DESeq2 Wald statistic (the `stat` column); it is skipped for `"LRT"` (whose statistic is unsigned). For single cell it runs only when `Seurat::FindMarkers` keeps the full gene universe (`sc_logfc_threshold = 0`, `sc_min_pct = 0`, `sc_only_pos = FALSE`) and keeps ranking genes by `sign(log2FC) * -log10(p_val)`. When GSEA is skipped a warning is emitted and the report omits the GSEA tab for that level.

- GSEA ranking now resolves duplicated gene symbols by keeping the entry with the strongest absolute statistic (previously the arbitrary first occurrence), drops non-finite ranks, and errors early if the marker table has no `p_val` column.

- Pseudobulk DE now checks, before aggregation, that each `pseudobulk_unit` belongs to a single comparison group and errors listing the offending units otherwise (the DESeq2 design assumes one group per replicate).

- Pseudobulk `DESeq2::results()` is now called with `alpha` set to the per-level pseudobulk adjusted-p cutoff (from the object's `padj_cutoffs`, default 0.05) instead of DESeq2's default 0.1, so independent filtering is optimised for the FDR threshold actually used downstream.

- Hardened the pseudobulk DESeq2 input: the `colData` rows are aligned and named to the count-matrix columns before `DESeqDataSetFromMatrix()`; only the comparison samples (group1 vs the reference) are kept, rather than relying on Seurat to drop non-comparison cells; and design covariates are joined to pseudobulk samples by `pseudobulk_unit` (one value per replicate), erroring if a covariate is not constant within a unit.

- Fixed `gsea_barplot()` top-N selection: gene sets are now chosen by significance (adjusted p, then p-value) rather than by the display metric, so the `"signed_pval"` metric no longer shows the least significant sets.

- Removed the `pb_covariate_key` argument from `run_dea()`. Design covariates (referenced by `pb_design`) are now always joined to pseudobulk samples by `pseudobulk_unit` — one value per biological replicate — which is the only correct key for a DESeq2 replicate-level design.

- The xlsx export (`write_dea_xlsx` / `dea_write_xlsx`) now uses `openxlsx2` instead of the unmaintained `openxlsx`. The "summary" sheet's "Comparison summary" block renames `levels` to `computed_levels` (all levels computed for the comparison) and adds `level_in_document` (the level written in this file); the `reason` row is shown only when non-empty. Per-level files no longer prefix sheet names with `sc_`/`pb_` (the prefix is kept only when several levels share one workbook).

- Report generation: new `dea_comparisons_lines()` groups every comparison's `dea_report_lines()` block under a single **"DEA Comparisons"** tab (each comparison a nested child tab), as a sibling of the cell-state-composition / testability tabs, instead of looping `dea_report_lines()` by hand. The comparison "Comparison summary" tab now renders its cell-count tables with `DT2`. `composition_plot()` drops its (large) legend and scales the canvas width with the number of samples (the interactive `ggiraph` tooltip carries the cell state); `composition_boxplot()` lays its facets out in two columns and scales height with the number of cell states. `gsea_barplot()`'s subtitle is more explicit ("gene sets with adjusted p-value < ...").

- GO and GSEA result tables are now stored in full inside the `scitargets_dea` object (the GO table was previously truncated to `top_nodes` at compute time). `go_table()` and `gsea_table()` gain an `n_terms` argument: it caps the returned rows to the object's `top_nodes` by default (used by the report and dashboards) or returns the entire table with `n_terms = Inf`. The xlsx export (`dea_write_xlsx`) now writes the complete tables.

- GO bar/network plots and the genes-terms list show the top-N terms by adjusted p-value (not only the cutoff-passing ones), so a plot is produced whenever terms exist. `go_barplot()` draws a dashed reference line at the adjusted-p cutoff. GO and GSEA bar plots drop terms / gene sets with an adjusted p-value of 1 (their `-log10` bar has zero height).

- The report's GO and GSEA tables (`dea_report_lines()`) drop the per-term `Genes` (GO) and `leadingEdge` (GSEA) columns for readability; those columns remain in the object and in the xlsx export.


# scitargets 1.4.0

For this release `Claude Code (Opus 4.8)` was used to speed up development. The code was manually checked and tested in a work project.

## Seurat 5.5.0 compatibility

- Changes made to adjust to latest `Seurat` modification (version 5.5.0):

  1.  inside `Seurat::CreateSeuratObject` made it unreliable on workers thread. Considering the low requirement in term of time for `scitargets::load_seurat_data_10X` the `scitargets::tar_demultiplex_hto` `seurat_obj_<RUN_ID>_raw` step was set to always be deployed on the `main` R process.
  2.  Updated some filter_cell_and_run_reduction to determine the number of PC to use when computing UMAP and tSNE using `SeuratObject::Stdev(obj, reduction="pca")` as `Seurat::ElbowPlot` return a ggplot2 obj in the latest versions of Seurat.
  3.  Seurat requirment was raised to the version `5.5.0`.

## New `scitargets_dea` differential-expression subsystem

- Added the `scitargets_dea` S7 class and its analysis functions (moved in from the PRETERRAH project). It holds, for one clinical-group comparison within one cluster, the differential-expression results (single-cell `Seurat::FindMarkers` and/or pseudobulk `DESeq2`), GO enrichment (one or more of BP/CC/MF, adjusted per ontology), GSEA results (MSigDB Hallmark / <GO:BP> / C7 ImmuneSigDB via `fgsea`), the pseudobulk PCA / outlier coordinates and the per-level p-value cutoffs.

  - New S7 generics: `markers_table()`, `go_table()`, `go_plot_data()`, `go_genes_html()`, `volcano_plot()`, `pca_plot()`, `go_barplot()`, `go_cnetplot()`, `gsea_table()`, `gsea_barplot()`, `dea_write_xlsx()`. S7 methods are registered in `.onLoad()`.

  - New functions: `run_dea()` (single-entry dispatcher), `dea_comparisons()` (enumerate comparisons, with `groups` and `clusters` selectors), `get_msigdbr_pathways()`, `write_dea_xlsx()`, `normalize_dea_list()`, `is_scitargets_dea()`.

- `run_dea()` analysis controls (pseudobulk / `DESeq2` unless noted):

  - **DESeq2 model:** `pb_test` (`"Wald"` / `"LRT"`), `pb_design` and `pb_reduced` design formulas (character, passed via `as.formula`) and `pb_low_count_filter = c(min_count, min_samples)` low-count pre-filter. Covariates referenced by `pb_design` are looked up one value per biological replicate (keyed on `pseudobulk_unit`). Per-group mean normalized counts are added to the result table.
  - **Replicate gate:** `min_replicates` (default 3) — minimum biological replicates per group for a pseudobulk comparison to run.
  - **LFC shrinkage:** `pb_lfc_shrink` toggle; an `apeglm` -> `ashr` -> `normal` cascade adds a shrunken log2 fold-change, and a message reports which method was used.
  - **Pseudobulk PCA + outlier detection:** `pb_pca`, `pca_n_top_genes`, `pca_outlier_conf` — `DESeq2` VST -> PCA -> Mahalanobis outlier flagging (`mt::pca.outlier`). `pb_remove_outliers` drops the flagged replicates from the DESeq2 fit (re-checking the replicate gate afterwards).
  - **Per-level adjusted-p cutoffs:** `padj_cutoff_single_cell` / `padj_cutoff_pseudobulk`, stored on the object; each drives that level's volcano cutoff line and GO foreground-gene selection.
  - **Species:** `species` (`"human"` / `"mouse"`) is the single source of truth for the GO annotation DB (`org.Hs.eg.db` / `org.Mm.eg.db`) and the GSEA `msigdbr` species.

### Interactive plots

- `volcano_plot()`, `pca_plot()`, `go_barplot()` and `gsea_barplot()` gained an `interactive` argument. The default is a static `ggplot2` figure (the static volcano labels significant genes via `ggrepel`); when interactive, each is rendered with `ggiraph` and per-element tooltips.

- `gsea_barplot()` now defaults to displaying the **FDR-adjusted** enrichment significance of each gene set (`metric = "signed_nlog10_padj"`, i.e. `sign(NES) * -log10(adjusted p)`) and selects pathways by `padj < cutoff`. The raw-p (`"signed_nlog10_pval"`) and signed-p metrics remain available via `metric=`. The gene-level ranking fed to GSEA is unchanged (raw per-gene p-value).

### Quarto report generation

- `dea_report_lines()` builds the full Quarto/knitr child block for one comparison (comparison summary, volcano, DE table, GO and GSEA tabsets). It gained: `levels` (split the report into single-cell / pseudobulk documents), `interactive_plots` (per-plot `ggiraph` opt-in), patient cell-state composition options, the `clusters_to_show` / `groups_to_show` display filters (+ `available_*` for validation), `gsea_metric` (defaults to the FDR display), a top-DE-gene heatmap (named-list of Quarto chunk options + a `Seurat` object), and `DT2` tables with captions for the DE and GSEA tables.

- New cross-comparison companion generators (each rendered once per document):

  - `composition_plot()` / `composition_boxplot()` / `composition_lines()` — patient cell-state composition: a stacked per-patient state-proportion barplot and a per-state proportion boxplot with jittered per-patient points, faceted by clinical group.
  - `dea_upset_plot()` / `dea_upset_lines()` — UpSet plots (via `ggupset`) of DE-gene overlap across comparisons, for down-regulated / up-regulated / all DE genes.
  - `dea_testability_heatmap()` / `dea_testability_lines()` — testability heatmap (cell state x comparison) of the minimum replicates / cells per group and whether each comparison is testable, shown before and after outlier removal (plus an outliers-removed count heatmap).
  - `top_de_heatmap()` / `dea_top_de_lines()` — top DE-gene heatmap across all cells (`Seurat::DoHeatmap`), using the top genes by |fold-change| of the shown comparisons.
  - `dea_gsea_yield_heatmap()` / `dea_gsea_yield_lines()` — GSEA-yield heatmap: the number of significantly enriched gene sets (`padj < cutoff`) per cell state x comparison, one heatmap per MSigDB collection.

## Dependencies

- The DEA subsystem adds `mt` (PCA outlier detection) and `DT2` (interactive tables) to `Imports`, and `ggupset`, `apeglm` and `ashr` to `Suggests`. `org.Hs.eg.db` and `org.Mm.eg.db` are in `Suggests` and attached on demand (via `require()`) when the matching species is used.

## Other

- Modified the shape of `pattern` to hexagon for better visual of the graph in `scitargets::tar_visnetwork_enhanced`.

# scitargets 1.3.0

- Added a new `plot_resolution_tree` function that create a direcred tree graph representing how clusters split across multiple clustering resolutions.

- Modified the assay used by Azimuth to the RNA assay instead of the SCT assay as Azimuth normalize the data.

# scitargets 1.2.0

- Modified the QMD templates. Now the quality control violins plot display lines representing the min/max features and mitochondrial cutoff set on tar_demultiplex_hto.

- Correction of tar_demultiplex_hto that previously used run_id instead of run_path for a path check.

# scitargets 1.1.0

- Modified the extract_singlets function to add a parameter to indicate a numeric vector representing the clustering resolutions the pipeline compute (default to 0.2 to 1.5 with a step of 0.1).

- Modification of the tar_demultiplex_hto factory to use the new parameter of extract_singlets.

- Modification of tar_visnetwork_enhanced function to allow user to filter the report and specify which dataset to plot.

- Modification of the report to use the new version of tar_visnetwork_enhanced function.

# scitargets 1.0.1

- QMD template bugfix

- tar_demultiplex_hto bugfix: previous use of "run_azimuth = T" was bugged when removing features

# scitargets 1.0.0

Initial release
