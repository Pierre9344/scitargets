---

editor_options: 
  markdown: 
    wrap: 72
---

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
