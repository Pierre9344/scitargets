# scitargets 0.3.0

-   Correction of a bug in the `demultiplex_cell` function
-   Correction of bugs concerning the removal of features in the `tar_demultiplex_hto` function
-   Modification of the quarto template document "02_run_1.qmd".
-   Added a run_path input to `tar_demultiplex_hto` allowing to indicatesubdirectory when indicating the runs inside the data/cellranger_output/ folder.

# scitargets 0.2.1

-   Modification of README.md
-   Correction of vignettes
-   Correction of autometric setup in _targets.R by removing the deprecated pids call from crew controller (code is commented by default)



# scitargets 0.2.0

-   Initial github commit
-   Added support for realize an automatic annotation of cell types using azimuth and for the identification of markers that differentiate cell clusters. Both can be activated by indicating a valid metadata field of the singlets to the new "singlets_clusters_to_use" argument of "tar_demultiplex_hto".
-   Modification of the template:
    -   The template now use "quarto::write_yaml_metadata_block" and "targets::tar_exist_objects" to check if the users removed singlets in the pipeline. If this is the case, it will show PCA and UMAP of this object.
    -   The same functions are used to check if there are objects with name based on "seurat_obj_RUN_ID_azimuth" or "markers_RUN_ID" and show the results when needed.
-   The template now need the quarto R package version to be 1.5.1 or higher.

# scitargets 0.1.0

-   Initial version.
