# scitargets 0.2.0

-   Initial github commit
-   Added support for realize an automatic annotation of cell types using azimuth and for the identification of markers that differentiate cell clusters. Both can be activated by indicating a valid metadata field of the singlets to the new "singlets_clusters_to_use" argument of "tar_demultiplex_hto".
-   Modification of the template:
    -   The template now use "quarto::write_yaml_metadata_block" and "targets::tar_exist_objects" to check if the users removed singlets in the pipeline. If this is the case, it will show PCA and UMAP of this object.
    -   The same functions are used to check if there are objects with name based on "seurat_obj_RUN_ID_azimuth" or "markers_RUN_ID" and show the results when needed.
-   The template now need the quarto R package version to be 1.5.1 or higher.

# scitargets 0.1.0

-   Initial version.
