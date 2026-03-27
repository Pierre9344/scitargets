# scitargets 1.3.0

- Added a new `plot_resolution_tree` function that create a direcred tree graph representing how clusters split across multiple clustering resolutions.

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
