# scitargets 1.1.0

Modified the extract_singlets function to add a parameter to indicate a numeric vector representing the clustering resolutions the pipeline compute (default to 0.2 to 1.5 with a step of 0.1).

Modification fo the tar_demultiplex_hto factory to use the new parameter of extract_singlets.

# scitargets 1.0.1

QMD template bugfix

tar_demultiplex_hto bugfix: previous use of "run_azimuth = T" was bugged when removing features

# scitargets 1.0.0

Initial release
