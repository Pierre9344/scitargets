# Register S7 methods when the package is loaded.
# See https://rconsortium.github.io/S7/articles/packages.html
.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}

# Column names referenced via non-standard evaluation in dplyr/ggplot2 pipelines
# (.build_cnet_data, .compute_topgo_result; the composition / testability / UpSet /
# GSEA-yield plot aes()); declare them so R CMD check does not flag
# "no visible binding for global variable".
utils::globalVariables(c(
  "GO.ID", "Term", "GO_label", "adjpval", "Significant", "Genes",
  "xpos", "ypos", "Enrichment",
  "comparison", "value", "tooltip", "count",
  "group", "proportion", "patient", "state"
))
