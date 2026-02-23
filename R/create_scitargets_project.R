create_quarto_yaml <- function(path) {
  project_name <- base::strsplit(path, "/") %>%
    base::unlist() %>%
    utils::tail(n = 1)
  template_path <- system.file("rstudio/templates/project/files", package = "scitargets")
  file.copy(
    from = paste0(template_path, c("/_quarto.yml", "/_targets.R", "/style.css", "/QMD")),
    to = path,
    recursive = T,
    overwrite = T
    )
  file.copy(
    from = paste0(template_path, "/rprofile.txt"),
    to = paste0(path, "/.Rprofile"),
    overwrite = T
  )
}


create_scitargets_project <- function(path, overwrite = F, ...) {
  if (!overwrite && dir.exists(path)) {
    stop("Directory already exists.")
  }
  dir.create(paste0(path, "/data/cellranger_output"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(path, "/R"), recursive = TRUE, showWarnings = FALSE)
  create_quarto_yaml(path = path)
}
