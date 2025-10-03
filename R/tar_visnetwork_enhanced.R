utils::globalVariables(c("type", "from", "to", "description", "has_function",
                         "has_stem", "function_parent"))


#' Enhanced visNetwork dependency graph
#'
#' Enhanced version of the `targets::tar_visnetwork` function.
#'
#' Visualize the pipeline dependency graph with a visNetwork HTML widget.
#'
#' @param color_stem Color for stems (targets steps) that depends of other steps
#' @param color_standalone Color for stems that don't depends of other steps (input files)
#' @param color_function Color for the functions used by the pipeline
#' @param color_pattern Color for the steps that use dynamic/static branching
#' @param redirect_stem_to_function If set to true and a stem depends of both a function and another stem, then we redirect the edge to the function.
#'
#' @returns A HTML widget with a visNetwork graph representing the project targets pipeline.
#' @export
#'
#' @examples
#' \dontrun{
#' tar_visnetwork_enhanced()
#' }
tar_visnetwork_enhanced <- function(
    color_stem = "#75AADB",       # Light Blue for stems that depends of other
    color_standalone = "#F4A261", # Orange for standalone stems
    color_function = "#9CCC65",   # Green for functions
    color_pattern = "#E5989B",    # Soft pink for pattern nodes
    redirect_stem_to_function = TRUE
) {
  # Generate the network structure
  network <- targets::tar_network()

  # Include stem, function, and pattern nodes (ignore 'object' nodes)
  nodes <- network$vertices %>%
    dplyr::filter(type %in% c("stem", "function", "pattern"))

  # Prepare nodes
  nodes$id <- nodes$name
  nodes$label <- nodes$name

  #tooltip_column <- if ("command" %in% colnames(nodes)) "command" else "type"
  nodes %<>%
    dplyr::mutate(title = base::ifelse(base::is.na(description), type, base::paste0(type, ": ", description))
    )

  # Identify dependent nodes
  dependent_nodes <- base::unique(network$edges$to)

  # Assign colors
  nodes$color <- dplyr::case_when(
    nodes$type == "stem" & !(nodes$name %in% dependent_nodes) ~ color_standalone,
    nodes$type == "stem" ~ color_stem,
    nodes$type == "function" ~ color_function,
    nodes$type == "pattern" ~ color_pattern
  )

  # Assign shapes
  nodes$shape <- dplyr::case_when(
    nodes$type == "function" ~ "box",
    nodes$type == "pattern" ~ "database",  # Pattern nodes can be visually distinct
    TRUE ~ "dot"
  )

  # Filter edges to only include connections among relevant nodes
  edges <- network$edges %>%
    dplyr::filter(from %in% nodes$name & to %in% nodes$name)

  # Redirect logic: stem → function if a stem depends on both stem & function
  if (redirect_stem_to_function) {
    stem_dependencies <- edges %>%
      dplyr::inner_join(nodes, by = c("to" = "name")) %>%
      dplyr::group_by(to) %>%
      dplyr::summarize(
        has_function = base::any(from %in% nodes$name[nodes$type == "function"]),
        has_stem = base::any(from %in% nodes$name[nodes$type == "stem"]),
        function_parent = dplyr::first(from[from %in% nodes$name[nodes$type == "function"]]),
        .groups = "drop"
      ) %>%
      dplyr::filter(has_function & has_stem) %>%
      dplyr::select(to, function_parent)

    edges <- edges %>%
      dplyr::left_join(stem_dependencies, by = "to") %>%
      dplyr::mutate(to = ifelse(!is.na(function_parent) & from %in% nodes$name[nodes$type == "stem"], function_parent, to)) %>%
      dplyr::select(from, to)
  }

  # Style edges
  edges$arrows <- "to"
  edges$color <- "black"

  # Render the network
  visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visInteraction(hover = TRUE)
}
