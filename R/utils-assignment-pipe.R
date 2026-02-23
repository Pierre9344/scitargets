#' Assignment pipe operator
#'
#' See \code{magrittr::\link[magrittr:compound]{\%<>\%}} for details.
#'
#' @name %<>%
#' @rdname compound
#' @aliases %<>%
#' @export
#' @family exported functions
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
#' @param lhs An object which serves both as the initial value and as target.
#' @param rhs A function call using the magrittr semantics.
#' @details The assignment pipe, %<>%, is used to update a value by first piping it into one or more rhs expressions, and then assigning the result. For example, some_object %<>% foo %>% bar is equivalent to some_object <- some_object %>% foo %>% bar. It must be the first pipe-operator in a chain, but otherwise it works like %>%.
NULL
