#' Calculate the degree of each node in the network
#' 
#' @param adjacency an adjacency matrix
#' @return a named vector
#' @export
degree <- function(adjacency) {
  sort(colSums(adjacency) + rowSums(adjacency), decreasing=TRUE)
}

#' Identify sink nodes in the network
#' 
#' @param adjacency an adjacency matrix
#' @return a vector of node names
#' @export
sinks <- function(adjacency) {
  d <- degree(adjacency)
  names(d)[(d - colSums(adjacency)[names(d)]) == 0]
}

#' Identify source nodes in the network
#' 
#' @param adjacency an adjacency matrix
#' @return a vector of node names
#' @export
sources <- function(adjacency) {
  d <- degree(adjacency)
  names(d)[(d - rowSums(adjacency)[names(d)]) == 0]
}

#' Identify mediator nodes in the network
#' 
#' @param adjacency an adjacency matrix
#' @return a vector of node names
#' @export
mediators <- function(adjacency) {
  sinks <- sinks(adjacency)
  sources <- sources(adjacency)
  colnames(adjacency)[!(colnames(adjacency) %in% c(sinks, sources))]
}
