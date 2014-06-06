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

#' Retrieve all nodes in the Neighborhood of a given node
#' 
#' This function mimics the behaviour of \code{neighborhood} in the 
#' \code{igraph} package, but is backwards compatible with \code{R 2.5.1}
#' 
#' @details
#' The neighborhood of a given order ‘o’ of a vertex ‘v’ includes all
#' vertices which are closer to ‘v’ than the order. Ie. order 0 is
#' always ‘v’ itself, order 1 is ‘v’ plus its immediate neighbors,
#' order 2 is order 1 plus the immediate neighbors of the vertices in
#' order 1, etc.
#' 
#' @param node Character constant, Node you want to retrieve neighborhood of.
#' @param adjacency an adjacency matrix
#' @param order Integer giving the order of the neighborhood.
#' @param mode Character constant, it specifies how to use the direction of
#'            the edges if a directed graph is analyzed. For ‘out’ only the
#'            outgoing edges are followed, so all vertices reachable from
#'            the source vertex in at most ‘order’ steps are counted. For
#'            ‘"in"’ all vertices from which the source vertex is reachable
#'            in at most ‘order’ steps are counted. ‘"all"’ ignores the
#'            direction of the edges.
#' @return a vector of node names in the requested neighborhood.
#' @export
neighborhood <- function(node, adjacency, order, mode="all") {
  stopifnot(order >= 0)
  stopifnot(mode %in% c("in", "out", "all"))
  stopifnot(class(adjacency) == "matrix")
  stopifnot(class(node) == "character")
  if (order == 0) {
    node
  } else {
    if (mode == "in") {
      neighbors <- names(adjacency[,node][adjacency[,node] == 1])
    } else if (mode == "out") {
      neighbors <- names(adjacency[node,][adjacency[node,] == 1])
    } else {
      neighbors <- c(
        names(adjacency[,node][adjacency[,node] == 1]),
        names(adjacency[node,][adjacency[node,] == 1])
      )
    }
    unique(c(neighbors, 
             unlist(sapply(neighbors, neighborhood, adjacency, order-1, mode)),
             node
           ))
  }
}
