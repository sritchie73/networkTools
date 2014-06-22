#' @name networkTools-package
#' @aliases networkTools-package networkTools
#' @docType package
#' 
#' @title Tools for analysing and visualising networks.
#' @description
#' This R package provides functions for analysisng and visualising large 
#' networks, primarily those generated from the Horvath et al suite of network
#' packages and tools, e.g. WGCNA and the NEO software 
#' \link{http://labs.genetics.ucla.edu/horvath/aten/NEO/}.
#' 
#' This package is backwards compatible with older versions of R (i.e. 2.5.1).
#' 
#' @details
#' \tabular{ll}{
#'  Package: \tab networkTools\cr
#'  Type: \tab Package\cr
#'  Version: \tab 1.0\cr
#'  Date: \tab 2014-02-19\cr
#'  License: \tab GPL (>= 2)\cr
#' }
#' 
#' @author 
#' Scott Ritchie \email{sritchie73@@gmail.com}
#' 
#' @rdname networkTools-package
#' 
NULL

#' Convert an adjacency matrix to an edge list
#'   
#' @param adj Adjacency matrix.
#' @param noEdge values in the adjacency matrix to consider as "no edge". 
#'  Defaults to 0 and NA.
#' @param diag Logical; Include the diagonals?
#' @return NULL
#' 
#' @export
adj2edge <- function(adj, noEdge=c(0, NA), diag=FALSE) {
  nrow <- sum(!(adj %in% noEdge))
  if (!diag) {
    nrow <- nrow - sum(!(diag(adj) %in% noEdge))
  }
  edgeList <- data.frame(Source=rep("", nrow), Target="", weight=0)
  e <- 1
  for (i in 1:nrow(adj)) {
    for (j in 1:ncol(adj)) {
      if ((i != j || diag) && (adj[i, j] %nin% noEdge)) {
        rname <- ifelse(is.null(rownames(adj)), i, rownames(adj)[i])
        cname <- ifelse(is.null(colnames(adj)), j, colnames(adj)[j])
        edgeList[e,] <- list(rname, cname, adj[i,j])
        e <- e + 1
      }
    }
  }
  edgeList
}


#' Converts an edge table to an adjacency matrix
#' 
#' @note
#' Does work on networks with multiple edges between nodes.
#' 
#' @param edgetable assume the first two columns are the edge from the node 
#'   in column 1 to the node in column 2.
#' @param weights.col if specified, the adjacency matrix is filled with the 
#'   weights from the given column number/name.
#' @return
#'   an adjacency matrix
#'   
#' @export
edge2adj <- function(edgetable, weights.col=NULL) {
  nodes <- unique(c(edgetable[,1], edgetable[,2]))
  adj <- matrix(0, length(nodes), length(nodes), dimnames=list(nodes, nodes))
  for (i in 1:nrow(edgetable)) {
    if (!is.null(weights.col)) {
      adj[edgetable[i,1], edgetable[i,2]] = adj[edgetable[i,1], edgetable[i,2]] + edgetable[i,weights.col]
    } else {
      adj[edgetable[i,1], edgetable[i,2]] = adj[edgetable[i,1], edgetable[i,2]] + 1
    }
  }
  adj
}

`%nin%` <- function(...) !(`%in%`(...))

