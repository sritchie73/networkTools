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

#' Write a network edge CSV file from an adjacency matrix.
#'
#' The Cytoscape network format expects a csv file, where each row contains an
#' edge pair, 'source' and 'target', along with an edge weight. This function
#' will convert an adjacency matrix to this 'long' format.
#' 
#' @details 
#'  In most cases a tab is the best separator between the elements in each row,
#'  but if for some reason the node names have tabs in them, this can be 
#'  changed.
#'   
#' @param adj Adjacency matrix.
#' @param filename filename to write to
#' @param noEdge values in the adjacency matrix to consider as "no edge". 
#'  Defaults to 0 and NA.
#' @param diag Logical; Include the diagonals?
#' @param sep Defaults to tab. See details.
#' @return NULL
#' 
#' @export
writeEdgeTableFromAdj <- function(adj, filename, noEdge=c(0, NA), diag=FALSE,
                                  sep="\t") {
  sink(filename)
  for (i in 1:nrow(adj)) {
    for (j in 1:ncol(adj)) {
      if ((i != j || diag) && (adj[i, j] %nin% noEdge)) {
        rname <- ifelse(is.null(rownames(adj)), i, rownames(adj)[i])
        cname <- ifelse(is.null(colnames(adj)), j, colnames(adj)[j])
        cat(paste(rname, cname, adj[i,j], sep=sep), "\n")
      }
    }
  }
  sink()
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

