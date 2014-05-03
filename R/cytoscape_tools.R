#' @name networkTools-package
#' @aliases networkTools-package networkTools
#' @docType package
#' 
#' @title Tools for analysing and visualising networks.
#' @description
#' This R package provides functions for analysisng and visualising large 
#' networks, primarily those generated from the R package WGCNA.
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
<<<<<<< HEAD

`%nin%` <- function(...) !(`%in%`(...))

=======
>>>>>>> 619addf57a7b0a65e546d08c0832bcd1a6933d21
