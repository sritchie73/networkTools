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

`%nin%` <- function(...) !(`%in%`(...))

#' Converts the output of \code{NEO} into something sane.
#' 
#' \code{NEO} outputs a massive CSV file that is not easily parseable, nor
#' loadable into cytoscape. This function will parse this CSV file and write
#' two tables to file: 
#' 
#' An edge table, filled with relevant columns, notably
#' splitting out the \code{edge} column into Source and Target nodes, and 
#' splitting the \code{Final.SNPs.LEO.NB.OCA} column into a vector of SNPs for 
#' the target and source nodes respectively.
#' 
#' A second table, filled with SNP to Node relationships will also be output.
#' 
#' @param neofilename path to the neo csv file.
#' @param outfile file prefix to output the .edge.table.csv and .snp.table.csv
#'        files
#' @export
parseNEOcsv <- function(neofilename, outfile) {
  neo.csv <- read.csv(neofilename, stringsAsFactors=FALSE)
  snps <- strsplit(neo.csv$Final.SNPs.LEO.NB.OCA, ";")
  df <- data.frame(
    stringsAsFactors=FALSE,
    Source = sapply(strsplit(neo.csv$edge, " -> "), `[`, 1),
    Target = sapply(strsplit(neo.csv$edge, " -> "), `[`, 2),
    CPA.Anchors = sapply(snps, function(x) gsub("fsnp.cpa=", "", x[1])),
    OCA.Anchors = sapply(snps, function(x) gsub("fsnp.oca=", "", x[2])),
    LEO.NB.OCA = neo.csv$LEO.NB.OCA,
    LEO.NB.CPA = neo.csv$LEO.NB.CPA,
    P.Model.AB = neo.csv$Model.P.value.AtoB,
    P.weighted.LEO.NB.OCA = neo.csv$P.weighted.LEO.NB.OCA,
    Z.WaldPath.AB = neo.csv$ZtestPathCoefficientAB,
    P.WaldPath.AB = pnorm(neo.csv$ZtestPathCoefficientAB, lower.tail=FALSE),
    PearsonCor = neo.csv$PearsonCor
  )
  
  write.csv(df[!is.na(CPA.Anchors)], file=paste0(outfile, ".edge.table.csv"), 
            quote=FALSE)
  write.csv(df[is.na(CPA.Anchros), c("Source", "Target", "PearsonCor")],
            file=paste0(outfile, ".snp.table.csv"), quote=FALSE)
}

