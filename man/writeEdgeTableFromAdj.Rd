\name{writeEdgeTableFromAdj}
\alias{writeEdgeTableFromAdj}
\title{Write a network edge CSV file from an adjacency matrix.}
\usage{
  writeEdgeTableFromAdj(adj, filename, noEdge = c(0, NA),
    diag = FALSE, sep = "\t")
}
\arguments{
  \item{adj}{Adjacency matrix.}

  \item{filename}{filename to write to}

  \item{noEdge}{values in the adjacency matrix to consider
  as "no edge".  Defaults to 0 and NA.}

  \item{diag}{Logical; Include the diagonals?}

  \item{sep}{Defaults to tab. See details.}
}
\value{
  NULL
}
\description{
  The Cytoscape network format expects a csv file, where
  each row contains an edge pair, 'source' and 'target',
  along with an edge weight. This function will convert an
  adjacency matrix to this 'long' format.
}
\details{
  In most cases a tab is the best separator between the
  elements in each row, but if for some reason the node
  names have tabs in them, this can be changed.
}

