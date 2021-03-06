% Generated by roxygen2 (4.0.0): do not edit by hand
\name{edge2adj}
\alias{edge2adj}
\title{Converts an edge table to an adjacency matrix}
\usage{
edge2adj(edgetable, weights.col = NULL)
}
\arguments{
\item{edgetable}{assume the first two columns are the edge from the node
in column 1 to the node in column 2.}

\item{weights.col}{if specified, the adjacency matrix is filled with the
weights from the given column number/name.}
}
\value{
an adjacency matrix
}
\description{
Converts an edge table to an adjacency matrix
}
\note{
Does work on networks with multiple edges between nodes.
}

