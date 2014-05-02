\name{PlotNetworkHeatmap}
\alias{PlotNetworkHeatmap}
\title{Plot network edge weights as a heatmap.

Plots the edge weights between nodes as a heatmap, optionally drawing boxes
around pre-determined clusters. The \code{xlab} and \code{ylab} arguments are
useful when comparing the edge weights}
\usage{
PlotNetworkHeatmap(edge.matrix, network.labels = NULL, xlab = "",
  ylab = "")
}
\arguments{
  \item{edge.matrix}{a square matrix containing the edge
  weights between each pair of nodes.}

  \item{network.labels}{an optional named vector assigning
  each node to a network/cluster.}

  \item{xlab}{}

  \item{ylab}{}
}
\description{
Plot network edge weights as a heatmap.

Plots the edge weights between nodes as a heatmap,
optionally drawing boxes around pre-determined clusters.
The \code{xlab} and \code{ylab} arguments are useful when
comparing the edge weights
}
\note{
It expects nodes to be pre-ordered within the edge weight
matrix: I.e. clusters should be in continuous blocks.
Otherwise, garbage-in, garbage-out.
}

