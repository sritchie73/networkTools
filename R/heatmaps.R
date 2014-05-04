#' Order nodes and clusters in a network
#' 
#' Nodes are ordered by degree, and clusters are ordered by similarity if some
#' underlying summary profile is provided in \code{cluster.summaries}
#' 
#' @param edge.matrix Matrix of edge weights 
#' @param node.labels Optional vector of cluster assignments for each node.
#' @param cluster.summaries Optional data.frame of summary measurements to 
#'                          compare clusters with.
#' @return a vector of ordered nodes. 
#'
OrderNetwork <- function(edge.matrix, node.labels=NULL, cluster.summaries=NULL) {
  # Cluster gene networks by Eigengene similarity
  order <- NULL
  
  if (!is.null(node.labels)) {
    if (!is.null(cluster.summaries)) {
      h <- hclust(as.dist(1 - abs(cor(cluster.summaries))))
      cluster.order <- h$label[h$order]
    } else {
      cluster.order <- unique(node.labels)
    }
    for (network in cluster.order) {
      net.genes <- names(node.labels[node.labels %in% network])
      # Order genes within each network by their connectivity to all other genes
      order <- c(order, names(sort(colSums(abs(edge.matrix[,net.genes])), 
                                   decreasing=TRUE)))
    }
  }
  
  # Handle the unassigned nodes (if present)
  `%NIN%` <- function(a, b) a[!(`%in%`(a, b))]
  unassigned <- rownames(edge.matrix) %NIN% order
  order <- c(order, names(sort(colSums(abs(edge.matrix[,unassigned])))))
  
  return(order)
}

#' Combine edge weight matrices for a network derived from two different
#' data sources.
#'
#' Primarily for use with \code{\link{PlotNetworkHeatmap}}.
#'
#' @param edge.matrix1
#' @param edge.matrix2
#' @return
#'   A edge.matrix matrix where the upper triangle corresponds to 
#'   \code{edge.matrix1}, and the lower triangle corresponds to i
#'   \code{edge.matrix2}.
#'
ComparativeEdgeMatrix <- function(edge.matrix1, edge.matrix2) {
  # Check the two edge.matrix matrices are the same size, and contains
  # the same nodes (and that each matrix is well formed)
  stopifnot(dim(edge.matrix1) == dim(edge.matrix2))
  stopifnot(
    length(
      intersect(
        unlist(dimnames(edge.matrix1)), 
        unlist(dimnames(edge.matrix2))
      )   
    ) == nrow(edge.matrix1)
  )
  
  order <- colnames(edge.matrix1)
  combined <- edge.matrix1
  combined[lower.tri(combined)] <- edge.matrix2[order, order][lower.tri(combined)]
  return(combined)
}

#' @title Plot network edge weights as a heatmap.
#'
#' Plots the edge weights between nodes as a heatmap, optionally drawing boxes 
#' around pre-determined clusters. The \code{xlab} and \code{ylab} arguments are
#' useful when comparing the edge weights.
#'
#' @param edge.matrix a square matrix containing the edge weights between each
#'                    pair of nodes.
#' @param network.labels an optional named vector assigning each node to a 
#'                       network/cluster.
#' @param xlab x-axis heatmap label.
#' @param ylab y-axis heatmap label.
#' @param legend.main title to give to the legend denoting heatmap colors.
#' @param cluster.legend.main title to give to the legend denoting cluster 
#'        assignment.
#' @param heatmap.bins Number of bins to partition edge weights into for 
#'        plotting.
#' @param edge.weight.range The range of values edge weights can take. Defaults
#'        to \code{c(-1, 1)}.
#' @param heatmap.gradient A vector of colors to use for (or interpolate over)
#'        the heatmap.bins.
#' @param cluster.cols an optional vector of colors for denoting each cluster 
#'        on a legend.
#' @note 
#'  It expects nodes to be pre-ordered within the edge weight matrix:
#'  I.e. clusters should be in continuous blocks. Otherwise, garbage-in, 
#'  garbage-out.
#'
#' @export
PlotNetworkHeatmap <- function(edge.matrix, network.labels=NULL, 
                               xlab="", ylab="", legend.main="", 
                               cluster.legend.main="", heatmap.bins=11,
                               edge.weight.range=c(-1,1),
                               heatmap.gradient, cluster.cols) {
  if(missing(heatmap.gradient)) {
    # "RdYlBu" RColorBrewer palette, with the middle value replaced with white:
    # this gives a nicer contrast than he RdBu palette.
    heatmap.gradient <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
                          "#FFFFFF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", 
                          "#A50026")
  }
  
  # Find where one network cluster changes to the next:
  if (!is.null(network.labels)) {
    breaks <- rle(network.labels[colnames(edge.matrix)])
    break.points <- sapply(1:length(breaks$length), function(i) {
      sum(breaks$length[1:i])
    })  
    break.points <- c(0, break.points)
    if (missing(cluster.cols)) {
      cluster.cols <- 1:length(unique(network.labels))
    }
  }
  
  # Render edge.matrix heatmap
  if(!is.null(network.labels)) {
    layout(matrix(1:4, ncol=2, byrow=TRUE), widths=c(0.8, 0.2), 
           heights=c(0.8, 0.2))
  } else {
    layout(matrix(1:2, ncol=2, byrow=TRUE), widths=c(0.8, 0.2))
  }
  
  par(mar=c(4.1,6.1,4.1,4.1))
  image(x      = 0:ncol(edge.matrix),
        y      = 0:ncol(edge.matrix),
        z      = edge.matrix, 
        col    = colorRampPalette(heatmap.gradient)(heatmap.bins), 
        asp    = 1,  
        breaks = seq(edge.weight.range[1], edge.weight.range[2], 
                     length.out=heatmap.bins+1),
        axes   = FALSE,
        xlab   = "", 
        ylab   = "") 
  abline(0, 1, col="black")
  mtext(xlab, side=1, line=1)
  mtext(ylab, side=2, line=1)
  
  # Draw boxes around network clusters
  if(!is.null(network.labels)) {
    lines(rect(0, 0, ncol(edge.matrix), ncol(edge.matrix), lwd=1))  # heatmap border
    for (i in 2:length(break.points)) {
      lines(rect(xleft   = break.points[i-1],
                 ybottom = break.points[i-1],
                 xright  = break.points[i],
                 ytop    = break.points[i],
                 lwd     = 1 
      ))  
    }
  }
  
  # Plot edge.matrix Legend
  par(mar=c(12.1, 1.1, 8.1, 6.6))
  gradient.bar(nBins        = 11, 
               direction    = "y", 
               lines        = "black", 
               col.gradient = rev(heatmap.gradient), 
               bin.lab      = format(seq(edge.weight.range[2], 
                                         edge.weight.range[1], 
                                         length=heatmap.bins+1), digits=3)
  )   
  mtext(legend.main, side=3, line=1)
  
  # Plot indication of network clusters underneath heatmap
  if (!is.null(network.labels)) {
    par(mar=c(7.1, 4.6, 2.1, 2.6))
    gradient.bar(range        = c(0, ncol(edge.matrix)), 
                 break.points = break.points,
                 col          = cluster.cols,
                 bin.lab      = breaks$values,
                 main         = cluster.legend.main,
                 direction    = "x",
                 lines        = NA)   
  }
}
