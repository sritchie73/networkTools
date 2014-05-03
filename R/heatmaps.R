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
      net.genes <- names(node.labels[which(node.labels == network)])
      # Order genes within each network by their connectivity to all other genes
      order <- c(order, names(sort(colSums(abs(edge.matrix[,net.genes])))))
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
  combined[lower.tri(combined)] <- edge.matrix2[order, order][lower.tri(edge.matrix2)]
  return(combined)
}

#' @title Plot network edge weights as a heatmap.
#'
#' Plots the edge weights between nodes as a heatmap, optionally drawing boxes 
#' around pre-determined clusters. The \code{xlab} and \code{ylab} arguments are
#' useful when comparing the edge weights 
#'
#' @param edge.matrix a square matrix containing the edge weights between each
#'                    pair of nodes.
#' @param network.labels an optional named vector assigning each node to a 
#'                       network/cluster.
#' @param xlab 
#' @param ylab
#' @note 
#'  It expects nodes to be pre-ordered within the edge weight matrix:
#'  I.e. clusters should be in continuous blocks. Otherwise, garbage-in, 
#'  garbage-out.
#'
#' @export
PlotNetworkHeatmap <- function(edge.matrix, network.labels=NULL, xlab="", ylab="") {
  # "RdYlBu" ColorBrewer palette, with the middle value replaced with white:
  # this gives a nicer contrast than he RdBu palette.
  gradient <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFFF", 
                "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  
  # Find where one network cluster changes to the next:
  if (!is.null(network.labels)) {
    breaks <- rle(network.labels[colnames(edge.matrix)])
    break.points <- sapply(1:length(breaks$length), function(i) {
      sum(breaks$length[1:i])
    })  
    break.points <- c(0, break.points)
  }
  
  # Render edge.matrix heatmap
  if(!is.null(network.labels)) {
    layout(matrix(1:4, ncol=2, byrow=TRUE), widths=c(0.8, 0.2), heights=c(0.8, 0.2))
  } else {
    layout(matrix(1:2, ncol=2, byrow=TRUE), widths=c(0.8, 0.2))
  }
  
  par(mar=c(4.1,6.1,4.1,4.1))
  image(x      = 0:ncol(edge.matrix),
        y      = 0:ncol(edge.matrix),
        z      = edge.matrix, 
        col    = gradient, 
        asp    = 1,  
        breaks = seq(-1, 1, length.out=12),
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
               col.gradient = rev(gradient), 
               bin.lab      = format(seq(1, -1, length=12), digits=3)
  )   
  mtext("Edge Weight", side=3, line=1)
  
  # Plot indication of network clusters underneath heatmap
  if (!is.null(network.labels)) {
    par(mar=c(7.1, 4.6, 2.1, 2.6))
    gradient.bar(range        = c(0, ncol(edge.matrix)), 
                 break.points = break.points,
                 col          = .ClusterColors[1:length(breaks$values)],
                 bin.lab      = breaks$values,
                 main         = "Modules",
                 direction    = "x",
                 lines        = NA)   
  }
}

# Taken from the WGCNA package
.ClusterColors <- c("turquoise", "blue", "brown", "yellow", "green", "red", "black", 
                    "pink", "magenta", "purple", "greenyellow", "tan", "salmon", 
                    "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen", 
                    "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", 
                    "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", 
                    "steelblue", "paleturquoise", "violet", "darkolivegreen", "darkmagenta", 
                    "sienna3", "yellowgreen", "skyblue3", "plum1", "orangered4", 
                    "mediumpurple3", "lightsteelblue1", "lightcyan1", "ivory", "floralwhite", 
                    "darkorange2", "brown4", "bisque4", "darkslateblue", "plum2", 
                    "thistle2", "thistle1", "salmon4", "palevioletred3", "navajowhite2", 
                    "maroon", "lightpink4", "lavenderblush3", "honeydew1", "darkseagreen4", 
                    "coral1", "antiquewhite4", "coral2", "mediumorchid", "skyblue2", 
                    "yellow4", "skyblue1", "plum", "orangered3", "mediumpurple2", 
                    "lightsteelblue", "lightcoral", "indianred4", "firebrick4", "darkolivegreen4", 
                    "brown2", "blue2", "darkviolet", "plum3", "thistle3", "thistle", 
                    "salmon2", "palevioletred2", "navajowhite1", "magenta4", "lightpink3", 
                    "lavenderblush2", "honeydew", "darkseagreen3", "coral", "antiquewhite2", 
                    "coral3", "mediumpurple4", "skyblue4", "yellow3", "sienna4", 
                    "pink4", "orangered1", "mediumpurple1", "lightslateblue", "lightblue4", 
                    "indianred3", "firebrick3", "darkolivegreen2", "blueviolet", 
                    "blue4", "deeppink", "plum4", "thistle4", "tan4", "salmon1", 
                    "palevioletred1", "navajowhite", "magenta3", "lightpink2", "lavenderblush1", 
                    "green4", "darkseagreen2", "chocolate4", "antiquewhite1", "coral4", 
                    "mistyrose", "slateblue", "yellow2", "sienna2", "pink3", "orangered", 
                    "mediumpurple", "lightskyblue4", "lightblue3", "indianred2", 
                    "firebrick2", "darkolivegreen1", "blue3", "brown1", "deeppink1", 
                    "powderblue", "tomato", "tan3", "royalblue3", "palevioletred", 
                    "moccasin", "magenta2", "lightpink1", "lavenderblush", "green3", 
                    "darkseagreen1", "chocolate3", "aliceblue", "cornflowerblue", 
                    "navajowhite3", "slateblue1", "whitesmoke", "sienna1", "pink2", 
                    "orange4", "mediumorchid4", "lightskyblue3", "lightblue2", "indianred1", 
                    "firebrick", "darkgoldenrod4", "blue1", "brown3", "deeppink2", 
                    "purple2", "tomato2", "tan2", "royalblue2", "paleturquoise4", 
                    "mistyrose4", "magenta1", "lightpink", "lavender", "green2", 
                    "darkseagreen", "chocolate2", "antiquewhite", "cornsilk", "navajowhite4", 
                    "slateblue2", "wheat3", "sienna", "pink1", "orange3", "mediumorchid3", 
                    "lightskyblue2", "lightblue1", "indianred", "dodgerblue4", "darkgoldenrod3", 
                    "blanchedalmond", "burlywood", "deepskyblue", "red1", "tomato4", 
                    "tan1", "rosybrown4", "paleturquoise3", "mistyrose3", "linen", 
                    "lightgoldenrodyellow", "khaki4", "green1", "darksalmon", "chocolate1", 
                    "antiquewhite3", "cornsilk2", "oldlace", "slateblue3", "wheat1", 
                    "seashell4", "peru", "orange2", "mediumorchid2", "lightskyblue1", 
                    "lightblue", "hotpink4", "dodgerblue3", "darkgoldenrod1", "bisque3", 
                    "burlywood1", "deepskyblue4", "red4", "turquoise2", "steelblue4", 
                    "rosybrown3", "paleturquoise1", "mistyrose2", "limegreen", "lightgoldenrod4", 
                    "khaki3", "goldenrod4", "darkorchid4", "chocolate", "aquamarine", 
                    "cyan1", "orange1", "slateblue4", "violetred4", "seashell3", 
                    "peachpuff4", "olivedrab4", "mediumorchid1", "lightskyblue", 
                    "lemonchiffon4", "hotpink3", "dodgerblue1", "darkgoldenrod", 
                    "bisque2", "burlywood2", "dodgerblue2", "rosybrown2", "turquoise4", 
                    "steelblue3", "rosybrown1", "palegreen4", "mistyrose1", "lightyellow4", 
                    "lightgoldenrod3", "khaki2", "goldenrod3", "darkorchid3", "chartreuse4", 
                    "aquamarine1", "cyan4", "orangered2", "snow", "violetred2", "seashell2", 
                    "peachpuff3", "olivedrab3", "mediumblue", "lightseagreen", "lemonchiffon3", 
                    "hotpink2", "dodgerblue", "darkblue", "bisque1", "burlywood3", 
                    "firebrick1", "royalblue1", "violetred1", "steelblue1", "rosybrown", 
                    "palegreen3", "mintcream", "lightyellow3", "lightgoldenrod2", 
                    "khaki1", "goldenrod2", "darkorchid2", "chartreuse3", "aquamarine2", 
                    "darkcyan", "orchid", "snow2", "violetred", "seashell1", "peachpuff2", 
                    "olivedrab2", "mediumaquamarine", "lightsalmon4", "lemonchiffon2", 
                    "hotpink1", "deepskyblue3", "cyan3", "bisque", "burlywood4", 
                    "forestgreen", "royalblue4", "violetred3", "springgreen3", "red3", 
                    "palegreen1", "mediumvioletred", "lightyellow2", "lightgoldenrod1", 
                    "khaki", "goldenrod1", "darkorchid1", "chartreuse2", "aquamarine3", 
                    "darkgoldenrod2", "orchid1", "snow4", "turquoise3", "seashell", 
                    "peachpuff1", "olivedrab1", "maroon4", "lightsalmon3", "lemonchiffon1", 
                    "hotpink", "deepskyblue2", "cyan2", "beige", "cadetblue", "gainsboro", 
                    "salmon3", "wheat", "springgreen2", "red2", "palegreen", "mediumturquoise", 
                    "lightyellow1", "lightgoldenrod", "ivory4", "goldenrod", "darkorchid", 
                    "chartreuse1", "aquamarine4", "darkkhaki", "orchid3", "springgreen1", 
                    "turquoise1", "seagreen4", "peachpuff", "olivedrab", "maroon3", 
                    "lightsalmon2", "lemonchiffon", "honeydew4", "deepskyblue1", 
                    "cornsilk4", "azure4", "cadetblue1", "ghostwhite", "sandybrown", 
                    "wheat2", "springgreen", "purple4", "palegoldenrod", "mediumspringgreen", 
                    "lightsteelblue4", "lightcyan4", "ivory3", "gold3", "darkorange4", 
                    "chartreuse", "azure", "darkolivegreen3", "palegreen2", "springgreen4", 
                    "tomato3", "seagreen3", "papayawhip", "navyblue", "maroon2", 
                    "lightsalmon1", "lawngreen", "honeydew3", "deeppink4", "cornsilk3", 
                    "azure3", "cadetblue2", "gold", "seagreen", "wheat4", "snow3", 
                    "purple3", "orchid4", "mediumslateblue", "lightsteelblue3", "lightcyan3", 
                    "ivory2", "gold2", "darkorange3", "cadetblue4", "azure1", "darkorange1", 
                    "paleturquoise2", "steelblue2", "tomato1", "seagreen2", "palevioletred4", 
                    "navy", "maroon1", "lightsalmon", "lavenderblush4", "honeydew2", 
                    "deeppink3", "cornsilk1", "azure2", "cadetblue3", "gold4", "seagreen1", 
                    "yellow1", "snow1", "purple1", "orchid2", "mediumseagreen", "lightsteelblue2", 
                    "lightcyan2", "ivory1", "gold1")
