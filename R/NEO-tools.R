#' Converts the output of \code{NEO} into something sane.
#' 
#' \code{NEO} outputs a massive CSV file that is not easily parseable, nor
#' loadable into cytoscape. This function will parse this CSV file and write
#' two tables to file in tab separatted format: 
#' An edge table, filled with relevant columns, notably
#' splitting out the \code{edge} column into Source and Target nodes, and 
#' splitting the \code{Final.SNPs.LEO.NB.OCA} column into a vector of SNPs for 
#' the target and source nodes respectively, and a second table with SNP to Node 
#' relationships.
#' 
#' @param neofilename path to the neo csv file.
#' @param outfile file prefix to output the .edge.table.tsv and .snp.table.tsv
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
    P.WaldPath.AB = 2*pnorm(abs(neo.csv$ZtestPathCoefficientAB), lower.tail=FALSE),
    BLV = neo.csv$BLV.or.BilayerZscore,
    PearsonCor = neo.csv$PearsonCor
  )
  
  write.table(df[!is.na(df$CPA.Anchors),], 
              file=paste0(outfile, ".edge.table.tsv"), 
              quote=FALSE, sep="\t")
  write.table(df[is.na(df$CPA.Anchors), c("Source", "Target", "PearsonCor")],
              file=paste0(outfile, ".snp.table.tsv"),
              quote=FALSE, sep="\t")
}

#' Get recommended thresholds for filtering NEO edges
#'
#' @export
getNEOthresholds <- function() {
  thresholds=list(
    LEO.NB.OCA > 0.3,
    LEO.NB.CPA > 0,
    P.Model.AB > 0.05, # Verify
    P.weighted.LEO.NB.OCA > 0.05, # Verify
    P.WaldPath.AB < 0.05,
    BLV > 0
  )
}
