#' Converts the output of \code{NEO} into something sane.
#' 
#' \code{NEO} outputs a massive CSV file that is not easy to threshold nor
#' loadable into cytoscape. This function will parse this CSV file and write
#' out an edge table to file, filled with relevant columns from the CSV, along
#' with additional consistency checks. Some columns are themselves parse and 
#' split to be useful, namely the \code{edge} column (into Source and Target 
#' nodes), and the \code{Final.SNPs.LEO.NB.OCA} column into a vector of SNPs 
#' for the Target and Source nodes respectively.
#' 
#' @param neofilename path to the neo csv file.
#' @param outfile file prefix to output the .edge.table.tsv and .snp.table.tsv
#'        files
#' @param datC data matrix given to NEO as its \code{datC} argument
#' @export
parseNEOcsv <- function(neofilename, outfile, datC) {
  neo.csv <- read.csv(neofilename, stringsAsFactors=FALSE)
  snps <- strsplit(gsub('"', '', neo.csv$Final.SNPs.LEO.NB.OCA), ";")
  if (!("edge" %in% colnames(neo.csv))) {
    warning("parsing incomplete neo log file. Your run may have failed, or is",
            " still in progress.")
    colnames(neo.csv)[1] <- "edge"
  }
  dt <- data.frame(
    Source = sapply(strsplit(neo.csv$edge, " -> "), `[`, 1),
    Target = sapply(strsplit(neo.csv$edge, " -> "), `[`, 2),
    CPA.Anchors = sapply(snps, function(x) gsub("fsnp.cpa=", "", x[1])),
    OCA.Anchors = sapply(snps, function(x) gsub("fsnp.oca=", "", x[2])),
    LEO.NB.OCA = neo.csv$LEO.NB.OCA,
    LEO.NB.CPA = neo.csv$LEO.NB.CPA,
    P.Model.AB = neo.csv$Model.P.value.AtoB,
    # This is LEO.NB.OCA * P.Model.AB, and how NEO internally ranks edges.
    P.weighted.LEO.NB.OCA = neo.csv$P.weighted.LEO.NB.OCA, 
    Z.WaldPath.AB = neo.csv$ZtestPathCoefficientAB,
    # Useful for Cytoscape: A negative sign indicates an inhibitory relationship 
    Wald.Path.Sign = sign(neo.csv$ZtestPathCoefficientAB), 
    P.WaldPath.AB = 2*pnorm(abs(neo.csv$ZtestPathCoefficientAB), lower.tail=FALSE),
    BLV = neo.csv$BLV.or.BilayerZscore,
    RMSEA = neo.csv$RMSEA.2m,
    PearsonCor = neo.csv$PearsonCor
  )
  
  # ignore the SNP-trait rows, we calculate more useful information.
  edge.table <- evaluateAnchors(
                  dt[!is.na(dt$CPA.Anchors) | !is.na(dt$OCA.Anchors),], 
                  datC
                )
  write.table(edge.table, 
              file=paste(outfile, ".edge.table.tsv", sep=""), 
              quote=FALSE, sep="\t")
}

#' Filter a NEO edge table
#' 
#' Uses a combination of recommended thresholds, along with other consistency
#' checks, to confidently determine which edges are significant.
#' 
#' @details
#' {
#' NEO uses partial correlation to assess the probability of each model between
#' the Source and Target nodes. The closer the partial correlation for a model
#' is to 0, the more evidence there is for that model (1). However, the partial
#' correlation can also be close to 0 when the CPA and OCA SNPs NEO chooses are 
#' not truly causal of their respective traits. Although one would expect the
#' \code{LEO.NB.OCA} score to be close to 0, we find that many edges will pass
#' the recommended threshold even when the SNPs are not at all causal of the 
#' traits (NEO simply choses the best SNPs, but imposes no threshold on what 
#' "best" means. In other words, Garbage-in, garbage-out). Therefore, we impose
#' additional consistency checks on the relationship between the causal anchors
#' NEO choses and their trait, based on the fundamental literature on 
#' Instrumental Variable analysis utilised in the context of genetics 
#' (Mendelian Randomisation) (3):
#' 
#' If the Source is truly causal of the Target, we expect the CPA SNPs to be
#' significantly associated with the Source and Target nodes (the relationship
#' with the Target node should be entirely mediated by the Source node, and thus
#' weaker). We also expect to see the OCA SNPs to NOT be significantly 
#' associated with the Source node; the flow of information should prevent any
#' correlation (unless they are pleiotropicly associated with both nodes, in
#' which case we are less confident in the results, but you may wish to relax
#' this constraint) (3, 4). 
#' 
#' If the CPA SNPs are not significantly associated with the Source node OR the
#' OCA SNPs are not significantly associated with the Target node (or both), the
#' result should be discarded completely, and you should review your pre-NEO 
#' SNP selection and genotype loading code. In the case only one set of SNPs is
#' significant, the Next Best score will be unreliable, because NEO will be 
#' biased towards the models which use the non-significant SNPs.
#' }
#' 
#' @references
#' \enumerate{
#'  \item{
#'    Aten, J. E., Fuller, T. F., Lusis, A. J. & Horvath, S. Using genetic 
#'    markers to orient the edges in quantitative trait networks: the NEO 
#'    software. BMC Syst. Biol. 2, 34 (2008). 
#'  }
#'  \item{
#'    Aten, J. E. Causal not Confounded Gene Networks: Inferring Acyclic and 
#'    Non-acyclic Gene Bayesian Networks in mRNA Expression Studies. (2008). 
#'    at <http://labs.genetics.ucla.edu/horvath/aten/NEO/jason_e_aten_dissertation2008.pdf> 
#'  }
#'  \item{
#'    Lawlor, D. A., Harbord, R. M., Sterne, J. A. C., Timpson, N. & Davey 
#'    Smith, G. Mendelian randomization: using genes as instruments for making 
#'    causal inferences in epidemiology. Stat. Med. 27, 1133–1163 (2008). 
#'  }
#'  \item{
#'    Pearl, J. Causality: models, reasoning and inference. 
#'    (Cambridge Univ Press, 2000). 
#'  }
#' }
#' @param edge.table edge.table generated by \code{\link{parseNEOcsv}}.
#' @param min.OCA Primary recommended threshold for significant edges by NEO (1).
#' @param min.CPA The OCA score is more reliable than the CPA score, but they
#'                should be in agreement (1).
#' @param min.model.P Each model is assessed by a partial correlation: the
#'                    closer to 0, the more support for that model.
#'                    \code{P.Model.AB} is the p-value for the partial 
#'                    correlation for the causal model A --> B, as a consistency
#'                    check, this should be non-significant.
#' @param max.path.P P-value from the standardised Wald Path Coefficient, which
#'                   supports the existence of an edge (1). This should be
#'                   significant.
#' @param min.BLV The Bi-layer V score measures the presence of confounding
#'                paths. Recommended to be above 0 (2).
#' @param max.RMSEA The RMSEA is a structural edge model (sem) model-fitting 
#'                  indice (See the \code{\link{sem}} package \code{NEO} 
#'                  depends on). This should be below 0.05 (1).
#' @param max.CPA.Source.P Significance threshold for CPA SNPs on the Source 
#'                         node.
#' @param max.CPA.Target.P Significance threshold for CPA SNPs on the Target
#'                         node. Note: We expect this p-value to be much higher
#'                         than the \code{CPA.Source.P}.
#' @param max.OCA.Target.P Significance threshold for the OCA SNPs on the Target
#'                         node.
#' @param min.OCA.Source.P The relationship between the OCA SNPs and the Source
#'                         node should be non-significant if the Source is 
#'                         causal of the Target.
#' @return a subset of the edge table containing only significant edges.
#' @export
filterNEO <- function(edge.table, 
                      min.OCA = 0.3, 
                      min.CPA = 0,
                      min.model.P = 0.05,
                      max.path.P = 0.05,
                      min.BLV = 0,
                      max.RMSEA = 0.05,
                      max.CPA.Source.P = 5e-8,
                      max.CPA.Target.P = 0.05,
                      max.OCA.Target.P = 5e-8,
                      min.OCA.Source.P = 0.05
                      ) {
  edge.table[
    # Recommended Thresholds from Aten et al's paper (1)
    edge.table$LEO.NB.OCA > min.OCA & 
    edge.table$LEO.NB.CPA > min.CPA & 
    # P-value weighted by LEO.NB.OCA. Not clear to me this is a good filter?
    #edge.table$P.weighted.LEO.NB.OCA > 0.05 & 
    edge.table$P.WaldPath.AB < max.path.P &
    edge.table$RMSEA < max.RMSEA & 
    # Recommended Threshold from Aten's dissertation (2)
    edge.table$BLV > min.BLV &
    # Consistency checks with MR theory (3)
    edge.table$CPA.Source.P < max.CPA.source.P &
    edge.table$CPA.Target.P < max.CPA.target.P &
    edge.table$P.Model.AB > min.model.P & 
    # Consistency Checks with knowledge about Causality (4)
    edge.table$CPA.Target.P > edge.table$CPA.Source.P &
    edge.table$OCA.Target.P < max.OCA.target.P &
    edge.table$OCA.Source.P > min.OCA.source.P &
    # Sanity Check
    edge.table$CPA.Anchors != "c()" &
    edge.table$OCA.Anchors != "c()",
  ]
}

#' Parse an anchor column to obtain all unique SNPs
#' 
#' @param anchors the CPA.Anchors or OCA.Anchors column to parse
#' @return a list of all SNPs that appear in the column
#' @export
parseAnchors <- function(anchors) {
  unique(unlist(lapply(anchors, parseAnchor)))
}

#' Parse the entry of an anchor column
parseAnchor <- function(anchor) {
  x <- gsub("c\\(", "", anchor)
  x <- gsub("\\)", "", x)
  strsplit(x, ",")[[1]]
}

#' Evalute the Quality of the Causal and Orthoganl Anchors chosen by NEO.
#' 
#' @param dt edge table to be written out by \code{\link{parseNEOcsv}}.
#' @param datC matrix given to NEO
#' @param qtls a combined table of QTL associations generated by \code{plink} 
#'             for each trait. Table should have a column named \code{PHEN}
#' @return 
#'   \code{edge.table} with the following columns added:
#'   \itemize{
#'     \item{CPA.P}{
#'      Model p-value for a linear model of the CPA anchors on the Source node,
#'      in the subset of individuals in \code{datC}.
#'     }
#'     \item{CPA.R2}{
#'      Adjusted R-squared for a linear model of the CPA anchors on the Source
#'      node, in the subset of individuals in \code{datC}.
#'     }
#'     \item{OCA.P}{
#'      Model p-value for a linear model of the OCA anchors on the Target node,
#'      in the subset of individuals in \code{datC}.
#'     }
#'     \item{OCA.R2}{
#'      Adjusted R-squared for a linear model of the OCA anchors on the Target
#'      node, in the subset of individuals in \code{datC}.
#'     }
#'   }
#' 
evaluateAnchors <- function(edge.table, datC) {
  CPA.Source.P  <- rep(NA, nrow(edge.table))
  CPA.Source.R2 <- rep(NA, nrow(edge.table))
  CPA.Target.P  <- rep(NA, nrow(edge.table))
  CPA.Target.R2 <- rep(NA, nrow(edge.table))
  OCA.Source.P  <- rep(NA, nrow(edge.table))
  OCA.Source.R2 <- rep(NA, nrow(edge.table))
  OCA.Target.P  <- rep(NA, nrow(edge.table))
  OCA.Target.R2 <- rep(NA, nrow(edge.table))
  for (i in 1:nrow(edge.table)) {
    cpa.to.source = anchor.lm(datC, edge.table[i, "Source"], edge.table[i,"CPA.Anchors"])
    cpa.to.target = anchor.lm(datC, edge.table[i, "Target"], edge.table[i,"CPA.Anchors"])
    oca.to.source = anchor.lm(datC, edge.table[i, "Source"], edge.table[i,"OCA.Anchors"])
    oca.to.target = anchor.lm(datC, edge.table[i, "Target"], edge.table[i,"OCA.Anchors"])
    if (!is.null(cpa.to.source)) {
      f <- summary(cpa.to.source)$fstat
      CPA.Source.P[i]  = pf(f[1], f[2], f[3], lower.tail=F)
      CPA.Source.R2[i] = summary(cpa.to.source)$adj.r.squared
    }
    if (!is.null(cpa.to.target)) {
      f <- summary(cpa.to.target)$fstat
      CPA.Target.P[i]  = pf(f[1], f[2], f[3], lower.tail=F)
      CPA.Target.R2[i] = summary(cpa.to.target)$adj.r.squared
    }
    if (!is.null(oca.to.source)) {
      f <- summary(oca.to.source)$fstat
      OCA.Source.P[i]  = pf(f[1], f[2], f[3], lower.tail=F)
      OCA.Source.R2[i] = summary(oca.to.source)$adj.r.squared
    }
    if (!is.null(oca.to.target)) {
      f <- summary(oca.to.target)$fstat
      OCA.Target.P[i]  = pf(f[1], f[2], f[3], lower.tail=F)
      OCA.Target.R2[i] = summary(oca.to.target)$adj.r.squared
    }
  }
  data.frame(
    edge.table,
    CPA.Source.P = CPA.Source.P,
    CPA.Target.P = CPA.Target.P,
    OCA.Target.P = OCA.Target.P,
    OCA.Source.P = OCA.Source.P,
    CPA.Source.R2 = CPA.Source.R2,
    CPA.Target.R2 = CPA.Target.R2,
    OCA.Target.R2 = OCA.Target.R2,
    OCA.Source.R2 = OCA.Source.R2
  )
}

# Run a linear model between a trait and its causal anchors.
anchor.lm <- function(datC, trait, anchor) {
  anchors <- parseAnchor(anchor)
  if (length(anchors) > 0) {
    call <- paste("lm(", trait, " ~ ", paste(anchors, collapse=" + "), 
                   ", dat=as.data.frame(datC))", sep="")
    eval(parse(text=call))
  } else {
    NULL
  }
}

#' Get a list of relationships that NEO was unable to orient.
#' 
#' In some cases, NEO is unable to make an attempt at orienting an edge.
#' These edges will be missing completely from the edge.table, but it is useful
#' to know which edges could not be determined. 
#' 
#' @param edge.table The full edge table output by \code{\link{parseNEOcsv}}.
#' @return an edge.table of non-orientable edges.
#' @export
getNonOrientableEdges <- function(edge.table) {
  if (all(edge.table[,"LEO.NB.OCA"] > 0)) {
    stop("It appears you have provided a filtered edge.table.")
  }
  all.nodes <- unique(c(edge.table$Source, edge.table$Target))
  all.rels <- expand.grid(all.nodes, all.nodes)
  all.rels <- all.rels[all.rels[,1] != all.rels[,2],]
  all.rels <- apply(all.rels, 1, paste, collapse=" ")
  edge.rels <- apply(edge.table[,c("Source", "Target")], 1, paste, collapse=" ")
  missing <- all.rels[!(all.rels %in% edge.rels)]
  data.frame(
    Source = sapply(strsplit(missing, " "), `[`, 1),
    Target = sapply(strsplit(missing, " "), `[`, 2)
  )
}

#' Identify Nodes for whom no edges could be confidently called by NEO
#' 
#' @param edge.table The full edge.table
#' @param filtered The edge.table filtered by \code{\link{filterNEO}}.
#' @return a list of nodes with no edges
#' @export
getOrphans <- function(edge.table, filtered) {
  all.nodes <- unique(edge.table$Source)
  filtered.nodes <- unique(c(filtered$Source, filtered$Target))
  all.nodes[!(all.nodes %in% filtered.nodes)]
}

#' Compare two filtered edge tables for differences
#' 
#' Provides a list of edges that are only present in one filtered table, along
#' with a list of edges whose direction conflicts between both.
#' 
#' @param filtered1 filtered table of edges.
#' @param filtered2 filtered table of edges.
#' @return
#'   A list with the following elements:
#'   \itemize{
#'    \item{only.in.first}{A data frame of the edges only in \code{filtered1}}.
#'    \item{only.in.second}{A data frame of the edges only in \code{filtered2}}.
#'    \item{conflicting}{
#'      A data frame of edges whose direction is opposite between 
#'      \code{filtered1} and \code{filtered2}.
#'    }
#'   }
#' @export
compareFiltered <- function(filtered1, filtered2) {
  ret <- list()
  e1 <- apply(filtered1[,c("Source", "Target")], 1, paste, collapse=" ")
  e2 <- apply(filtered2[,c("Source", "Target")], 1, paste, collapse=" ")
  oe1 <- strsplit(e1[!(e1 %in% e2)], " ")
  oe2 <- strsplit(e2[!(e2 %in% e1)], " ")
  ret$only.in.first <- data.frame(Source = sapply(oe1, `[`, 1), 
                                  Target = sapply(oe1, `[`, 2))
  ret$only.in.second <- data.frame(Source = sapply(oe2, `[`, 1), 
                                   Target = sapply(oe2, `[`, 2))
  e2r <- apply(filtered2[,c("Target", "Source")], 1, paste, collapse=" ")
  conf <- strsplit(e1[e1 %in% e2], " ")
  
  ret$conflicting <- data.frame(Node1 = sapply(conf, `[`, 1),
                                Node2 = sapply(conf, `[`, 2))
  ret
}

#' Estimate the effect size of a causal edge.
#' 
#' Uses the Two-stage least squares method to estimate the causal effect size
#' of each Source of each Target using the CPA SNPs.
#' 
#' @references
#' \enumerate{
#'  \item{
#'    Angrist, J. D. & Imbens, G. W. Two-Stage Least Squares Estimation of 
#'    Average Causal Effects in Models with Variable Treatment Intensity. 
#'    J. Am. Stat. Assoc. 90, 431–442 (1995). 
#'  }
#'  \item{
#'    Lawlor, D. A., Harbord, R. M., Sterne, J. A. C., Timpson, N. & Davey 
#'    Smith, G. Mendelian randomization: using genes as instruments for making 
#'    causal inferences in epidemiology. Stat. Med. 27, 1133–1163 (2008). 
#'  }
#' }
#' 
#' @param filtered A table of significant edges to estimate the effect size for.
#' @param datC A \code{data.frame} containing the observations for each SNP and
#'             trait.
#' @return 
#' Returns the \code{data.frame} given in to the \code{filtered} parameter with
#' three additional columns: one for the estimated effect size, one with the
#' p-value, and one with the r-squared (proportion of variance explained).
#' @export
estimateCausalEffectSize <- function(filtered, datC) {
  effectSize <- rep(NA, nrow(filtered))
  effectPval <- rep(NA, nrow(filtered))
  effectR2 <- rep(NA, nrow(filtered))
  for (i in 1:nrow(filtered)) {
    stage1 <- anchor.lm(datC, filtered[i, "Source"], filtered[i,"CPA.Anchors"])
    stage2 <- lm(datC[,filtered[i, "Target"]] ~ stage1$fitted.values)
    effectSize[i] <- coef(summary(stage2))[2, "Estimate"]
    effectPval[i] <- coef(summary(stage2))[2, "Pr(>|t|)"]
    effectR2[i] <- summary(stage2)$r.squared
  }
  data.frame(filtered, 
             Causal.Effect.Size = effectSize, 
             Causal.Effect.P.value = effectPval,
             Causal.Effect.R2 = effectR2)
}
