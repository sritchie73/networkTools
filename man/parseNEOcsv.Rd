% Generated by roxygen2 (4.0.0): do not edit by hand
\name{parseNEOcsv}
\alias{parseNEOcsv}
\title{Converts the output of NEO into something sane.}
\usage{
parseNEOcsv(neofilename, outfile, datC)
}
\arguments{
\item{neofilename}{path to the neo csv file.}

\item{outfile}{file prefix to output the .edge.table.tsv and .snp.table.tsv
files}

\item{datC}{data matrix given to NEO as its \code{datC} argument}
}
\description{
\code{NEO} outputs a massive CSV file that is not easy to threshold nor
loadable into cytoscape. This function will parse this CSV file and write
out an edge table to file, filled with relevant columns from the CSV, along
with additional consistency checks. Some columns are themselves parse and
split to be useful, namely the \code{edge} column (into Source and Target
nodes), and the \code{Final.SNPs.LEO.NB.OCA} column into a vector of SNPs
for the Target and Source nodes respectively.
}

