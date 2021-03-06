% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concoct.R
\name{concoct}
\alias{concoct}
\title{concoct: DRUID wrapper}
\usage{
concoct(dge_matrix, tfidf_matrix, tfidf_crossproduct, num_random,
  druid_direction, fold_thr, pvalue_thr, entrez)
}
\arguments{
\item{dge_matrix}{This is a 2 column matrix for gene expression changes, where column 1 is the gene fold change and column 2 is the corresponding p-value for the fold change. NOTE: Use log2 of the fold changes as output, for example, from `limma` or `DESeq2`.}

\item{tfidf_matrix}{tf-idf matrix drug-gene matrix. Column names are Entrez IDs. Computed with \code{\link{ctfidf}}}

\item{tfidf_crossproduct}{Cross-product of the tf-idf matrix. Computed with \code{\link{crossprod_matrix}}.}

\item{num_random}{Number of random sets to be generated to calculate significance of enrichment.  Defaults to 1,000.}

\item{druid_direction}{Desired effect for DRUID to run on: "pos" mimics query phenotype, "neg" reverts query phenotype. Defaults to "neg".}

\item{fold_thr}{Threshold for the fold change to be considered. Defaults to 0 (i.e., log2(1), where fold change is not used as filter)}

\item{pvalue_thr}{Threshold for the p-value of the fold change to be considered. Defaults to 0.05.}

\item{entrez}{EntrezIDs for genes in differentially expressed set. Must be same order as the input matrix.}
}
\value{
A data frame.
}
\description{
This function is a wrapper for a start to finish to identify drugs/frug profiles that mimic or revert a given profile of interest.
For a given gene set of interest, the algorithm will compute the cosine similarity between this vector and the gene sets in a given collection of pathways (for example, KEGG pathways.)
Additionally, the algorithm will generate a set of random gene sets that will be used to determine the probability of a given enrichment being random.
The output of the wrapper is a tidy data frame with all the enrichment scores that are better than 0.
}
