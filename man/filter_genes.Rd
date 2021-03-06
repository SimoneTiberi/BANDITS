% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_Genes.R
\name{filter_genes}
\alias{filter_genes}
\title{Filter lowly abundant genes.}
\usage{
filter_genes(BANDITS_data, min_counts_per_gene = 10)
}
\arguments{
\item{BANDITS_data}{a 'BANDITS_data' object, created with the \code{\link{create_data}} function.}

\item{min_counts_per_gene}{the minimum number of counts compatible with a gene (across all samples).}
}
\value{
A \code{\linkS4class{BANDITS_data}} object.
}
\description{
\code{filter_genes} filters genes, according to the overall number of counts (across all samples) compatible with the gene.
The filtering also applies to groups of genes with reads/fragments compatible with >1 gene;
in this case, the number of counts considered is across all genes in the group.
}
\details{
The function inputs a 'BANDITS_data' object, and returns again a 'BANDITS_data' object after filtering genes and groups of genes.
}
\examples{
# load the pre-computed data:
data("input_data", package = "BANDITS")
input_data

# Filter lowly abundant genes:
input_data = filter_genes(input_data, min_counts_per_gene = 20)

}
\seealso{
\code{\link{filter_transcripts}}, \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}
}
\author{
Simone Tiberi \email{simone.tiberi@uzh.ch}
}
