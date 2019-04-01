\name{BANDITS-package}
\alias{BANDITS-package}
\alias{BANDITS}
\docType{package}
\title{
\packageTitle{BANDITS}
}
\description{
\packageDescription{BANDITS}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BANDITS}
\packageIndices{BANDITS}
}
\author{
\packageAuthor{BANDITS}

Maintainer: \packageMaintainer{BANDITS}
}
\keyword{
differential splicing, differential transcript usage, DTU, Bayesian hierarchical modelling
}
\examples{
## Preliminary information

# specify the directory of the internal data:
data_dir = system.file("extdata", package = "BANDITS")
data_dir

# load gene_to_transcript matching:
data("gene_tr_id", package = "BANDITS")
# gene_tr_id contains transcripts ids on the first column
# and the corresponding gene ids on the second column:
head(gene_tr_id)

# Specify the directory of the transcript level estimated counts.
sample_names = paste0("sample", seq_len(4))
quant_files = file.path(data_dir, sample_names, "quant.sf")
file.exists(quant_files)

# Load the transcript level estimated counts via tximport:
library(tximport)
txi = tximport(files = quant_files, type = "salmon", txOut = TRUE)
counts = txi$counts
head(counts)

# We define the design of the study: in our case we have 2 groups, 
# that we call "A" and "B" of 2 samples each.
samples_design = data.frame(sample_id = sample_names,
                            group = c("A", "A", "B", "B"))
samples_design

# The groups are defined in:
levels(samples_design$group)



## Optional (recommended): transcript pre-filtering

transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
                                         transcript_counts = counts,
                                         min_transcript_proportion = 0.01,
                                         min_transcript_counts = 10,
                                         min_gene_counts = 20)
head(transcripts_to_keep)



## Load the data:

# compute the Median estimated effective length for each transcript:
eff_len = eff_len_compute(x_eff_len = txi$length)

# specify the path to the equivalence classes:
equiv_classes_files = file.path(data_dir, sample_names, "aux_info", "eq_classes.txt")
file.exists(equiv_classes_files)

# Warning: the sample names in equiv_classes_files must have the same order
# as those in the design object, containted in samples_design.
equiv_classes_files
samples_design$sample_id

# create data and filter internally lowly abundant transcripts:
BANDITS_data = create_data(gene_to_transcript = gene_tr_id,
                           path_to_eq_classes = equiv_classes_files, eff_len = eff_len, 
                           n_cores = 2,
                           transcripts_to_keep = transcripts_to_keep)

# If transcripts pre-filtering is not wanted, 
# do not specify \\code{transcripts_to_keep} parameter.

# Filter lowly abundant genes:
BANDITS_data = filter_genes(BANDITS_data, min_counts_per_gene = 20)



## Optional (recommended): infer an informative prior for the precision parameter

# Use the same filtering criteria as in \\code{\\link{filter_transcripts}}; 
# if transcript pre-filtering is not performed, set \\code{min_transcript_proportion},
# \\code{min_transcript_counts} and \\code{min_gene_counts} to 0.

#set.seed(61217)
#precision = prior_precision(gene_to_transcript = gene_tr_id, transcript_counts = counts,
#                       min_transcript_proportion = 0.01, min_transcript_counts = 10,
#                       min_gene_counts = 20, n_cores = 2)

# load the pre-computed precision estimates:
data(precision, package = "BANDITS")

# Plot the histogram of the genewise log-precision estimates.
# The black solid line represents the normally distributed prior distribution 
# for the log-precision parameter.
plot_precision(precision)



## Test for DTU
#set.seed(61217)
#results = test_DTU(BANDITS_data = BANDITS_data,
#             prior_precision = precision$prior,
#             samples_design = samples_design,
#             R = 10^4, burn_in = 2*10^3, n_cores = 2,
#             gene_to_transcript = gene_tr_id)

# load the pre-computed results:
data("results", package = "BANDITS")
results

# Visualize the most significant Genes, sorted by gene level significance.
head(top_genes(results))

# Alternatively, gene-level results can also be sorted according to DTU_measure, 
# which is a measure of the strength of the change between the 
# average relative abundances of the two groups.
head(top_genes(results, sort_by = "DTU_measure"))

# Visualize the most significant transcripts, sorted by transcript level significance.
head(top_transcripts(results, sort_by = "transcript"))

# Visualize the convergence output for the most significant genes, 
# sorted by gene level significance.
head(convergence(results))

# We can further use the \\code{gene} function to gather all output for a specific gene:
# gene level, transcript level and convergence results.
top_gene = top_genes(results, n = 1)
gene(results, top_gene$Gene_id)

# Similarly we can use the \\code{transcript} function to gather all output 
# for a specific transcript.
top_transcript = top_transcripts(results, n = 1)
transcript(results, top_transcript$Transcript_id)

#Finally, we can plot the estimated average transcript relative expression 
# in the two groups for a specific gene via \\code{plot_proportions}.
library(ggplot2)
plot_proportions(results, top_gene$Gene_id)
}