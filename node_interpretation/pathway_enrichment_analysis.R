################################################################################
# This script does pathway enrichment analysis for one ADAGE/eADAGE model.
#
# Usage:
#     Rscript pathway_enrichment_analysis.R netfile pathway_file data_file
#     output_prefix output_all HW_cutoff
#
#     netfile: the network file of an ensemble model
#     pathway_file: file path to 'pseudomonas_KEGG_terms.txt' or
#                   'pseudomonas_GO_terms.txt'
#     data_file: the training compendium
#     output_prefix: prefix to output file
#     output_all: whether the full results including non-significant enrichments
#                 will be written into files
#     HW_cutoff: number of standard deviations from the mean to be counted as
#                high-weight
################################################################################

#source(file.path("..", "netsize_evaluation", "pathway_enrichment.R"))
source('node_interpretation/pathway_enrichment.R')
############# load in arguments
#print(commandArgs)
netfile <- commandArgs(trailingOnly = TRUE)[1]
print(netfile)
pathway_file <- commandArgs(trailingOnly = TRUE)[2]
data_file <- commandArgs(trailingOnly = TRUE)[3]
output_prefix <- commandArgs(trailingOnly = TRUE)[4]
output_all <- as.logical(commandArgs(trailingOnly = TRUE)[5])
HW_cutoff <- as.numeric(commandArgs(trailingOnly = TRUE)[6])

###### load in constant

sig_cutoff <- 0.05

############# read in data

# count the number of columns in the data file
col_n <- count.fields(data_file, sep = "\t")[1]
# read in the gene IDs from data file
geneID <- read.table(data_file, sep = ",", header = T,
                     colClasses = c("character", rep("NULL", col_n - 1)),
                     skip = 0)
colnames(geneID)[1]
# read in the pathway file
pathway <- read.table(pathway_file, sep = "\t", header = F, row.names = 1,
                      stringsAsFactors = F)
# create names of the output files
outfile1 <- paste0(output_prefix, "_sigPathway.txt")
outfile2 <- paste0(output_prefix, "_allPathway.txt")


############# run analysis

one.pathway.analysis(netfile, geneID, pathway, outfile1, outfile2, HW_cutoff,
                     sig_cutoff, output_all)