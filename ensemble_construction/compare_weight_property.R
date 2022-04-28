###############################################################################
# This script compares weight properties between eADAGE models and ADAGE models
# including the distribution of weight vector correlation, the genome coverage
# of HW genes, the range of weight vectors, and number of HW genes.
#
# Usage:
#     Rscript compare_weight_property.R ensem_weight_folder indi_weight_folder
#     geneN
#
#     ensem_weight_folder: the folder that stores eADAGE models
#     indi_weight_folder: the folder that stores ADAGE models
#     geneN: number of genes in the compendium

###############################################################################
pacman::p_load("reshape2", "ggplot2", "readr", "gridExtra")
source(file.path("..", "netsize_evaluation", "weight_property.R"))

######## load command arguments

commArgs <- commandArgs(trailingOnly = TRUE)
ensem_weight_folder <- commArgs[1]
indi_weight_folder <- commArgs[2]
geneN <- as.numeric(commArgs[3])

######## load constants

property.plot <- './weight_property_comparison.pdf'
ensem.cor.plot <- "ensemble_weight_cor_dist.pdf"
indi.cor.plot <- "individual_weight_cor_dist.pdf"

######## weight properties of eADAGE models

ensem_weight_files <- list.files(
  ensem_weight_folder,
  pattern = '*_ClusterByweighted_avgweight_network_ADAGE.txt')
ensem_weight_files <- file.path(ensem_weight_folder, ensem_weight_files)
ensem_weight_properties <- get.weight.properties(ensem_weight_files,
                                                 "ensemble", geneN,
                                                 plot.cor = TRUE,
                                                 plot.path = ensem.cor.plot)
ensem_range_table <- ensem_weight_properties$range
ensem_HWcount_table <- ensem_weight_properties$HWcount
ensem_HWG_coverage <- ensem_weight_properties$coverage

######### weight properties of individual ADAGE models

indi_weight_files <- list.files(indi_weight_folder, pattern = '*_network_SdA.txt')
# randomly choose the same number of individual models as eADAGE models
indi_weight_files_chosen <- sample(indi_weight_files, length(ensem_weight_files))
indi_weight_files_chosen <- file.path(indi_weight_folder,
                                      indi_weight_files_chosen)
indi_weight_properties <- get.weight.properties(indi_weight_files_chosen,
                                                "individual", geneN,
                                                plot.cor = TRUE,
                                                plot.path = indi.cor.plot)
indi_range_table <- indi_weight_properties$range
indi_HWcount_table <- indi_weight_properties$HWcount
indi_HWG_coverage <- indi_weight_properties$coverage

######### combine results

range_table <- rbind(ensem_range_table, indi_range_table)
HWcount_table <- rbind(ensem_HWcount_table, indi_HWcount_table)
HWG_coverage <- rbind(ensem_HWG_coverage, indi_HWG_coverage)

######## make plots

range_frame <- data.frame(range_table)
range_frame$range <- as.numeric(as.character(range_frame$range))
range_p <- ggplot(range_frame, aes(x = method, y = range)) +
  geom_boxplot() + labs(x = 'Model Type', y = 'Range of Each Weight Vector')

HWcount_frame <- data.frame(HWcount_table)
HWcount_frame$HWcount <- as.numeric(as.character(HWcount_frame$HWcount))
HWcount_p <- ggplot(HWcount_frame, aes(x = method, y = HWcount)) +
  geom_boxplot() + labs(x = 'Model Type', y = 'Number of HW Genes in Each Node')

HWG_coverage_frame <- data.frame(HWG_coverage, check.names = FALSE)
HWG_coverage_melted <- melt(HWG_coverage_frame, id = "method")
HWG_coverage_melted$value <- as.numeric(as.character(HWG_coverage_melted$value))
coverage_p <- ggplot(HWG_coverage_melted, aes(x = method, y = value, colour = variable)) +
  geom_boxplot() + labs(x = 'Model Type', y = 'Genome Coverage of HW Genes',
                          colour = 'Cutoff')

pdf(property.plot, height = 10, width = 5)
do.call(grid.arrange, c(list(range_p, HWcount_p, coverage_p),
                        list(nrow = 3, ncol = 1, top = "weight properties")))
dev.off()
