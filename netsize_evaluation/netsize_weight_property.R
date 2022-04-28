#############################################################################
# This script evaluates the properties of weight vectors at each model size.
#
# Usage:
#     Rscript weight_property_netsize.R modelDir model_sizes plot_folder
#     sampleN geneN subsample modelDir_subsample
#
#     modelDir: the folder that stores models
#     model_sizes: list of model sizes separated by comma, e.g. "50,300"
#     plot_folder: the folder to store output plots
#     sampleN: number of samples in the compendium
#     geneN: number of genes in the compendium
#     subsample: whether add models from subsampling analysis
#     modelDir_subsample: the folder that stores models built using all samples
#                         in the compendium
#############################################################################


# use the pacman to install and load required packages
pacman::p_load("reshape2", "ggplot2", "plyr", "readr", "gridExtra")
source("weight_property.R")

############ load in command arguments and constants

comArgs <- commandArgs(trailingOnly = TRUE)
modelDir <- comArgs[1]
model_sizes <- comArgs[2]
plot_folder <- comArgs[3]
sampleN <- comArgs[4]
geneN <- as.numeric(comArgs[5])
subsample <- as.logical(comArgs[6])

if (subsample) {
  modelDir_subsample <- comArgs[7]
  property.plot <- "samplesize_weight_property.pdf"
} else {
  property.plot <- "netsize_weight_property.pdf"
}
dir.create(plot_folder)

############ read in each model and calculate its weight properties

model_sizes <- unlist(strsplit(model_sizes,','))
netfolders <- file.path(modelDir, model_sizes, sep = '')

if (subsample) {
  # if the goal here is to evaluate sample size, then combine models built with
  # full compendium with models built in the subsampling analysis
  netfolders_subsample <- file.path(modelDir_subsample, model_sizes, sep = '')
  netfolders <- c(netfolders, netfolders_subsample)
}

range_table <- c()
HWcount_table <- c()
HWG_coverage <- c()

for (netfolder in netfolders) {
  # list all files that end with 'network_SdA.txt'
  netfiles <- list.files(netfolder, pattern = "*_network_SdA.txt")
  netfiles <- file.path(netfolder, netfiles)
  modelsize <- basename(netfolder)  # get the model size

  weight_properties <- get.weight.properties(netfiles, "netsize", geneN, sampleN,
                                             modelsize)
  range_table <- rbind(range_table, weight_properties$range)
  HWcount_table <- rbind(HWcount_table, weight_properties$HWcount)
  HWG_coverage <- rbind(HWG_coverage, weight_properties$coverage)

}

# build data.frame and make values numeric
range_frame <- data.frame(range_table, stringsAsFactors = FALSE)
range_frame$netsize <- as.numeric(range_frame$netsize)
range_frame$samplesize <- as.numeric(range_frame$samplesize)
range_frame$range <- as.numeric(range_frame$range)

HWcount_frame <- data.frame(HWcount_table, stringsAsFactors = FALSE)
HWcount_frame$netsize <- as.numeric(HWcount_frame$netsize)
HWcount_frame$samplesize <- as.numeric(HWcount_frame$samplesize)
HWcount_frame$HWcount <- as.numeric(HWcount_frame$HWcount)

HWG_coverage <- data.frame(HWG_coverage, stringsAsFactors = FALSE)
HWG_coverage_melted <- melt(HWG_coverage,
                            id = c("netsize", "samplesize", "seed"))
HWG_coverage_melted$netsize <- as.numeric(as.character(
  HWG_coverage_melted$netsize))
HWG_coverage_melted$samplesize <- as.numeric(as.character(
  HWG_coverage_melted$samplesize))
HWG_coverage_melted$value <- as.numeric(as.character(
  HWG_coverage_melted$value))

############ make plots to summarize weight properties

if (subsample) {
  # plot the relationship between weight properties and sample size
  range_p <- ggplot(range_frame, aes(x = factor(samplesize), y = range,
                                     fill = factor(netsize))) +
    geom_boxplot() + theme_bw() +
    labs(x = "Sample Size", y = "Range of Each Weight Vector") 

  HWcount_p <- ggplot(HWcount_frame, aes(x = factor(samplesize), y = HWcount,
                                         fill = factor(netsize))) +
    geom_boxplot() + theme_bw() +
    labs(x = "Sample Size", y = "Number of HW Genes in Each Node")

  coverage_p <- ggplot(HWG_coverage_melted, aes(x = factor(samplesize),
                                                y = value, colour = variable)) +
    geom_boxplot() + facet_wrap(~netsize) + theme_bw() +
    labs(x = "Sample Size", y = "HW Genes Coverage Rate", colour = "Cutoff")

} else{
  # plot the relationship between weight properties and model size
  range_p <- ggplot(range_frame, aes(x = netsize, y = range,
                                     group = round_any(netsize, 50, floor))) +
    geom_boxplot() + labs(x = "Model Size", y = "Range of Each Weight Vector")

  HWcount_p <- ggplot(HWcount_frame, aes(x = netsize, y = HWcount,
                                         group = round_any(netsize, 50, floor))) +
    geom_boxplot() + labs(x = "Model Size", y = "Number of HW Genes in Each Node")

  coverage_p <- ggplot(HWG_coverage_melted, aes(x = netsize, y = value,
                                                group = round_any(netsize, 50, floor))) +
    geom_boxplot() + labs(x = "Model Size", y = "HW Genes Coverage Rate",
                          colour = "Cutoff") + facet_grid(~variable)
}

pdf(file.path(plot_folder, property.plot), height = 10, width = 5)
do.call(grid.arrange, c(list(range_p, HWcount_p, coverage_p),
                        list(nrow = 3, ncol = 1, top = "weight properties")))
dev.off()
