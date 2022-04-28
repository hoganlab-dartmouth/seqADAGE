################################################################
# This script makes a plot that summarizes the training cost at different
# corruption level.
#
# Usage:
#     Rscript corruption_cost_plots.R cost_file plot_folder
#
#     cost_file:  the file that stores the reconstruction cost for each model
#     plot_folder: the folder to save plots
################################################################

pacman::p_load("ggplot2")

########### load command arguments

cost_file <- commandArgs(trailingOnly = TRUE)[1]
plot_folder <- commandArgs(trailingOnly = TRUE)[2]
dir.create(plot_folder)

########### read in data
cost <- read.table(cost_file, header = T, sep = "\t")

########### make plot
pdf(file.path(plot_folder, "corruption_train_cost.pdf"), height = 4, width = 15)
ggplot(data = cost, aes(x = CorruptLevel, y = TrainCost)) +
  geom_boxplot() + labs(x = "Corruption level", y = "Reconstruction Error") +
  stat_summary(fun.y = median, geom = "line", aes(group = 1), colour = "red") +
  theme_bw()
dev.off()
