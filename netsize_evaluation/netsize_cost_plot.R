################################################################
# This script makes a plot that summarizes the training cost at different
# model size.
#
# Usage:
#     Rscript netsize_cost_plots.R cost_file plot_folder
#
#     cost_file:  the file that stores the reconstruction cost at each model
#                 size
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
pdf(file.path(plot_folder, "netsize_train_cost.pdf"), height = 3, width = 5)
ggplot(data = cost, aes(x = NetworkStructure, y = TrainCost,
                        group = round_any(cost$NetworkStructure, 50, floor))) +
  geom_boxplot() + labs(x = "Network size", y = "Reconstruction Error") +
  stat_summary(fun.y = median, geom = "line", aes(group = 1), colour = "red") +
  theme_bw()
dev.off()
