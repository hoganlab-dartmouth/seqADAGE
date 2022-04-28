################################################################
# This script makes plots that summarize the training and testing cost in
# sample size evaluation. In the evaluation, samples were randomly divided
# into a training set with sizes of 100, 200, 500, and 800 and a testing set
# of 200 unused samples or whatever samples left that were not used in the
# training.
#
# Usage:
#     Rscript subsample_cost_plot.R cost_file test_size plot_folder
#
#     cost_file:  the file that stores the reconstruction cost at each sample
#                 size and each network size
#     test_size:  number of samples in the testing set, 200 in this case
#     plot_folder: the folder to save plots
################################################################

# use the pacman to install and load required packages
pacman::p_load("reshape2", "ggplot2", "plyr")

######### load in command arguments

cost_file <- commandArgs(trailingOnly = TRUE)[1]
test_size <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
plot_folder <- commandArgs(trailingOnly = TRUE)[3]
dir.create(plot_folder)

######### define constants

summary_file <- "./subsample_cost_plot_output.txt"

######### read in data

cost <- read.table(cost_file, header = TRUE, sep = "\t")
cost$diff <- cost$TestCost - cost$TrainCost
cost_fixed <- cost[cost$TestSize == test_size, ]
cost_left <- cost[cost$TestSize != test_size, ]
cost_melted <- melt(cost, id = c("NetworkStructure", "TrainSize", "TestSize",
                                 "Seed1", "Seed2", "diff"))
cost_melted_fixed <- cost_melted[cost_melted$TestSize == test_size, ]

######### make plots

pdf(file.path(plot_folder, "samplesize_train_test_cost_fixed.pdf"),
    height = 5, width = 8)
ggplot(data = cost_melted_fixed, aes(x = factor(TrainSize), y = value,
                                     fill = variable)) + geom_boxplot() +
  facet_grid(. ~ NetworkStructure) + guides(fill = guide_legend(title = "")) +
  labs(x = "Training size", y = "Cost") + theme_bw()
dev.off()

pdf(file.path(plot_folder,
              "samplesize_train_test_cost_groupby_initialization_fixed.pdf"),
    height = 5, width = 8)
ggplot(data = cost_melted_fixed, aes(x = factor(TrainSize), y = value,
                                     colour = factor(Seed1), fill = variable)) +
  geom_boxplot() + facet_grid(. ~ NetworkStructure) +
  guides(fill = guide_legend(title = "")) +
  labs(x = "Training size", y = "Cost") + theme_bw()
dev.off()

pdf(file.path(plot_folder,
              "samplesize_train_test_cost_groupby_subsample_fixed.pdf"),
    height = 5, width = 8)
ggplot(data = cost_melted_fixed, aes(x = factor(TrainSize), y = value,
                                     colour = factor(Seed2),
                                     fill = variable)) + geom_boxplot() +
  facet_grid(. ~ NetworkStructure) + guides(fill = guide_legend(title = "")) +
  labs(x = "Training size", y = "Cost") + theme_bw()
dev.off()

pdf(file.path(plot_folder, "samplesize_train_test_cost_diff_fixed.pdf"),
    height = 5, width = 8)
ggplot(data = cost_fixed, aes(x = TrainSize, y = diff)) +
  facet_grid(. ~ NetworkStructure) +
  stat_smooth(method = "lm", formula = y ~ log(x)) + geom_point()
dev.off()

######### fit an exponential model for differences in train and test costs

cost_fixed$NetworkStructure <- factor(cost_fixed$NetworkStructure)
for (netsize in levels(cost_fixed$NetworkStructure)) {
  cost_fixed_netsize <- cost_fixed[cost_fixed$NetworkStructure == netsize,
                               c("TrainSize", "diff")]
  cost_fixed_exp <- lm(log(cost_fixed_netsize$TrainSize) ~
                         cost_fixed_netsize$diff)
  intersect <- exp(cost_fixed_exp$coefficients[1])

  write(paste("fitting an exponential model for model size", netsize),
        summary_file, append = TRUE)
  write(capture.output(summary(cost_fixed_exp)), summary_file, append = TRUE)
  write(paste("Based on the exponential model, it would require approximately",
              round(intersect, 0), "samples to support a", netsize,
              "node model."), summary_file, append = TRUE)
}
