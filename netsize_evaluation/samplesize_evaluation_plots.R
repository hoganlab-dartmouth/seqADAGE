################################################################
# This script makes plots that summarize pathway coverage of each model in
# the subsampling analysis.
#
# Usage:
#     Rscript samplesize_evaluation_plots.R dataDir dataDir_full model_sizes
#     compendium_size plot_dir
#
#     dataDir: the folder that stores pathway enrichment results of the
#              subsampling analysis
#     dataDir_full: the folder that stores pathway enrichment results with
#                   the full sample size
#     model_sizes: list of network sizes separated by comma, e.g. "50,300"
#     compendium_size: number of samples in the compendium
#     plot_dir: the folder to save plots
################################################################

# use the pacman to install and load required packages
pacman::p_load("reshape2", "ggplot2", "plyr")

########## load command arguments

dataDir <- commandArgs(trailingOnly = TRUE)[1]
dataDir_full <- commandArgs(trailingOnly = TRUE)[2]
model_sizes <- commandArgs(trailingOnly = TRUE)[3]
compendium_size <- commandArgs(trailingOnly = TRUE)[4]
plot_dir <- commandArgs(trailingOnly = TRUE)[5]

########## load constant
# summary_file stores the output of this script, it will be overwritten
# in each run.
summary_file = "./samplesize_evaluation_plots_output.txt"

########## process command arguments

model_sizes <- unlist(strsplit(model_sizes,','))
subsample_files <- file.path(dataDir,paste("netsize", model_sizes,
                                           "_sigPathway.txt", sep = ''))
compendium_files <- file.path(dataDir_full,paste("netsize", model_sizes,
                                            "_sigPathway.txt", sep = ''))

########## loop through each model size and perform analyses

write("Here are pathway coverage results for the subsampling analysis:\n",
      summary_file)

netsize_samplesize_merged <- c()
for (i in 1:length(model_sizes)) {

  # read in pathway enrichment results of models built for the subsampling
  # analysis
  sig_path <- read.table(subsample_files[i], header = T, sep = "\t")
  size_seed_pathway <- melt(table(sig_path$samplesize, sig_path$seed2,
                                  sig_path$seed1, sig_path$pathway))
  size_seed_pathway_nonzero <- size_seed_pathway[which(
    size_seed_pathway$value != 0), ]
  # get pathway coverage of each model
  size_seed_coverage <- melt(table(size_seed_pathway_nonzero$Var1,
    size_seed_pathway_nonzero$Var2, size_seed_pathway_nonzero$Var3))
  colnames(size_seed_coverage) <- c("samplesize", "seed1", "seed2", "coverage")

  # plot the pathway coverage boxplot grouped by subsampling seed
  ggplot(data = size_seed_coverage, aes(x = factor(samplesize),
                                        y = coverage,
                                        fill = factor(seed1))) +
    geom_boxplot() + theme_bw() + labs(x = "Sample Size",
                                       y = "Pathway Coverage") +
    guides(fill = guide_legend(title = "subsample\ngroup"))
  ggsave(file.path(plot_dir,
                   paste0("net", model_sizes[i],
                         "_SampleSize_PathwayCoverage_groupby_subsample.pdf")),
      height = 5, width = 8)

  # plot the pathway coverage boxplot grouped by model initialization seed
  ggplot(data = size_seed_coverage, aes(x = factor(samplesize),
                                        y = coverage,
                                        fill = factor(seed2))) +
    geom_boxplot() + theme_bw() + labs(x = "Sample Size",
                                       y = "Pathway Coverage") +
    guides(fill = guide_legend(title = "initialization\ngroup"))
  ggsave(file.path(plot_dir,
                paste0("net", model_sizes[i],
                       "_SampleSize_PathwayCoverage_groupby_initialization.pdf")),
         height = 5, width = 8)

  # read in pathway enrichment results of models built with full sample size
  sig_path_compendium <- read.table(compendium_files[i], header = T, sep = "\t")
  # we built 1000 models at model size 300, only using the first 100 models is
  # enough here
  sig_path_compendium <- subset(sig_path_compendium,
    subset = sig_path_compendium$seed <= 100)
  seed_pathway <- melt(table(sig_path_compendium$seed,
    sig_path_compendium$pathway))
  seed_pathway_nonzero <- seed_pathway[which(seed_pathway$value != 0), ]
  seed_coverage <- melt(table(seed_pathway_nonzero$Var1))
  colnames(seed_coverage) <- c("seed2", "coverage")
  seed_coverage$samplesize <- rep(compendium_size, 100)
  seed_coverage$seed1 <- rep(1, 100)

  # combine the subsampling coverage and full coverage
  coverage_merged <- rbind(size_seed_coverage, seed_coverage)

  # test the significance in the pathway coverage increase
  aov.fit <- aov(coverage ~ factor(samplesize), data = coverage_merged)
  posthoc <- TukeyHSD(x = aov.fit, conf.level = 0.95)
  write(paste("Testing the significance of pathway coverage increase",
              "at model size", model_sizes[i], ":"),
        summary_file, append = TRUE)
  write(capture.output(posthoc), summary_file, append = TRUE)

  coverage_merged$netsize <- rep(model_sizes[i], nrow(coverage_merged))
  netsize_samplesize_merged <- rbind(netsize_samplesize_merged, coverage_merged)
}

netsize_samplesize_merged$samplesize <- as.numeric(
  netsize_samplesize_merged$samplesize)

# plot the relationship between sample size and pathway coverage
pdf(file.path(plot_dir, "SampleSize_PathwayCoverage.pdf"), height = 5,
  width = 8)
ggplot(data = netsize_samplesize_merged, aes(x = samplesize, y = coverage,
  group = round_any(netsize_samplesize_merged$samplesize, 100, floor))) +
  facet_grid(.~netsize) + geom_boxplot() + theme_bw() +
  labs(x = "Sample Size", y = "Pathway Coverage") +
  stat_summary(fun.y = median, geom = "line", aes(group = 1), colour = "red") +
  scale_x_continuous(breaks = c(100, 200, 500, 800, 1000))
dev.off()
