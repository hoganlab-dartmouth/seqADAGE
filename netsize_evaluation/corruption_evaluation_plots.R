###########################################################
# This script makes plots that summarize pathway enrichment results of models
# trained with different corruption levels
#
# Usage:
#     Rscript corruption_evaluation_plots.R enrichment_file plot_folder
#
#     enrichment_file: the file that stores pathway enrichment results
#     plot_folder: the folder to save plots
#
###########################################################

pacman::p_load("ggplot2", "dplyr", "readr")

# load command arguments --------------

commArgs <- commandArgs(trailingOnly = TRUE)
enrichment_file <- commArgs[1]
plot_folder <- commArgs[2]
dir.create(plot_folder)

# define constants  --------------

sig_cutoff <- 0.05

# pathway coverage --------------

# read in pathway enrichment result
pathway_result <- readr::read_tsv(enrichment_file)
# get significant result only
pathway_result_sig <- pathway_result[pathway_result$qvalue <= sig_cutoff, ]
# count the number of signatures a pathway is associated with
seed_corrupt_pathway_count <- pathway_result_sig %>%
  dplyr::group_by(seed, corruption, pathway) %>% dplyr::count(pathway)

# calculate pathway coverage of each model
seed_corrupt_coverage <- seed_corrupt_pathway_count %>%
  dplyr::group_by(seed, corruption) %>%
  dplyr::summarise(coverage = length(unique(pathway)))

pdf(file.path(plot_folder, "corruption_pathway_coverage.pdf"),
    width = 15, height = 3)
ggplot(data = seed_corrupt_coverage, aes(x = corruption, y = coverage)) +
  geom_boxplot() + theme_bw()
dev.off()

# pathway significance --------------

# get the highest -log10(qvalue) for each pathway within a model
seed_corrupt_pathway_lowq <- pathway_result %>%
  dplyr::group_by(seed, corruption, pathway) %>%
  dplyr::summarise(lowq = max(-log10(qvalue)))

# plot the lowest q value for each model per pathway
pdf(file.path(plot_folder, "corruption_mostsig_enrichment_per_pathway.pdf"),
    height = 200, width = 15)
ggplot(data = seed_corrupt_pathway_lowq, aes(x = factor(pathway), y = lowq,
                                             colour = corruption)) +
  geom_boxplot() + coord_flip() +
  geom_hline(aes(yintercept = -log10(sig_cutoff))) +
  labs(x = "", y = "-log10(qvalue)") + theme_bw()
dev.off()

# get the median(max(-log10(q value))) within a corruption level
corrupt_pathway_median_lowq <- seed_corrupt_pathway_lowq %>%
  dplyr::group_by(corruption, pathway) %>%
  dplyr::summarise(median_lowq = median(lowq))

# get the max(median(max(-log10(q value)))) across corruption levels
pathway_best_median_lowq <- corrupt_pathway_median_lowq %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(best_median_lowq = max(median_lowq))

# calculate the percentage of pathways that were best captured under each
# corruption level
corrupt_levels <- levels(factor(corrupt_pathway_median_lowq$corruption))
best_percentage <- sapply(corrupt_levels, function(x){
  sum(pathway_best_median_lowq$best_median_lowq ==
        corrupt_pathway_median_lowq[corrupt_pathway_median_lowq$corruption == x,
                                    "median_lowq"]) /
    nrow(pathway_best_median_lowq)
})

corruption_best_percentage <- dplyr::data_frame(corruption = corrupt_levels,
                                                best_percentage = best_percentage)

pdf(file.path(plot_folder, "corruption_most_significant_percentage.pdf"),
    height = 3, width = 15)
ggplot(data = corruption_best_percentage,
       aes(x = corruption, y = best_percentage)) + geom_bar(stat = "identity") +
  labs(x = "corruption level", y = "percentage of pathways that were best
captured")
dev.off()