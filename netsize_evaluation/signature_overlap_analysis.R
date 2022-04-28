################################################################
# This script evaluates the relationship between corruption level and signature
# redundancy.
#
# Usage:
#     Rscript signature_redundancy_analysis.R compendium_file model_folder
#     model_size plot_folder
#
#     compendium_file: the training compendium, used to get gene IDs
#     model_folder: the folder that saves ADAGE models trained with different
#                   corruption levels
#     model_size: int, model size, can deal with only one model size
#     plot_folder: the folder to save plots
################################################################

pacman::p_load("ggplot2", "dplyr", "readr")
source(file.path("..", "node_interpretation", "signature_helper.R"))

# load command arguments --------------
commArgs <- commandArgs(trailingOnly = TRUE)
compendium_file <- commArgs[1]
model_folder <- commArgs[2]
model_size <- commArgs[3]
enrichment_file <- commArgs[4]
plot_folder <- commArgs[5]
dir.create(plot_folder)

# load in expression compendium --------------

compendium <- readr::read_tsv(compendium_file)
colnames(compendium)[1] <- "geneID"
geneID <- compendium[, 1]

# load in all models --------------

netfiles <- list.files(file.path(model_folder, model_size),
                       pattern = "*_network_SdA.txt")
# extract random seeds from file names
seeds <- sapply(netfiles, function(x) as.numeric(tail(unlist(strsplit(x, "_")),
                                                      3)[1]))
# there are 100 models at each corruption level. To save time, we will only
# examine the first 10 models at each corruption level
netfiles <- netfiles[seeds <= 10]
model_folder_netfiles <- file.path(model_folder, model_size, netfiles)
all_models <- lapply(model_folder_netfiles, function(x) load_model(x, geneID))
# extract corruption levels from file names
corrupt_levels <- sapply(netfiles, function(x) unlist(strsplit(x, "_"))[5])

# number of signatures each gene is in --------------

# extract signatures for each model
all_signatures <- lapply(all_models, function(x) extract_signatures(model = x))
# count the number of signatures a gene is in for each model
gene_sig_counts <- lapply(all_signatures,
                          function(x) data.frame(table(unlist(x))))
# get the distribution of signature counts
gene_sig_counts_freq <- lapply(gene_sig_counts,
                               function(x) data.frame(table(x$Freq)))
# join counting result of each model into one data.frame
gene_sig_counts_freq <- gene_sig_counts_freq %>% Reduce(function(dtf1, dtf2)
  dplyr::full_join(dtf1, dtf2, by = "Var1"), .)
# replace NA to 0 count
gene_sig_counts_freq[is.na(gene_sig_counts_freq)] <- 0
# rename the columns
colnames(gene_sig_counts_freq) <- c("signature_count", netfiles)
# convert signature count to numeric and sort by it
gene_sig_counts_freq$signature_count <- as.numeric(
  gene_sig_counts_freq$signature_count)
gene_sig_counts_freq <- gene_sig_counts_freq[order(
  gene_sig_counts_freq$signature_count), ]

# loop through each corruption level
Nsignature_freq_corrupt <- lapply(levels(factor(corrupt_levels)), function(x){
  # get the mean frequency of each signature count for a corruption level
  mean_freq <- gene_sig_counts_freq %>% dplyr::select(matches(x)) %>% rowMeans
  dplyr::data_frame(
    Nsignature = as.numeric(names(mean_freq)),
    freq = mean_freq,
    corrupt = rep(x, length(mean_freq))
  )
})
Nsignature_freq_corrupt <- dplyr::bind_rows(Nsignature_freq_corrupt)

pdf(file.path(plot_folder, "frequency_of_signatures_per_gene(complete).pdf"),
              width = 50, height = 4)
ggplot(data = Nsignature_freq_corrupt,
       aes(x = Nsignature, y = freq, fill = corrupt)) +
  geom_bar(stat = "identity", position = "dodge") + theme_bw()
dev.off()

pdf(file.path(plot_folder, "frequency_of_signatures_per_gene(partial).pdf"),
              width = 30, height = 4)
ggplot(data = Nsignature_freq_corrupt,
       aes(x = Nsignature, y = freq, fill = corrupt)) +
  geom_bar(stat = "identity", position = "dodge") + xlim(NA, 30) + theme_bw()
dev.off()

# number of signatures each pathway is associated with --------------

# read in pathway enrichment result
pathway_result <- readr::read_tsv(enrichment_file)
# only consider models with seeds <= 10
pathway_result <- pathway_result[pathway_result$seed <= 10, ]
# count the number of signatures a pathway is associated with
seed_corrupt_pathway_count <- pathway_result %>%
  dplyr::group_by(seed, corruption, pathway) %>% dplyr::count(pathway)
# calculate the mean and median counts
seed_corrupt_pathCount_summ <- seed_corrupt_pathway_count %>%
  dplyr::group_by(seed, corruption) %>% dplyr::summarise(meanN = mean(n),
                                                         medianN = median(n))

pdf(file.path(plot_folder, "pathway_#signature_mean.pdf"),
    width = 12, height = 4)
ggplot(data = seed_corrupt_pathCount_summ, aes(x = corruption, y = meanN)) +
  geom_boxplot() + labs(y = "average number of associated
signatures per pathway") + theme_bw()
dev.off()

# get the distribution of signature counts
seed_corrupt_pathway_count_freq <- seed_corrupt_pathway_count %>%
  dplyr::group_by(seed, corruption) %>% dplyr::count(n)

pdf(file.path(plot_folder, "pathway_#signature_freq.pdf"),
    width = 50, height = 4)
ggplot(data = seed_corrupt_pathway_count_freq, aes(x = factor(n), y = nn,
                                                   colour = corruption)) +
  geom_boxplot() + labs(x = "number of associated signatures",
                        y = "counts") + theme_bw()
dev.off()

# number of signature pairs that significantly overlap --------------

# count the number of signature pairs that have significant number of genes
# overlap
num_significant_overlap <- sapply(all_signatures, function(x)
  sum(test_signature_overlap(x) <= 0.05))
corrupt_sigOverlap <- data.frame(corrupt = corrupt_levels,
                                 sig_overlap_num = num_significant_overlap)

pdf(file.path(plot_folder, "significant_overlapped_pairs.pdf"),
    width = 12, height = 4)
ggplot(data = corrupt_sigOverlap, aes(x = corrupt, y = sig_overlap_num)) +
  geom_boxplot() + labs(x = "corruption",
                        y = "number of significantly overlapped
signature pairs") + theme_bw()
dev.off()
