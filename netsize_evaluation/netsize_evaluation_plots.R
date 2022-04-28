###########################################################
# This script makes plots that summarize netsize evaluation
#
# Usage:
#     Rscript netsize_evaluation_plots.R inputDir plotDir
#
#     inputDir: the folder that stores pathway enrichment results
#     plotDir: the folder to save plots
#
###########################################################

# use the pacman to install and load required packages
pacman::p_load("reshape2", "ggplot2", "plyr", "readr", "tools")

########### load in command arguments

inputDir <- commandArgs(trailingOnly = TRUE)[1]
plotDir <- commandArgs(trailingOnly = TRUE)[2]
dir.create(plotDir)
pathway_file <- commandArgs(trailingOnly = TRUE)[3]
pathway_file <- file_path_sans_ext(basename(pathway_file))

########### load in constants

# significance cutoff
sig_cutoff <- 0.05

########### read in pathway enrichment results

allPathfiles <- list.files(inputDir, pattern = "*_allPathway.txt")
allPathfiles <- allPathfiles[grep(pathway_file, allPathfiles)]

# all_path_lowq stores the lowest q value of each pathway's enrichment in all
# models
all_path_lowq <- c()
# sig_path stores only statistically significantly enriched pathways in all
# models
sig_path <- c()

for (allPathfile in allPathfiles) {
  this_all_path <- read_delim(file.path(inputDir, allPathfile),
                              delim = "\t", col_names = TRUE)
  # only consider the first 100 models if there are more than 100 models
  this_all_path <- this_all_path[this_all_path$seed <= 100, ]

  # get the network size from file name
  netsize <- gsub("netsize", "", unlist(strsplit(allPathfile, "_"))[1])

  # get the smallest q value per model per pathway
  this_all_path_lowq <- ddply(this_all_path, .(seed, pathway),
                              summarize, lowq = min(qvalue))
  this_all_path_lowq <- data.frame(netsize = rep(netsize,
                                                 nrow(this_all_path_lowq)),
                                   this_all_path_lowq)
  all_path_lowq <- rbind(all_path_lowq, this_all_path_lowq)

  # get only significantly enriched pathways
  this_sig_path <- this_all_path[this_all_path$qvalue <= sig_cutoff, ]
  this_sig_path <- data.frame(netsize = rep(netsize, nrow(this_sig_path)),
                              this_sig_path)
  sig_path <- rbind(sig_path, this_sig_path)
}

# reorder factor levels
all_path_lowq$netsize <- factor(all_path_lowq$netsize,
                                levels = as.character(sort(as.numeric(levels(
                                  all_path_lowq$netsize)))))
# convert to negative log10 scale
all_path_lowq$lowq <- -log10(all_path_lowq$lowq)

########### make summary plots

# plot the lowest q value for each model per pathway
pdf(file.path(plotDir, paste0(pathway_file, "_", "mostsig_node_per_pathway.pdf")),
    height = 200, width = 15)
ggplot(all_path_lowq, aes(x = factor(pathway), y = lowq, colour = netsize)) +
  geom_boxplot() + coord_flip() +
  geom_hline(aes(yintercept = -log10(sig_cutoff))) +
  labs(y = "-log10(qvalue)", x = "Pathway") + theme_bw()
dev.off()

# get number of nodes significantly associated with each pathway per model size
# per seed
netsize_seed_pathway <- melt(table(sig_path$netsize, sig_path$seed,
                                   sig_path$pathway))
colnames(netsize_seed_pathway) <- c("netsize", "seed", "pathway", "nodecount")

# plot the number of associated nodes per pathway
pdf(file.path(plotDir, paste0(pathway_file, "_",
                              "Number of associated nodes per pathway.pdf")),
    height = 200, width = 15)
ggplot(netsize_seed_pathway, aes(x = factor(pathway), y = nodecount,
                                 colour = factor(netsize))) + geom_boxplot() +
  coord_flip() + labs(x = "Pathway", y = "Number of nodes")
dev.off()

# filter out pathways not significantly associated with any nodes in a model
netsize_seed_pathway_nonzero <- netsize_seed_pathway[which(
  netsize_seed_pathway$nodecount != 0), ]

# plot how many times a pathway is covered (significantly assoicated)
pdf(file.path(plotDir, paste0(pathway_file, "_",
                              "Times of coverage per pathway.pdf")),
    height = 100, width = 15)
ggplot(netsize_seed_pathway_nonzero, aes(x = factor(pathway),
                                         fill = factor(netsize))) +
  geom_bar(position = "dodge") + coord_flip()
dev.off()

# get how many pathways a model covers (is significantly associated with)
netsize_seed_coverage <- melt(table(netsize_seed_pathway_nonzero$netsize,
                                    netsize_seed_pathway_nonzero$seed))
colnames(netsize_seed_coverage) <- c("netsize", "seed", "coverage")

# plot coverage of pathways at each model size for each random seed
pdf(file.path(plotDir, paste0(pathway_file, "_",
                              "Netsize_PathwayCoverage.pdf")),
    height = 5, width = 8)
ggplot(netsize_seed_coverage, aes(x = netsize_seed_coverage$netsize,
                                  y = netsize_seed_coverage$coverage,
                                  group = round_any(netsize_seed_coverage$netsize,
                                                    50, floor))) +
  geom_boxplot() + labs(x = "Model size", y = "Pathway coverage") +
  stat_summary(fun.y = median, geom = "line", aes(group = 1), colour = "red") +
  theme_bw()
dev.off()

# get significantly enriched pathway count per model size per node
netsize_seed_node <- melt(table(sig_path$netsize, sig_path$seed,
                                sig_path$node))
colnames(netsize_seed_node) <- c("netsize", "seed", "node", "pathwaycount")
# filter out nodes that are associated with no pathways
netsize_seed_node_nonzero <- netsize_seed_node[which(
  netsize_seed_node$pathwaycount != 0), ]
# get the pathway per node on average for each model
netsize_seed_PathPerNode <- ddply(netsize_seed_node_nonzero, .(netsize, seed),
                                  summarize, path.per.node = sum(pathwaycount)/
                                    length(pathwaycount),
                                  effective.node = length(pathwaycount))

# plot the number of pathways a node is significant to (exclude zero)
pdf(file.path(plotDir, paste0(pathway_file, "_", "Netsize_PathPerNode.pdf")),
    height = 5, width = 8)
ggplot(netsize_seed_PathPerNode, aes(x = netsize, y = path.per.node,
                                     group = round_any(netsize, 50, floor))) +
  geom_boxplot() + labs(x = "Model size",
                        y = paste0("Number of pathways a node significantly\n",
                                  "associated with (exclude zero)")) +
  theme_bw() + stat_summary(fun.y = median, geom = "line", aes(group = 1),
                            colour = "red")
dev.off()

# plot number of nodes associated with known pathways
pdf(file.path(plotDir, paste0(pathway_file, "_",
                              "Netsize_ExplainedNodeCount.pdf")),
    height = 5, width = 8)
ggplot(netsize_seed_PathPerNode, aes(x = netsize, y = effective.node,
                                     group = round_any(netsize, 50, floor))) +
  geom_boxplot() + labs(x = "Netsize", y = "Number of explained nodes") +
  theme_bw() + stat_summary(fun.y = median, geom = "line",
                            aes(group = 1), colour = "red")
dev.off()

# get the ratio of effective nodes in each model
netsize_seed_PathPerNode$ratio <- netsize_seed_PathPerNode$effective.node /
  netsize_seed_PathPerNode$netsize

# plot the percentage of effective nodes (nodes associated with known pathways)
pdf(file.path(plotDir, paste0(pathway_file, "_",
                              "Netsize_ExplainedNodeRatio.pdf")),
    height = 5, width = 8)
ggplot(netsize_seed_PathPerNode, aes(x = netsize, y = ratio,
                                     group = round_any(netsize, 50, floor))) +
  geom_boxplot() + labs(x = "Model size", y = "Ratio of explained nodes") +
  theme_bw() + stat_summary(fun.y = median, geom = "line", aes(group = 1),
                            colour = "red")
dev.off()
