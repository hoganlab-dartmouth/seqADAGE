###############################################################################
# This script compares the significance of each pathway association between
# ensemble ADAGE models and ADAGE models
#
# Usage:
#     Rscript analyze_ensemble_significance.R
#
###############################################################################
pacman::p_load("reshape2", "ggplot2", "plyr", "readr", "tools")

########### load in command arguments
pathway_file <- commandArgs(trailingOnly = TRUE)[1]
pathway_file <- file_path_sans_ext(basename(pathway_file))

######### define constants

ADAGEPathfolder <- file.path("..", "netsize_evaluation", "models")
eADAGEPathfile <- paste0(pathway_file, "_",
                         "diff100_eADAGE_allPathway.txt")
corADAGEPathfile <- paste0(pathway_file, "_",
                           "diff100_corADAGE_allPathway.txt")
# summary_file stores the output of this script, it will be overwritten
# in each run.
summary_file <- paste0(pathway_file, "_",
                       "analyze_ensemble_significance_output.txt")
sig_cutoff <- 0.05

####### load the pathway enrichment results of all individual ADAGE models

print("reading in pathway enrichment results of ADAGE models...")
allPathfiles <- list.files(ADAGEPathfolder, pattern = "*_allPathway.txt")
allPathfiles <- allPathfiles[grep(pathway_file, allPathfiles)]
all_path_lowq <- c()
for (allPathfile in allPathfiles) {
  this_all_path <- read_delim(file.path(ADAGEPathfolder, allPathfile),
                              delim = "\t", col_names = T)
  # get the network size from file name
  netsize <- gsub("netsize", "", unlist(strsplit(allPathfile, "_"))[1])
  this_all_path_lowq <- ddply(this_all_path, .(seed, pathway), summarize,
                              lowq = min(qvalue))
  this_all_path_lowq <- data.frame(netsize = rep(netsize,
                                                 nrow(this_all_path_lowq)),
                                   this_all_path_lowq)
  all_path_lowq <- rbind(all_path_lowq, this_all_path_lowq)
}

####### load the pathway enrichment results of corADAGE models

print("reading in pathway enrichment results of corADAGE models...")
netsize <- "corADAGE"
corADAGE_all_path <- read_delim(corADAGEPathfile, delim = "\t", col_names = T)
corADAGE_all_path_lowq <- ddply(corADAGE_all_path, .(model, pathway), summarize,
                                lowq = min(qvalue))
corADAGE_all_path_lowq <- data.frame(netsize = rep(netsize,
                                                   nrow(corADAGE_all_path_lowq)),
                                     corADAGE_all_path_lowq)
corADAGE_all_path_lowq <- setNames(corADAGE_all_path_lowq, names(all_path_lowq))

####### load the pathway enrichment results of eADAGE models

print("reading in pathway enrichment results of eADAGE models...")
netsize <- "eADAGE"
eADAGE_all_path <- read_delim(eADAGEPathfile, delim = "\t", col_names = T)
eADAGE_all_path_lowq <- ddply(eADAGE_all_path, .(model, pathway), summarize,
                              lowq = min(qvalue))
eADAGE_all_path_lowq <- data.frame(netsize = rep(netsize,
                                                 nrow(eADAGE_all_path_lowq)),
                                   eADAGE_all_path_lowq)
eADAGE_all_path_lowq <- setNames(eADAGE_all_path_lowq, names(all_path_lowq))

####### plot the pathway significance across models

print("making plots ...")

# plot the lowest q value for each model per pathway
all_path_lowq <- rbind(all_path_lowq, eADAGE_all_path_lowq, corADAGE_all_path_lowq)
all_path_lowq$netsize <- factor(all_path_lowq$netsize,
                                levels = c("10", "50", "100", "200", "300","500",
                                           "750", "1000", "corADAGE", "eADAGE"))
all_path_lowq$lowq <- -log10(all_path_lowq$lowq)

pdf("netsize_corADAGE_eADAGE_mostsig_node_per_pathway.pdf", height = 200,
    width = 15)
ggplot(all_path_lowq, aes(x = factor(pathway), y = lowq, colour = netsize)) +
  geom_boxplot() + coord_flip() + geom_hline(aes(yintercept = -log10(sig_cutoff))) +
  labs(x = "", y = "-log10(qvalue)") + theme_bw()
dev.off()

# plot the lowest q value of selected pathways for each model
selected_path <- c(paste("KEGG-Pathway-pae00350: Tyrosine metabolism -",
                         " Pseudomonas aeruginosa PAO1", sep = ""),
                   paste("KEGG-Pathway-pae02040: Flagellar assembly -",
                         " Pseudomonas aeruginosa PAO1", sep = ""),
                   paste("KEGG-Pathway-pae00190: Oxidative phosphorylation -",
                         " Pseudomonas aeruginosa PAO1", sep = ""))
selected_path_lowq <- all_path_lowq[all_path_lowq$pathway %in% selected_path, ]

pdf("netsize_corADAGE_eADAGE_mostsig_node_selected_pathway.pdf", height = 10,
    width = 8)
ggplot(selected_path_lowq, aes(x = factor(pathway), y = lowq, colour = netsize)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = -log10(sig_cutoff)), linetype = "dotted") +
  labs(x = "", y = "-log10(qvalue)") + theme_bw() +
  theme(axis.text.x = element_text(angle = 50, size = 5, vjust = 0.5))
dev.off()

####### compare enrichment significance between three types of models

# get the median of pathway significance
all_path_lowq_median <- ddply(all_path_lowq, .(netsize, pathway), summarize,
                              medianq = median(lowq))
net300_path_lowq_median <- all_path_lowq_median[all_path_lowq_median$netsize ==
                                                  "300", ]
eADAGE_path_lowq_median <- all_path_lowq_median[all_path_lowq_median$netsize ==
                                                  "eADAGE", ]
corADAGE_path_lowq_median <- all_path_lowq_median[all_path_lowq_median$netsize ==
                                                    "corADAGE", ]

# compare eADAGE with 300-node ADAGE
eADAGE_ADAGE_com <- round(sum(eADAGE_path_lowq_median$medianq >
                                net300_path_lowq_median$medianq) /
                            nrow(eADAGE_path_lowq_median) * 100, 2)
print(paste(eADAGE_ADAGE_com, "percentage of pathways are better captured in",
            "eADAGE than 300-node ADAGE."))
write("comparing eADAGE with 300-node ADAGE:", summary_file)
write(paste(eADAGE_ADAGE_com, "percentage of pathways are better captured in",
            "eADAGE than 300-node ADAGE."), summary_file, append = T)

# compare corADAGE with 300-node ADAGE
corADAGE_ADAGE_com <- round(sum(corADAGE_path_lowq_median$medianq >
                                  net300_path_lowq_median$medianq) /
                              nrow(corADAGE_path_lowq_median) * 100, 2)
print(paste(corADAGE_ADAGE_com, "percentage of pathways are better captured in",
            "corADAGE than 300-node ADAGE."))
write("comparing corADAGE with 300-node ADAGE:", summary_file, append = T)
write(paste(corADAGE_ADAGE_com, "percentage of pathways are better captured in",
            "corADAGE than 300-node ADAGE."), summary_file, append = T)

# compare eADAGE with corADAGE
eADAGE_corADAGE_com <- round(sum(eADAGE_path_lowq_median$medianq >
                                   corADAGE_path_lowq_median$medianq) /
                               nrow(eADAGE_path_lowq_median) * 100, 2)
print(paste(eADAGE_corADAGE_com, "percentage of pathways are better captured in",
            "eADAGE than corADAGE."))
write("comparing eADAGE with corADAGE:", summary_file, append = T)
write(paste(eADAGE_corADAGE_com, "percentage of pathways are better captured in",
            "eADAGE than corADAGE."), summary_file, append = T)

# compare eADAGE with ADAGE at every size
all_path_lowq_median_max <- ddply(all_path_lowq_median[all_path_lowq_median$netsize !=
                                                         "corADAGE", ],
                                  .(pathway), summarize, maxq = max(medianq))
eADAGE_everyADAGE_com <- round(sum(eADAGE_path_lowq_median$medianq ==
                                     all_path_lowq_median_max$maxq) /
                                 nrow(all_path_lowq_median_max) * 100, 2)
print(paste(eADAGE_everyADAGE_com, "percentage of pathways are better captured",
            "in eADAGE than ADAGE at sizes: 10, 50, 100, 200, 300, 500, 750, 1000."))
write("comparing eADAGE with ADAGE at every size:", summary_file, append = T)
write(paste(eADAGE_everyADAGE_com, "percentage of pathways are better captured",
            "in eADAGE than ADAGE at sizes: 10, 50, 100, 200, 300, 500, 750, 1000."),
      summary_file, append = T)