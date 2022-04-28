###############################################################################
# This script compares pathway coverage between ensemble models and individual
# models
#
# Usage:
#     Rscript analyze_ensemble_coverage.R
#
###############################################################################

pacman::p_load("reshape2", "ggplot2", "scales", "tools")

########### load in command arguments
pathway_file <- commandArgs(trailingOnly = TRUE)[1]
pathway_file <- file_path_sans_ext(basename(pathway_file))

######### load in constant
# summary_file stores the output of this script, it will be
# overwritten in each run.
summary_file <- paste0(pathway_file, "_",
                       "analyze_ensemble_coverage_output.txt")

######### load in function
path_coverage <- function(indiADAGE_sig_file, num_indi_models, corADAGE_sig_file,
                          eADAGE_sig_file, prefix, summary_file) {

  # This function compares the pathway coverage between ADAGE and ensemble ADAGEs
  # (corADAGE and eADAGE).

  # Arguments:
  # indiADAGE_sig_file: file that stores pathways significantly enriched
  #                     in individual ADAGE models
  # num_indi_models: number of individual models to include in the comparison
  # corADAGE_sig_file: file that stores pathways significantly enriched in
  #                    corADAGE models
  # eADAGE_sig_file: file that stores pathways significantly enriched in eADAGE
  #                  models
  # prefix: prefix to file name of output plots, must contain either 'same' or
  #         'diff'
  # summary_file: output file that saves the t.test results

  # Outputs:
  # five plots in pdf format
  # 1. plot the number of pathways covered by each model type
  # 2. plot the coverage rate per pathway between ADAGE models and corADAGE
  #    models
  # 3. plot the coverage rate per pathway between ADAGE models and eADAGE models
  # 4. plot the distribution of pathway coverage rate between ADAGE models
  #    and corADAGE models
  # 5. plot the distribution of pathway coverage rate between ADAGE models
  #    and eADAGE models
  # text output will be written into summary_file

  # load pathways covered by each individual model
  indi_sig_path <- read.table(indiADAGE_sig_file, header = T, sep = "\t",
                              stringsAsFactors = F)
  # filter the indi_sig_path to only contain the needed individual models
  indi_sig_path <- subset(indi_sig_path,
                          subset = indi_sig_path$seed <= num_indi_models)
  indi_pathway <- melt(table(indi_sig_path$seed, indi_sig_path$pathway))
  colnames(indi_pathway) <- c("seed", "pathway", "count")
  indi_pathway_nonzero <- indi_pathway[which(indi_pathway$count != 0), ]
  indi_coverage <- melt(table(indi_pathway_nonzero$seed))
  colnames(indi_coverage) <- c("seed", "coverage")

  # load pathways covered by each corADAGE model
  corADAGE_sig_path <- read.table(corADAGE_sig_file, header = T, sep = "\t",
                                  stringsAsFactors = F)
  if (grepl("same", prefix)) {
    corADAGE_pathway <- melt(table(corADAGE_sig_path$seed,
                                   corADAGE_sig_path$pathway))
  } else if (grepl("diff", prefix)) {
    corADAGE_pathway <- melt(table(corADAGE_sig_path$model,
                                   corADAGE_sig_path$pathway))
  }
  colnames(corADAGE_pathway) <- c("seed", "pathway", "count")
  corADAGE_pathway_nonzero <- corADAGE_pathway[which(corADAGE_pathway$count !=
                                                       0), ]
  corADAGE_coverage <- melt(table(corADAGE_pathway_nonzero$seed))
  colnames(corADAGE_coverage) <- c("seed", "coverage")

  # load pathways covered by each eADAGE model
  eADAGE_sig_path <- read.table(eADAGE_sig_file, header = T, sep = "\t",
                                stringsAsFactors = F)
  if (grepl("same", prefix)) {
    eADAGE_pathway <- melt(table(eADAGE_sig_path$seed, eADAGE_sig_path$pathway))
  } else if (grepl("diff", prefix)) {
    eADAGE_pathway <- melt(table(eADAGE_sig_path$model, eADAGE_sig_path$pathway))
  }
  colnames(eADAGE_pathway) <- c("seed", "pathway", "count")
  eADAGE_pathway_nonzero <- eADAGE_pathway[which(eADAGE_pathway$count != 0), ]
  eADAGE_coverage <- melt(table(eADAGE_pathway_nonzero$seed))
  colnames(eADAGE_coverage) <- c("seed", "coverage")

  # compare the number of pathways covered by each model type
  merge_coverage <- rbind(cbind(indi_coverage,
                                method = rep("individual", nrow(indi_coverage))),
                          cbind(corADAGE_coverage,
                                method = rep("corADAGE", nrow(corADAGE_coverage))),
                          cbind(eADAGE_coverage,
                                method = rep("eADAGE", nrow(eADAGE_coverage))))
  ggplot(merge_coverage, aes(x = method, y = coverage)) +
    geom_boxplot() + labs(y = "Pathway Coverage") + theme_bw()
  ggsave(paste(prefix, "pathway coverage comparison.pdf", sep = "_"), height = 5,
         width = 5)

  # t test on comparing the number of pathways covered by ADAGE models with the
  # ensemble models
  write("comparing corADAGE with ADAGE:", summary_file, append = TRUE)
  t.result <- t.test(corADAGE_coverage$coverage, indi_coverage$coverage,
                     alternative = "greater")
  write(capture.output(t.result), summary_file, append = TRUE)
  write("comparing eADAGE with ADAGE:", summary_file, append = TRUE)
  t.result <- t.test(eADAGE_coverage$coverage, indi_coverage$coverage,
                     alternative = "greater")
  write(capture.output(t.result), summary_file, append = TRUE)
  write("comparing eADAGE and corADAGE:", summary_file, append = TRUE)
  t.result <- t.test(eADAGE_coverage$coverage, corADAGE_coverage$coverage,
                     alternative = "greater")
  write(capture.output(t.result), summary_file, append = TRUE)

  # compare the coverage rate per pathway between ADAGE models and corADAGE models
  indi_pathway_count <- melt(table(indi_pathway_nonzero$pathway))
  corADAGE_pathways_count <- melt(table(corADAGE_pathway_nonzero$pathway))
  corADAGE_merged_pathways <- merge(indi_pathway_count, corADAGE_pathways_count,
                                    by = 1, all = T)
  colnames(corADAGE_merged_pathways) <- c("Pathway", "Individual", "corADAGE")
  corADAGE_merged_pathways$Individual <- corADAGE_merged_pathways$Individual /
    nrow(indi_coverage)
  corADAGE_merged_pathways$corADAGE <- corADAGE_merged_pathways$corADAGE /
    nrow(corADAGE_coverage)
  corADAGE_merged_pathways_molten <- melt(corADAGE_merged_pathways)

  ggplot(corADAGE_merged_pathways_molten,
         aes(x = Pathway, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(y = "coverage rate") + coord_flip() + theme_bw()
  ggsave(paste(prefix, "compare pathway coverage between ADAGE and corADAGE.pdf",
               sep = "_"), height = 100, width = 15, limitsize = FALSE)

  # compare the coverage rate per pathway between ADAGE models and eADAGE models
  indi_pathway_count <- melt(table(indi_pathway_nonzero$pathway))
  eADAGE_pathways_count <- melt(table(eADAGE_pathway_nonzero$pathway))
  eADAGE_merged_pathways <- merge(indi_pathway_count, eADAGE_pathways_count,
                                  by = 1, all = T)
  colnames(eADAGE_merged_pathways) <- c("Pathway", "Individual", "eADAGE")
  eADAGE_merged_pathways$Individual <- eADAGE_merged_pathways$Individual /
    nrow(indi_coverage)
  eADAGE_merged_pathways$eADAGE <- eADAGE_merged_pathways$eADAGE /
    nrow(eADAGE_coverage)
  eADAGE_merged_pathways_molten <- melt(eADAGE_merged_pathways)

  ggplot(eADAGE_merged_pathways_molten,
         aes(x = Pathway, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(y = "coverage rate") + coord_flip() + theme_bw()
  ggsave(paste(prefix, "compare pathway coverage between ADAGE and eADAGE.pdf",
               sep = "_"), height = 100, width = 15, limitsize = FALSE)

  # compare the distribution of pathway coverage rates between ADAGE and corADAGE
  corADAGE_merged_pathways[is.na(corADAGE_merged_pathways)] <- 0
  corADAGE_d <- density(corADAGE_merged_pathways$corADAGE)
  Individual_d <- density(corADAGE_merged_pathways$Individual)
  pdf(paste(prefix, "corADAGE pathway coverage rate.pdf", sep = "_"), height = 5,
      width = 8)
  plot(corADAGE_d, xlab = "coverage rate", main = "", xlim = c(0, 1), xaxs = "i")
  polygon(corADAGE_d, col = alpha("#00BFC4", 0.8), border = "#00BFC4")
  polygon(Individual_d, col = alpha("#F8766D", 0.8), border = "#F8766D")
  legend("topleft", c("corADAGE", "Individual"), lty = c(1, 1), lwd = c(2.5, 2.5),
         col = c("#00BFC4", "#F8766D"), cex = 0.75)
  dev.off()

  # compare the distribution of pathway coverage rates between ADAGE and eADAGE
  eADAGE_merged_pathways[is.na(eADAGE_merged_pathways)] <- 0
  eADAGE_d <- density(eADAGE_merged_pathways$eADAGE)
  Individual_d <- density(eADAGE_merged_pathways$Individual)
  pdf(paste(prefix, "eADAGE pathway coverage rate.pdf", sep = "_"), height = 5,
      width = 8)
  plot(eADAGE_d, xlab = "coverage rate", main = "", xlim = c(0, 1), xaxs = "i")
  polygon(eADAGE_d, col = alpha("#00BFC4", 0.8), border = "#00BFC4")
  polygon(Individual_d, col = alpha("#F8766D", 0.8), border = "#F8766D")
  legend("topleft", c("eADAGE", "Individual"), lty = c(1, 1), lwd = c(2.5, 2.5),
         col = c("#00BFC4", "#F8766D"), cex = 0.75)
  dev.off()
}

####################################################################
# compare ensemble models built from same 100 individual models with 100
# individual models

print(paste("comparing ensemble models built from the same 100 individual",
            "models with 100 individual ADAGE models..."))
write(paste("comparing ensemble models built from the same 100 individual",
            "models with 100 individual ADAGE models:"),
      summary_file)

indiADAGE_sig_file <- file.path("..", "netsize_evaluation", "models",
                                paste0("netsize300_", pathway_file,
                                       "_sigPathway.txt"))

num_indi_models <- 100
corADAGE_sig_file <- paste0(pathway_file, "_",
                            "same100_corADAGE_sigPathway.txt")
eADAGE_sig_file <- paste0(pathway_file, "_",
                          "same100_eADAGE_sigPathway.txt")
prefix <- paste0(pathway_file, "_", "same")
path_coverage(indiADAGE_sig_file, num_indi_models, corADAGE_sig_file,
              eADAGE_sig_file, prefix, summary_file)

####################################################################
# compare ensemble models built from different 100 individual models with 1000
# individual models

print(paste("comparing ensemble models built from different 100 individual",
            "models with 1000 individual ADAGE models..."))
write(paste("comparing ensemble models built from different 100 individual",
            "models with 100 individual ADAGE models:"),
      summary_file, append = TRUE)

num_indi_models <- 1000
corADAGE_sig_file <- paste0(pathway_file, "_",
                            "diff100_corADAGE_sigPathway.txt")
eADAGE_sig_file <- paste0(pathway_file, "_",
                          "diff100_eADAGE_sigPathway.txt")
prefix <- paste0(pathway_file, "_", "diff")
path_coverage(indiADAGE_sig_file, num_indi_models, corADAGE_sig_file,
              eADAGE_sig_file, prefix, summary_file)