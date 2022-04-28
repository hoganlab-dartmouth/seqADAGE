###########################################################
# This script does pathway enrichment analysis for all the ensemble models.
#
# Usage:
#     Run inside ensemble_construction folder
#     Rscript eADAGE_pathway_enrichment.R netfolder pathway_file data_file
#     output_all N_cores
#
#     netfolder: the directory that stores ensemble models
#     pathway_file: file path to 'pseudomonas_KEGG_terms.txt'
#     data_file: the training compendium
#     output_all: whether the full results including non-significant enrichments
#                 will be written into files
#     N_cores: number of cores to use in parallel
#     HW_cutoff: number of standard deviations from the mean to be counted as
#                high-weight
###########################################################

pacman::p_load("doParallel", "tools")
#source(file.path("..", "netsize_evaluation", "pathway_enrichment.R"))
source(file.path("..", "node_interpretation", "pathway_enrichment.R"))

####### load command arguments
commArgs <- commandArgs(trailingOnly = TRUE)
netfolder <- commArgs[1]
pathway_file <- commArgs[2]
data_file <- commArgs[3]
output_all <- as.logical(commArgs[4])
N_cores <- as.numeric(commArgs[5])
HW_cutoff <- as.numeric(commArgs[6])

# register the number of cores to use in parallel
registerDoParallel(cores = N_cores)

########### load constant

# FDR significance cutoff
sig_cutoff <- 0.05

######## pathway enrichment analysis for 4 types of models

pathway_file_base <- file_path_sans_ext(basename(pathway_file))

# for ensemble models built from the same 100 individual models, corADAGE method
#print(paste("doing pathway enrichment analysis for corADAGE models built from",
#            "the same 100 individual ADAGE models..."))
#netfiles <- list.files(netfolder,
#  pattern = glob2rx("*_1_100_*_ClusterByweight_avgweight_network_ADAGE.txt"))
#netfiles <- file.path(netfolder, netfiles)
#outfile1 <- paste0(pathway_file_base, "_", "same100_corADAGE_sigPathway.txt")
#outfile2 <- paste0(pathway_file_base, "_", "same100_corADAGE_allPathway.txt")
#multi.pathway.analysis(netfiles, data_file, pathway_file, outfile1, outfile2,
#                       output_all, "ensemble", HW_cutoff, sig_cutoff)

# for ensemble models built from the same 100 individual models, eADAGE method
print(paste("doing pathway enrichment analysis for eADAGE models built from",
            "the same 100 individual ADAGE models..."))
netfiles <- list.files(netfolder,
  pattern = glob2rx("*_660_735_*_ClusterByweighted_avgweight_network_ADAGE.txt"))
netfiles <- file.path(netfolder, netfiles)
outfile1 <- paste0(pathway_file_base, "_", "same100_eADAGE_sigPathway.txt")
outfile2 <- paste0(pathway_file_base, "_", "same100_eADAGE_allPathway.txt")
multi.pathway.analysis(netfiles, data_file, pathway_file, outfile1, outfile2,
                       output_all, "ensemble", HW_cutoff, sig_cutoff)

# for ensemble models built from different 100 individual models, corADAGE method
print(paste("doing pathway enrichment analysis for corADAGE models built from",
            "different 100 individual ADAGE models..."))
netfiles <- list.files(netfolder,
  pattern = "*_seed=1_ClusterByweight_avgweight_network_ADAGE.txt")
netfiles <- file.path(netfolder, netfiles)
outfile1 <- paste0(pathway_file_base, "_", "diff100_corADAGE_sigPathway.txt")
outfile2 <- paste0(pathway_file_base, "_", "diff100_corADAGE_allPathway.txt")
multi.pathway.analysis(netfiles, data_file, pathway_file, outfile1, outfile2,
                       output_all, "ensemble", HW_cutoff, sig_cutoff)

# for ensemble models built from different 100 individual models, eADAGE method
print(paste("doing pathway enrichment analysis for eADAGE models built from",
            "different 100 individual ADAGE models..."))
netfiles <- list.files(netfolder,
  pattern = "*_seed=1_ClusterByweighted_avgweight_network_ADAGE.txt")
netfiles <- file.path(netfolder, netfiles)
outfile1 <- paste0(pathway_file_base, "_", "diff100_eADAGE_sigPathway.txt")
outfile2 <- paste0(pathway_file_base, "_", "diff100_eADAGE_allPathway.txt")
multi.pathway.analysis(netfiles, data_file, pathway_file, outfile1, outfile2,
                       output_all, "ensemble", HW_cutoff, sig_cutoff)
