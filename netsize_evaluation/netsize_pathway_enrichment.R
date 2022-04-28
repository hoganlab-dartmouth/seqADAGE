###########################################################
# This script carries out pathway enrichment analysis for each model inside
# a folder for netsize evaluation.
#
# Usage:
#     Run inside netsize_evaluation folder
#     Rscript netsize_pathway_enrichment.R modelDir model_sizes pathway_file
#     data_file output_all analysis_type N_cores HW_cutoff
#
#     modelDir: the folder that stores models
#     model_sizes: list of network sizes separated by comma, e.g. "50,300"
#     pathway_file: 'pseudomonas_KEGG_terms.txt'
#     data_file: the training compendium, used to get gene IDs
#     output_all: 'TRUE' or 'FALSE', whether the full results including
#                 non-significant enrichments will be written into files
#     analysis_type: can be "netsize", "corruption", or "subsample"
#     N_cores: number of cores to use in parallel
#     HW_cutoff: number of standard deviations from the mean to be counted as
#                high-weight
###########################################################


# use the pacman to install and load required packages
pacman::p_load("doParallel", "tools")
source("pathway_enrichment.R")

########### load command arguments
commArgs <- commandArgs(trailingOnly = TRUE)
modelDir <- commArgs[1]
model_sizes <- commArgs[2]
pathway_file <- commArgs[3]
data_file <- commArgs[4]
output_all <- as.logical(commArgs[5])
analysis_type <- commArgs[6]
N_cores <- as.numeric(commArgs[7])
HW_cutoff <- as.numeric(commArgs[8])

# register the number of cores to use in parallel
registerDoParallel(cores = N_cores)

########### load constants

# FDR significance cutoff
sig_cutoff <- 0.05

########### perform pathway enrichment analysis for each model in netfolders

model_sizes <- unlist(strsplit(model_sizes,','))
netfolders <- file.path(modelDir, model_sizes, sep = '')

for (netfolder in netfolders) {

  if (analysis_type == "subsample") {
    # list all files that end with 'network_SdA.txt' and also have "subsize"
    netfiles <- list.files(netfolder,
                           pattern = glob2rx("*_subsize_*_network_SdA.txt"))
  } else{
    # list all files that end with 'network_SdA.txt'
    netfiles <- list.files(netfolder, pattern = "*_network_SdA.txt")
  }

  netfiles <- file.path(netfolder, netfiles)
  netsize <- basename(netfolder)  #get the model size
  outfile1 <- file.path(modelDir, paste0("netsize", netsize, "_",
                        file_path_sans_ext(basename(pathway_file)), "_sigPathway.txt"))
  outfile2 <- file.path(modelDir, paste0("netsize", netsize, "_",
                        file_path_sans_ext(basename(pathway_file)), "_allPathway.txt"))
  print(paste("...processing model size ", netsize))
  flush.console()

  multi.pathway.analysis(netfiles = netfiles, data_file = data_file,
                         pathway_file = pathway_file, outfile1 = outfile1,
                         outfile2 = outfile2, output_all = output_all,
                         type = analysis_type, HW_cutoff = HW_cutoff,
                         sig_cutoff = sig_cutoff)

}
