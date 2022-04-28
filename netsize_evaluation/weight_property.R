# ################################################################################
# This script provides several helper functions to calculate weight properties.
#
# Usage:
#        It should be imported by other scripts.
################################################################################

pacman::p_load("tools")

get_HWG_index <- function(weight_matrix, node, HW_cutoff){
  # This function returns the index to HW genes on both positive and negative
  # side of a node.
  #
  # Inputs:
  # weight_matrix: a data.matrix that stores the weight matrix of an ADAGE model
  # node: int, index to a node
  # HW_cutoff: double, number of standard deviations from mean to be counted as
  #            hight-weight

  pos_cutoff <- mean(weight_matrix[, node]) +
    HW_cutoff * sd(weight_matrix[, node])
  neg_cutoff <- mean(weight_matrix[, node]) -
    HW_cutoff * sd(weight_matrix[, node])
  HWG_index <- seq(1, nrow(weight_matrix))
  this_node_HWG <- c(HWG_index[weight_matrix[, node] >= pos_cutoff],
                     HWG_index[weight_matrix[, node] <= neg_cutoff])
  return(this_node_HWG)
}


get.weight.properties <- function(netfiles, type, geneN, sampleN = NA,
                                  modelsize = NA, plot.cor = FALSE,
                                  plot.path = NA){
  # This function calculates several weight properties for multiple ADAGE models,
  # including the distribution of weight vector correlation, the genome coverage
  # of HW genes, the range of weight vectors, and number of HW genes.
  #
  # Inputs:
  # netfiles: a vector that stores file paths to multiple ADAGE models
  # type: a character among "netsize", "ensemble", "individual", it determines
  #       the type of meta data to store with the weight property results
  # geneN: int, number of genes in the compendium
  # sampleN: character, number of samples in the compendium
  # modelsize: character, number of nodes in these ADAGE models
  # plot.cor: boolean, whether to plot the distribution of weight correlations
  # plot.path: file path to save the weight correlation distribution plot
  #
  # Return:
  # A list that stores results of weight range, HW gene count, and HW gene
  # coverage

  # make sure type is one of "netsize", "ensemble", and "individual"
  if (!type %in% c("netsize", "ensemble", "individual")) {
    stop("type must be netsize, ensemble, or individual")
  }

  range_table <- list()
  HWcount_table <- list()
  HWG_coverage <- list()
  list_index1 <- 1
  list_index2 <- 1

  if (plot.cor) {
    pdf(plot.path, height = 3 * length(netfiles), width = 5)
    par(mfrow = c(length(netfiles), 1))
  }

  for (netfile in netfiles) {

    seed <- 1
    if (type == "netsize") {
      # get the parameter info from file name
      parameter_info <- unlist(strsplit(basename(netfile), "_"))
      seed <- parameter_info[8]
      if (grepl("subsize", netfile)) {
        samplesize <- parameter_info[12]
      } else {
        samplesize <- sampleN
      }
      metadata <- c(modelsize, samplesize, seed)
      metaname <- c("netsize", "samplesize", "seed")
    } else if (type == "ensemble") {
      metadata <- "ensemble"
      metaname <- "method"
    } else if (type == "individual") {
      metadata <- "individual"
      metaname <- "method"
    }

    # we have 1000 models for model size 300, to save time and keep the model
    # number same with other model size, we will only use the first 100 models
    if (as.numeric(seed) <= 100) {
      print(paste("...processing model", netfile))
      flush.console()

      # read in the model's weight matrix
      weight <- read_delim(netfile, delim = '\t', col_names = F, n_max = geneN,
                           skip = 2)
      weight_matrix <- data.matrix(weight)
      # number of nodes in the model
      nodeN <- ncol(weight_matrix)

      all_node_HWG_2std <- c()
      all_node_HWG_3std <- c()
      all_node_HWG_2.5std <- c()

      # loop through all nodes
      for (node in 1:nodeN) {

        # get weight range of the node
        this_node_range <- max(weight_matrix[, node]) -
          min(weight_matrix[, node])
        range_table[[list_index1]] <- c(metadata, node, this_node_range)

        # get node's high-weight genes using 2std as cutoff
        this_node_HWG <- get_HWG_index(weight_matrix, node, 2)
        all_node_HWG_2std <- c(all_node_HWG_2std, this_node_HWG)

        # get node's high-weight genes using 3std as cutoff
        this_node_HWG <- get_HWG_index(weight_matrix, node, 3)
        all_node_HWG_3std <- c(all_node_HWG_3std, this_node_HWG)

        # get node's high-weight genes using 2.5std as cutoff
        this_node_HWG <- get_HWG_index(weight_matrix, node, 2.5)
        all_node_HWG_2.5std <- c(all_node_HWG_2.5std, this_node_HWG)

        # node's HW gene count is based on using 2.5 std as cutoff
        this_node_HWcount <- length(this_node_HWG)
        HWcount_table[[list_index1]] <- c(metadata, node, this_node_HWcount)

        list_index1 <- list_index1 + 1
      }

      # the genome coverage of HW genes from all nodes
      all_node_HWG_2std <- unique(all_node_HWG_2std)
      all_node_HWG_3std <- unique(all_node_HWG_3std)
      all_node_HWG_2.5std <- unique(all_node_HWG_2.5std)
      HWG_2std_coverage <- length(unique(all_node_HWG_2std))/geneN
      HWG_3std_coverage <- length(unique(all_node_HWG_3std))/geneN
      HWG_2.5std_coverage <- length(unique(all_node_HWG_2.5std))/geneN
      HWG_coverage[[list_index2]] <- c(metadata, HWG_2std_coverage,
                                       HWG_3std_coverage,
                                       HWG_2.5std_coverage)
      list_index2  <- list_index2 + 1
    }

    if (plot.cor) {
      # the distribution of the correlation between any two genes' weight vectors
      weight_cor = cor(t(weight_matrix))
      weight_cor_vec = as.vector(weight_cor)
      plot(density(weight_cor_vec), main = basename(netfile), cex.main = 0.8)
    }
  }
  range_table <- do.call(rbind, range_table)
  colnames(range_table) <- c(metaname, "node", "range")
  HWcount_table <- do.call(rbind, HWcount_table)
  colnames(HWcount_table) <- c(metaname, "node", "HWcount")
  HWG_coverage <- do.call(rbind, HWG_coverage)
  colnames(HWG_coverage) <- c(metaname, "2std", "3std", "2.5std")

  if (plot.cor) {
    dev.off()
  }

  return(list(range = range_table, HWcount = HWcount_table,
              coverage = HWG_coverage))

}