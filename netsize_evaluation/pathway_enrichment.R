################################################################################
# This script provides several helper functions for pathway enrichment analysis
#
# Usage:
#        It should be imported by other scripts.
################################################################################

pacman::p_load("readr")
source(file.path("..", "node_interpretation", "signature_helper.R"))

enrichment.test <- function(gene.set, pathway.terms, geneN){
  # This function tests pathway enrichment for a gene set on a group of
  # pathway terms.
  #
  # Inputs:
  # gene.set: a vector contains genes in the gene set to test
  # pathway.terms: a data.frame read from the pseudomonas_KEGG_terms.txt
  # geneN: number of all possible genes

  pvalue_list <- c()
  # loop through pathway terms
  for (term in 1:nrow(pathway.terms)) {
    # get genes in one pathway term
    term_genes <- unlist(strsplit(pathway.terms[term, 2], ";"))
    test_result <- enrich_test(gene.set, term_genes, geneN)
    # get p value
    pvalue <- test_result$p.value
    pvalue_list <- c(pvalue_list, pvalue)
  }
  return(pvalue_list)
}


test.one.model <- function(weight_matrix, pathway, HW_cutoff) {
  # This function tests pathway enrichment for all nodes in an ADAGE model.
  #
  # Inputs:
  # weight_matrix: a data.matrix that stores weights of an ADAGE model
  # pathway: a data.frame read from the pseudomonas_KEGG_terms.txt
  # HW_cutoff: the number of standard deviations from mean to be considered as
  #            high-weight
  # Output:
  # pvalue_table: a data.frame that stores the enrichment results between each
  #               node and each pathway

  nodeN <- ncol(weight_matrix)
  geneN <- nrow(weight_matrix)

  pvalue_table <- list()
  for (node in 1:nodeN) {

    # get high-weight genes at pos and neg sides
    pos_cutoff <- mean(weight_matrix[, node]) +
      HW_cutoff * sd(weight_matrix[, node])
    neg_cutoff <- mean(weight_matrix[, node]) -
      HW_cutoff * sd(weight_matrix[, node])
    HW_genes_pos <- rownames(weight_matrix)[weight_matrix[, node] >= pos_cutoff]
    HW_genes_neg <- rownames(weight_matrix)[weight_matrix[, node] <= neg_cutoff]

    # pathway enrichment analysis for positive high-weight genes
    pvalue_list <- enrichment.test(HW_genes_pos, pathway, geneN)
    this_table_pos <- data.frame(node = paste("Node", rep(node, nrow(pathway))),
                                 side = rep("pos", nrow(pathway)),
                                 pathway = rownames(pathway),
                                 pvalue = pvalue_list)

    # pathway enrichment analysis for negative high-weight genes
    pvalue_list <- enrichment.test(HW_genes_neg, pathway, geneN)
    this_table_neg <- data.frame(node = paste("Node", rep(node, nrow(pathway))),
                                 side = rep("neg", nrow(pathway)),
                                 pathway = rownames(pathway),
                                 pvalue = pvalue_list)

    # combine pos and neg sides
    this_table <- rbind(this_table_pos, this_table_neg)

    # FDR correction
    this_table$pvalue <- as.numeric(as.character(this_table$pvalue))
    this_table$qvalue <- p.adjust(this_table$pvalue, method = "fdr")
    pvalue_table[[node]] <- this_table
  }
  pvalue_table <- do.call(rbind, pvalue_table)

  return (pvalue_table)
}


one.pathway.analysis <- function(netfile, geneID, pathway, outfile1,
                                 outfile2, HW_cutoff, sig_cutoff,
                                 component = FALSE, output_all = FALSE){
  # This function reads in an ADAGE model from a file, performs pathway
  # enrichment analysis, and write the output into files.
  #
  # Inputs:
  # netfile: file path to the network file of an ADAGE model
  # geneID: the gene IDs that correspond to each row of an ADAGE weight matrix
  # pathway: a data.frame read from the pseudomonas_KEGG_terms.txt
  # outfile1: file path to output file that stores only significant enrichment
  #           results
  # outfile2: file path to output file that stores all tested node-pathway
  #           associations
  # HW_cutoff: float, number of standard deviations from mean to be considered
  #            as high-weight
  # sig_cutoff: float, significance cutoff of FDR q-values
  # component: boolean, whether each column in the weight matrix represent a
  #            component in PCA or ICA model
  # output_all: boolean, whether to write outfile2 that stores every enrichment
  #             test

  # load the weight matrix
  geneN = nrow(geneID)
  weight = read_delim(netfile, delim = '\t', col_names = F, n_max = geneN,
                      skip = 2)
  weight_matrix <- data.matrix(weight[1:geneN, ])
  rownames(weight_matrix) = geneID$Gene_symbol

  # test pathway association between all nodes in a model and all pathways
  pvalue_table <- test.one.model(weight_matrix, pathway, HW_cutoff)

  if (component) {
    # rename node to component for PCA and ICA models
    colnames(pvalue_table)[1] = 'component'
  }

  # get significant results and write to file
  pvalue_table_sig = pvalue_table[which(pvalue_table$qvalue <= sig_cutoff), ]
  write.table(pvalue_table_sig, outfile1, col.names = T, row.names = F,
              quote = F, sep = '\t')

  if (output_all) {
    # output all results including non significant associations
    write.table(pvalue_table, outfile2, col.names = T, row.names = F, quote = F,
                sep = '\t')
  }
}


multi.pathway.analysis <- function(netfiles, data_file, pathway_file, outfile1,
                                   outfile2, output_all, type,
                                   HW_cutoff = 2.5, sig_cutoff = 0.05) {
  # This function handles multiple ADAGE network files together. It is designed
  # to perform pathway enrichment analyses for models with different model
  # sizes (type="netsize"), or trained with different number of samples
  # (type="subsample"), or multiple ensemble models (type="ensemble"). It will
  # assemble enrichment results with meta data and write to output files.
  #
  # Inputs:
  # netfiles: file paths to network files of ADAGE models to be analyzed
  # data_file: file path to the input data compendium that are used to build
  #            those ADAGE models
  # pathway_file: file path to the pathway annotation file
  #               pseudomonas_KEGG_terms.txt
  # outfile1: file path to output file that stores only significant enrichment
  #           results
  # outfile2: file path to output file that stores all tested node-pathway
  #           associations
  # output_all: boolean, whether to write outfile2 that stores every enrichment
  #             test
  # type: can be "netsize", "corruption", "subsample", "ensemble", or
  #       "ensembleSize", it determines how the meta data are recorded
  # HW_cutoff: float, number of standard deviations from mean to be considered
  #            as high-weight
  # sig_cutoff: float, significance cutoff of FDR q-values


  # check input analysis type
  if (!type %in% c("netsize", "corruption", "subsample", "ensemble",
                   "ensembleSize")){
    stop("Analysis type unknown! Must be netsize, corruption, subsample,
         ensemble, or ensembleSize.")
  }

  # count the number of columns in the data file
  col_n <- count.fields(data_file, sep = "\t")[1]
  # read in the gene IDs from data file
  geneID <- read.table(data_file, sep = "\t", header = T,
                       colClasses = c("character", rep("NULL", col_n - 1)))
  geneN <- nrow(geneID)
  # read in pathway terms
  pathway <- read.table(pathway_file, sep = "\t", header = F, row.names = 1,
                        stringsAsFactors = F)

  # the pathway enrichment analysis for each network file runs in parallel
  all_path <- foreach(i = 1:length(netfiles), .combine = "rbind") %dopar% {

    weight = read_delim(netfiles[i], delim = '\t', col_names = F, n_max = geneN,
                        skip = 2)
    weight_matrix <- data.matrix(weight[1:geneN, ])
    rownames(weight_matrix) <- geneID$Gene_symbol

    pvalue_table <- test.one.model(weight_matrix, pathway, HW_cutoff)

    # get the parameter info from file name
    parameter_info <- unlist(strsplit(basename(netfiles[i]), "_"))
    if (type == "subsample") {
      seed1 <- parameter_info[8]
      seed2 <- parameter_info[10]
      samplesize <- parameter_info[12]
      pvalue_table <- data.frame(samplesize = rep(samplesize, nrow(pvalue_table)),
                                 seed2 = rep(seed1, nrow(pvalue_table)),
                                 seed1 = rep(seed2, nrow(pvalue_table)),
                                 pvalue_table)
    } else if (type == "netsize") {
      seed <- parameter_info[8]
      pvalue_table <- data.frame(seed = rep(seed, nrow(pvalue_table)),
                                 pvalue_table)
    } else if (type == "corruption") {
      seed <- parameter_info[8]
      corruption <- parameter_info[5]
      pvalue_table <- data.frame(seed = rep(seed, nrow(pvalue_table)),
                                 corruption = rep(corruption, nrow(pvalue_table)),
                                 pvalue_table)
    } else if (type == "ensemble") {
      model  <- parameter_info[3]
      seed <- parameter_info[6]
      pvalue_table <- data.frame(model = rep(model, nrow(pvalue_table)),
                                 seed = rep(seed, nrow(pvalue_table)),
                                 pvalue_table)
    } else if (type == "ensembleSize") {
      ensembleSize  <- parameter_info[4]
      seed <- parameter_info[6]
      pvalue_table <- data.frame(ensembleSize = rep(ensembleSize,
                                                    nrow(pvalue_table)),
                                 seed = rep(seed, nrow(pvalue_table)),
                                 pvalue_table)
    }

    pvalue_table
  }

  all_path$qvalue <- as.numeric(all_path$qvalue)
  sig_path <- all_path[which(all_path$qvalue <= sig_cutoff), ]
  # write the significant result table
  write.table(sig_path, outfile1, col.names = TRUE, row.names = FALSE,
              quote = FALSE, sep = "\t")
  if (output_all) {
    # write the complete result table
    write.table(all_path, outfile2, col.names = TRUE, row.names = FALSE,
                quote = FALSE, sep = "\t")
  }
}

