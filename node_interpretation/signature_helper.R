#pacman::p_load("readr", "dplyr")
library(readr)
library(dplyr)

load_model <- function(model_file, geneID){
  # This function reads in the ADAGE model file and return the weight matrix
  # of the model
  # 
  # model_file: file path to an ADAGE network file
  # geneID: gene identifiers correspond to the model
  # return: a data.frame with the first column being geneID and the rest
  # columns being node weights.

  model <- readr::read_delim(model_file, delim = "\t", col_names = FALSE,
                             n_max = nrow(geneID), skip = 2)
  colnames(model) <- paste0("Node", seq(1, ncol(model)))
  model <- dplyr::bind_cols(geneID, model)
  return(model)
}


one_signature <- function(node_weight, geneID, side, HW_cutoff = 2.5){
  # This function extracts a single gene signature from an ADAGE model.
  # 
  # node_weight: a vector storing weight values of each gene to a node
  # geneID: gene identifiers correspond to the weight vector
  # side: character, "pos" or "neg"
  # HW_cutoff: number of standard deviations from mean in a node's weight
  # distribution to be considered as high-weight (default to 2.5).
  # return: a character vector storing genes in the signature defined by the
  # weight vector and the side.

  if (side == "pos") {
    
    pos_cutoff <- mean(node_weight) + HW_cutoff * sd(node_weight)
    print(mean(node_weight))
    # order positive HW genes from high to low weight
    HW_order <- order(node_weight, decreasing = TRUE)
    ordered_node_weight <- node_weight[HW_order]
    ordered_geneID <- geneID[HW_order]
    HWG <- ordered_geneID[ordered_node_weight >= pos_cutoff]
    
  } else if (side == "neg") {
    
    neg_cutoff <- mean(node_weight) - HW_cutoff * sd(node_weight)
    # order negative HW genes from low to high weight
    HW_order <- order(node_weight, decreasing = FALSE)
    ordered_node_weight <- node_weight[HW_order]
    ordered_geneID <- geneID[HW_order]
    HWG <- ordered_geneID[ordered_node_weight <= neg_cutoff]
    
  } else {
    stop("side can only be pos or neg.")
  }
  
  return(HWG)
}

extract_signatures <- function(model, HW_cutoff = 2.5){
  # This function extracts all signatures of an ADAGE model.
  #
  # model: the ADAGE model to be used for extracting signatures
  # HW_cutoff: number of standard deviations from mean in a node's weight
  # distribution to be considered as high-weight (default to 2.5). Only
  # high-weight genes are included in signatures.
  # return: a named list with each element being a gene signature

  weight_matrix <- as.matrix(model[, -1])
  model_size <- ncol(weight_matrix)
  geneID <- as.data.frame(model)[, 1]
  
  pos_signatures <- lapply(1:model_size, function(x)
    one_signature(weight_matrix[, x], geneID, "pos", HW_cutoff))
  neg_signatures <- lapply(1:model_size, function(x)
    one_signature(weight_matrix[, x], geneID, "neg", HW_cutoff))
  
  signature_list <- c(pos_signatures, neg_signatures)
  names(signature_list) <- c(paste0("Node", seq(1, model_size), "pos"),
                             paste0("Node", seq(1, model_size), "neg"))
  
  return(signature_list)
}


enrich_test <- function(set1, set2, length_all){
  # This function performs a one-side fisher exact test to determine whether
  # two sets of genes have significant overlap.
  # 
  # set1: character vector storing genes in the first gene set.
  # set2: character vector storing genes in the second gene set.
  # length_all: the number of all possible genes to choose from.
  # return: object returned by fisher.test() function, which is a list with
  # class "htest" containing p.value, conf.int, estimate, and so on.

  set_overlap <- length(intersect(set1, set2))
  set1_only <- length(set1) - set_overlap
  set2_only <- length(set2) - set_overlap
  others <- length_all - set_overlap - set2_only - set1_only
  contingency_table <- matrix(c(set_overlap, set1_only, set2_only, others),
                              nrow = 2)
  fisher_result <- fisher.test(contingency_table, alternative='greater')
  
  return(fisher_result)
}

test_signature_overlap <- function(signatures_genes){
  # This function tests how significant any two combinations of input signatures
  # overlap with each other in term of their gene compositions.
  #
  # model: an ADAGE model to extract signatures from
  # return: a named list storing adjusted p values in all signature overlap
  # tests.

  # get unique genes in all signatures
  all_genes <- unique(unlist(signatures_genes))

  # generate all possible combinations of both signatures' gene lists and names
  comb_sig <- combn(signatures_genes, 2)
  comb_name <- combn(names(signatures_genes), 2)
  
  # test the significance of overlap between every two signature combinations
  result_table <- mapply(enrich_test, comb_sig[1, ], comb_sig[2, ],
                         MoreArgs = list(length_all = length(all_genes)))
  
  # name each element with its signature combination
  colnames(result_table) <- sapply(1: ncol(comb_name),function(x)
    paste(comb_name[, x], collapse = "_"))
  
  # extract odds ratio from the enrichment results
  p_values <- result_table["p.value", ]
  
  # multiple hypothesis correction
  q_values <- p.adjust(p_values, method = "fdr")
  
  return(q_values)
  
}

