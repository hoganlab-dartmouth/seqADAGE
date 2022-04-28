"""
This script iterates through a directory of ADAGE model weight matrices and
determines the pathways significant in each model after crosstalk is removed.

Output:
    .txt tab-delimited files specifying the significant pathways
    found in the model.
    Columns:
      model size, node #, side ('pos' or 'neg'), pathway, p-value, q-value

Usage:
    pipeline.py <models-dir> <output-dir>
        <pathway-definitions> <gene-compendium>
        [--cores=<n-cores>] [--replace=<replace>]
        [--std=<std-signature>] [--alpha=<alpha>]
        [-a | --all-genes]
    pipeline.py -h | --help

Options:
    -h --help                   Show this screen.

    <models-dir>                Path to the directory containing ADAGE models
    <output-dir>                Path to the directory that will store the
                                output files. Will be created if no such path
                                exists.
    <pathway-definitions>       Path to the pathway definitions file.
                                Formatted as tab-delimited columns:
                                pathway, N (num. genes), geneA;geneB;...geneN
    <gene-compendium>           Path to the gene compendium file.

    --replace=<replace>         The resulting file keeps the naming convention
                                of the input file, with the exception of this
                                substring, which will be replaced by
                                'SigPathways'
                                [default: network]
    --cores=<n-cores>           Number of cores used to run crosstalk removal
                                on models in parallel
                                [default: num. available cores - 1]
    --std=<std-signature>       Number of standard deviations from the mean
                                gene expression value that defines the +/-
                                cutoff for an ADAGE node's high positive weight
                                and high negative weight gene signatures.
                                [default: 2.5]
    --alpha=<alpha>             Significance level for pathway enrichment.
                                [default: 0.05]

    -a --all-genes              Apply crosstalk removal to the full gene set
                                in each ADAGE node. By default, the procedure
                                is only applied to genes in the node's gene
                                signatures (genes with a high positive weight
                                or high negative weight).
                                [default: False]


"""
import multiprocessing
import os
import sys
from time import time

from docopt import docopt
from joblib import Parallel, delayed
import numpy as np
import pandas as pd

import crosstalk_removal
import crosstalk_removal_utils as utils

OUTPUT_FILE_NAMING = "SigPathway"


def process(gene_ids, pathway_definitions_map, union_pathway_genes,
            alpha, std, models_directory, current_model_filename):
    """
    Apply the crosstalk removal procedure and determine the pathway coverage
    of an ADAGE weight matrix.

    Parameters
    -----------
    gene_ids : list(str)
        The list of gene identifiers in the compendium
    pathway_definitions_map : dict(str -> list(str))
        Pathway definitions, pre-crosstalk-removal. A pathway (key) is defined
        by a list of genes (value).
    union_pathway_genes : set(str)
        The union of all genes in the list of pathway definitions
    alpha : float
        Significance level for pathway enrichment.
    std : float
        Only consider "high weight" genes: genes with weights +/- `std`
        standard deviations from the mean.
    models_directory : str
        The path to the directory containing the models
    current_model_filename : str
        The name of the model that is currently being "processed"
        (e.g. applying crosstalk removal and determining pathway enrichment)

    Returns
    -----------
    tuple(str, pandas.DataFrame), the model filename and the significant
                                  pathways dataframe
    """
    print("Processing model {0}".format(current_model_filename))
    full_filepath = os.path.join(models_directory, current_model_filename)
    weight_matrix = utils.load_weight_matrix(full_filepath, gene_ids)
    significant_pathways_df = pd.DataFrame(
        columns=["feature", "pathway", "side", "p-value", "padjust"])
    n_genes = len(gene_ids)
    for feature in weight_matrix:
        feature_df = crosstalk_removal.pathway_enrichment_without_crosstalk(
            weight_matrix[feature], alpha, std, n_genes, union_pathway_genes,
            pathway_definitions_map)
        if feature_df is not None:
            feature_df.loc[:, "feature"] = pd.Series(
                [feature] * len(feature_df.index), index=feature_df.index)
            significant_pathways_df = pd.concat(
                [significant_pathways_df, feature_df], axis=0)
    significant_pathways_df.reset_index(drop=True, inplace=True)
    return (current_model_filename, significant_pathways_df)


if __name__ == "__main__":
    arguments = docopt(
        __doc__, version="eADAGE pathway coverage without crosstalk, 0.0")
    models_directory = arguments["<models-dir>"]
    output_directory = arguments["<output-dir>"]
    pathway_definitions_file = arguments["<pathway-definitions>"]
    gene_compendium = arguments["<gene-compendium>"]
    std = float(arguments["--std"])
    alpha = float(arguments["--alpha"])
    substring_to_replace = arguments["--replace"]

    # create the output directory if it does not exist
    try:
        os.makedirs(arguments["<output-dir>"])
    except OSError:
        if not os.path.isdir(arguments["<output-dir>"]):
            raise

    if arguments["--cores"].isdigit():
        n_cores = int(arguments["--cores"])
    else:
        n_cores = multiprocessing.cpu_count() - 1

    gene_ids = utils.load_gene_identifiers(gene_compendium)
    all_defined_genes, pathway_definitions = utils.load_pathway_definitions(
        pathway_definitions_file)

    t_o = time()
    np.seterr(all="raise")
    with Parallel(n_jobs=n_cores) as parallel:
        results = parallel(
            delayed(process)(
                gene_ids, pathway_definitions, all_defined_genes,
                alpha, std, models_directory, model)
            for model in os.listdir(models_directory))
    t_f = time() - t_o
    print("{0} models took {1} seconds to run on {2} cores.".format(
        len(results), t_f, n_cores))

    for model, significant_pathways_df in results:
        output_filename = model.replace(
            substring_to_replace, OUTPUT_FILE_NAMING)
        output_path = os.path.join(output_directory, output_filename)
        significant_pathways_df.to_csv(
            path_or_buf=output_path, sep="\t", index=False)
