# Executed by ensemble_construction.sh. This shell script goes through the steps
# necessary to (1) for each model, identify pathways significant after
# crosstalk removal, and (2) count the number of unique pathways covered
cores_N=2
KEGG_terms="../netsize_evaluation/pseudomonas_KEGG_terms.txt"
sigpathways_ensemble_dir="./models_sigpathways_no_crosstalk/ensemble/"
sigpathways_individ_dir="./models_sigpathways_no_crosstalk/individual/"

mkdir -p $sigpathways_ensemble_dir
mkdir -p $sigpathways_individ_dir

python ./pathway_coverage_no_crosstalk/pathway_coverage_no_crosstalk.py \
./ensemble_models/ $sigpathways_ensemble_dir \
$KEGG_terms $data_compendium --replace="network_ADAGE" --cores=$cores_N

python ./pathway_coverage_no_crosstalk/pathway_coverage_no_crosstalk.py \
../netsize_evaluation/models/300/ $sigpathways_individ_dir \
$KEGG_terms $data_compendium --replace="network_SdA" --cores=$cores_N

python ./pathway_coverage_no_crosstalk/pathway_counts_no_crosstalk.py \
$sigpathways_ensemble_dir $sigpathways_individ_dir
