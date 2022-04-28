# file path to the data compendium built in the previous step
data_compendium="../data_collection/all-pseudomonas-gene-normalized.pcl"
# number of samples in the compendium
compendiumSize=$(head -n1 $data_compendium | grep -o "\t" | wc -l)
# number of genes in the compendium
N_genes=$(expr $(wc -l <$data_compendium) - 1)

# the KEGG pathways used in the paper were downloaded on Aug. 12 2015 and is
# included in the repository named as "pseudomonas_KEGG_terms.txt".
# can also provide GO terms "manual_GO_BP_terms.txt" downloaded from
# pseudomonas.com on Feb. 27 2017
pathway_file="./pseudomonas_KEGG_terms.txt"

# number of cores to use in parallel in netsize_pathway_enrichment.R,
# please modify this to accommodate your need.
N_cores=2
# number of standard deviations from the mean to be counted as high-weight
HW_cutoff=2.5

# folder to save output plots
plot_folder="./plots"
# folder to save ADAGE models trained with different model sizes
models_folder="./models"
# folder to save ADAGE models trained with different sample sizes
subsample_models_folder="./subsample_models"
# folder to save ADAGE models trained with different corruption levels
corruption_models_folder="./corruption_models"

##########################
# download KEGG/GO pathways
##########################
# You can use the following script to download an updated version of
# KEGG terms. If you want to use the updated KEGG, make sure you also change
# the file path to pathway_file above.
python get_annotations.py ./pseudomonas_KEGG_terms_updated.txt KEGG
# The following script is used to download the GO terms from
# pseudomonas.com and process into the right format
Rscript process_GO.R ./manual_GO_BP_terms_updated.txt

###########################################
# evaluate the influence of model size
###########################################
# build 100 ADAGE models at each model size 10,50,100,200,300,500,750,1000.
# The corruption level is set to 0.1.
# This step is very time-consuming. We highly recommend configuring the
# run_ADAGE_train.py so that the construction of each individual ADAGE model is
# run as one job and distributed across a computing cluster.
python run_ADAGE_train.py $data_compendium $models_folder \
10,50,100,200,300,500,750,1000 0.1 1,101

# pathway enrichment analysis for each model at each model size
# This evaluation could run overnight.
Rscript netsize_pathway_enrichment.R $models_folder \
10,50,100,200,300,500,750,1000 $pathway_file $data_compendium TRUE netsize \
$N_cores $HW_cutoff

# make plots that summarize model size evaluation results
Rscript netsize_evaluation_plots.R $models_folder $plot_folder $pathway_file

# evaluate the properties of weight vectors at each model size
Rscript netsize_weight_property.R $models_folder 10,50,100,200,300,500,750,1000 \
$plot_folder $compendiumSize $N_genes FALSE

# after fixing netsize at 300, build additional 900 ADAGE models at model
# size 300 for ensemble model construction
# Again, this step is very time-consuming. We highly recommend configuring the
# run_ADAGE_train.py so that the construction of each individual ADAGE model is
# run as one job and distributed across to a computing cluster.
model_size=300
python run_ADAGE_train.py $data_compendium $models_folder $model_size 0.1 \
101,1001

# re-do pathway enrichment analysis for each model at model size 300 to
# include the 900 additional models
# This evaluation could run overnight.
Rscript netsize_pathway_enrichment.R $models_folder 300 $pathway_file \
$data_compendium FALSE netsize $N_cores $HW_cutoff

# calculate the reconstruction cost for each model in the model size analysis
python run_ADAGE_cost.py $models_folder $data_compendium \
$models_folder/netsize_cost.txt $compendiumSize

# plot the reconstruction cost at each model size
Rscript netsize_cost_plot.R $models_folder/netsize_cost.txt $plot_folder

###########################################
# evaluate the influence of sample size
###########################################
# build 100 ADAGE models at sample size 100,200,500,800 for model sizes of
# 50 and 300 using 10 different initialization random seeds and 10 different
# reshuffle random seeds.
# This step is very time-consuming. We highly recommend configuring the
# following script so that the construction of each individual ADAGE model is
# run as one job and distributed across to a computing cluster.
python run_ADAGE_train.py $data_compendium $subsample_models_folder \
50,300 0.1 1,11 --reshuffle 1,11 --subsample 100,200,500,800

# pathway enrichment analysis for each model at each sample size
# This evaluation could run overnight.
Rscript netsize_pathway_enrichment.R $subsample_models_folder 50,300 \
$pathway_file $data_compendium FALSE subsample $N_cores $HW_cutoff

# make plots that summarize the pathway coverage in the subsample analysis
Rscript samplesize_evaluation_plots.R $subsample_models_folder $models_folder \
50,300 $compendiumSize $plot_folder

# evaluate the properties of weight vectors at each sample size
Rscript netsize_weight_property.R $models_folder 50,300 $plot_folder \
$compendiumSize $N_genes TRUE $subsample_models_folder

# evaluate the reconstruction cost at each sample size
# testing cost is calculated both on 200 randomly chosen samples unused in the
# training or the rest unused samples
test_size=200

python run_ADAGE_cost.py $subsample_models_folder $data_compendium \
$subsample_models_folder/subsample_train_test_cost.txt $compendiumSize \
--testsize $test_size

Rscript subsample_cost_plot.R \
$subsample_models_folder/subsample_train_test_cost.txt $test_size $plot_folder

#############################################
# evaluate the influence of corruption level
#############################################
# build 100 300-node ADAGE models at corruption level of
# 0.0, 0.05, 0.1, 0.15, ..., 0.95
# This step is very time-consuming. We highly recommend configuring the
# run_ADAGE_train.py so that the construction of each individual ADAGE model is
# run as one job and distributed across a computing cluster.
model_size=300
corrupt_levels=$(seq -s, 0 0.05 1 | sed 's/,$//')
python run_ADAGE_train.py $data_compendium $corruption_models_folder \
$model_size $corrupt_levels 1,101

# pathway enrichment analysis for each model at each corruption level
Rscript netsize_pathway_enrichment.R $corruption_models_folder $model_size \
$pathway_file $data_compendium TRUE corruption $N_cores $HW_cutoff

# make plots that summarize the pathway coverage and pathway significance in
# the corruption analysis
Rscript corruption_evaluation_plots.R \
$corruption_models_folder/netsize${model_size}_\
$(basename -s .txt $pathway_file)_allPathway.txt $plot_folder

# evaluate the reconstruction cost of each model at each corruption level
python run_ADAGE_cost.py $corruption_models_folder $data_compendium \
$corruption_models_folder/corruption_cost.txt $compendiumSize

Rscript corruption_cost_plot.R

# evaluate the relationship between signature overlap and corruption level
Rscript signature_overlap_analysis.R $data_compendium $corruption_models_folder \
300 $corruption_models_folder/netsize$(model_size)_\
$(basename -s .txt $pathway_file)_sigPathway.txt $plot_folder
