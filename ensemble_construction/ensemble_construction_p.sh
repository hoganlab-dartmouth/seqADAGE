# file path to the data compendium built in the previous step
#data_compendium="../data_collection/all-pseudomonas-gene-normalized.pcl"
data_compendium="../data_files/pao1_aligned_rnaseq_compendium_zp2_TPM_log_01_byall.csv"
# count number of genes
N_genes=$(expr $(wc -l <$data_compendium) - 1)
# the model size that best supports the current P.a expression compendium is 300
model_size=300
# number of standard deviations from the mean to be counted as high-weight
HW_cutoff=2.5
# number of cores to use in parallel in netsize_pathway_enrichment.R,
# please modify this to accommodate your need.
N_cores=1
# folder to save ensemble models
ensemble_folder="./ensemble_models"
# file path to the KEGG pathway file "pseudomonas_KEGG_terms.txt" or the GO
# pathway file "manual_GO_BP_terms.txt"
pathway_file="./pseudomonas_KEGG_terms.txt"

##########################
# build ensemble models
##########################

# consensus clustering 100 300-nodeƒ models using the same 100 individual models,
# repeat 10 times for both corADAGE method and eADAGE method

# for eADAGE method, calculate the weighted Pearson correlation between each
# weight vector pair in the first 100 individual models
# run overnight when using 4 cores
#python weighted_pearson_correlation_parallel.py \
#../netsize_evaluation/models/$model_size/ net${model_size}_weighted_cor.txt \
#1,100 $N_cores $N_genes

#python weighted_pearson_correlation_parallel_edit.py \
#../outputs/e_models/$model_size/ net${model_size}_weighted_cor.txt \
#660,735 $N_cores $N_genes

# Note that the ensemble construction is done on a 64-core machine
# the number of processors to use
processors=6 #60
# the folder to store some large temporary files
scratch_folder='/global/scratch/'
seeds=(1) # 2 3 4 5 6 7 8 9 10
for i in "${seeds[@]}";
do
    #mpiexec -np $processors Rscript ensemble_construction.R $data_compendium \
    #../outputs/e_models/$model_size/ $scratch_folder $model_size \
    #737 738 $i weight $ensemble_folder;
    Rscript ensemble_construction_p.R $data_compendium \
    ../outputs/e_models/$model_size/ $scratch_folder $model_size \
    660 735 $i weighted $ensemble_folder ./net${model_size}_weighted_cor.txt;
done

# build additional 9 ensemble models from different 100 individual models for
# both corADAGE method and eADAGE method
#starts=(101 201 301 401 501 601 701 801 901)
#for i in "${starts[@]}";
#do
    #python weighted_pearson_correlation_parallel.py \
    #../outputs/e_models/$model_size/ \
    #net${model_size}_weighted_cor_seed${i}_$((i+99)).txt $i,$((i+99));
    #mpiexec -np $processors Rscript ensemble_construction.R $data_compendium \
    #../outputs/e_models/$model_size/ $scratch_folder $model_size \
    #$i $((i+99)) 1 weight $ensemble_folder;
    #mpiexec -np $processors Rscript ensemble_construction.R $data_compendium \
    #../outputs/e_models/$model_size/ $scratch_folder $model_size \
    #$i $((i+99)) 1 weighted \
    #./net${model_size}_weighted_cor_seed${i}_$((i+99)).txt $ensemble_folder;
#done

##########################
# evaluate ensemble models
##########################

# perform pathway enrichment analysis for ensemble models
Rscript eADAGE_pathway_enrichment.R $ensemble_folder $pathway_file \
$data_compendium TRUE $N_cores $HW_cutoff

#export data_compendium
# remove pathway crosstalk effects and then perform pathway enrichment
# analysis for individual and ensemble models.
# return data about the # of pathways covered in each model (pathway coverage)
#source ./pathway_counts_no_crosstalk.sh

# compare pathway coverage between ensemble models and individual models
Rscript analyze_ensemble_coverage.R $pathway_file

# compare pathway coverage between ensemble models and individual models after
# removing crosstalk effects
#Rscript analyze_ensemble_coverage_no_crosstalk.R

# compare the significance of each pathway association between ensemble models
# and individual models
#Rscript analyze_ensemble_significance.R $pathway_file

# compare weight properties between eADAGE models and individual models
#Rscript compare_weight_property.R $ensemble_folder \
#../outputs/e_models/300/ $N_genes

###################################################
# build the final eADAGE model from the first 100 ADAGE models with a random
# seed 123, this model is used for medium analysis later and save it to its
# own folder
#mpiexec -np $processors Rscript ensemble_construction.R $data_compendium \
#../outputs/e_models/$model_size/ $scratch_folder $model_size \
#1 100 123 weighted ./net${model_size}_weighted_cor.txt ./ensemble_ADAGE
