source("../R_functions/clustering_evaluation.R")

### Step 1: Choose data set from Original Datasets:

# RNAMix1
dataset='RNAMix1'
predictions<-read.csv('../Results/SCANPY_RNAMix1Basic_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/rnamix1_original_labels.csv')

##RNAMix2
dataset='RNAMix2'

predictions <- read.csv(file ='../Results/SCANPY_RNAMix2Basic_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/rnamix2_original_labels.csv')

## TMLung
dataset='TMLung'

predictions <- read.csv(file = '../Results/SCANPY_TMLungBasic_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_labels.csv')

## Beta2
dataset='Beta2'

predictions <- read.csv(file = '../Results//SCANPY_Beta2Basic_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/beta_cell_groups_3_4_10_after_filtering.csv')
true_labels<-true_labels[,-c(1,2)]

## TMPanc
dataset='TMPanc'

predictions<- read.csv(file = '../Results/SCANPY_TMPancBasic_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_tabula_muris_labels.csv')

## BaronPancSCT
dataset='BaronPancSCT'

predictions <- read.csv(file = '../Results/SCANPY_BaronPancSCT_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_num_labels_filtered.csv')


## PBMC4kBASIC
dataset='PBMC4k'

predictions <- read.csv(file = '../Results/SCANPY_PBMC4kBASIC_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k_labels.csv')


## CellMixSngSCT
dataset='CellMixSngSCT'
predictions <- read.csv(file = '../Results/SCANPY_CellMixSngSCT_noDR_40PCs_for_nn_RESULTS.csv')
true_labels <- read.csv('../Data_after_Imputation/cellmix_sng_labels.csv')

###########################################################################\]
true_labels<-true_labels[,-1]
predictions<-predictions[,-1]

eval<-clustering_evaluation(t(predictions) ,t(true_labels))
eval= t(data.frame(unlist(eval)[1:3]))
rownames(eval)=dataset

# results=data.frame(ARI=c(0),ECP=c(0),ECA=c(0))
results= readRDS('../Results/scanpy_ecp_eca_ari_RESULTS.rds')

results= rbind(results,eval)

saveRDS(results,'../Results/scanpy_ecp_eca_ari_RESULTS.rds')

##############################################################################

# runtime<-read.csv('../Results/SCANPY_all_runtime.csv')
# runtime<-runtime[,-1]
# 
# 
# 
# results$complete_rt<-runtime$Complete_rt
# results$clustering_rt<-runtime$Clustering_rt
# 
# results$dataset<-rownames(results)
# library(writexl)
# write_xlsx(results,paste0('../Results/reproduce_results_scanpy.xlsx'))
# #saveRDS(results,'../Results/reproduce_results_scanpy.rds')

