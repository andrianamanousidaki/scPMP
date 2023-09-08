library(Seurat)
library(patchwork)
library(plyr)


# Functions needed -------------------------------------------------------------

source('../R_functions/processing_prior_pm_clustering.R')
source('../R_functions/Path_Metrics_Clustering_R.R')
source("../R_functions/clustering_evaluation.R")

#The pancreatic data set -------------------------------------------------------

original_data<-read.csv('../Data_after_Imputation/pancreatic_saver_no_normalization_filtered.csv')
rownames(original_data)<-original_data[,1]
original_data<-original_data[,-1]

num_label<-read.csv("../Data_after_Imputation/pancreatic_num_labels_filtered.csv")
num_label<-num_label[,-1]

labels<-revalue(as.character(num_label), replace=c('5'='4',
                                                   '6'='5',
                                                   '8'='6',
                                                   '9'='7',
                                                   '10'='8'))
labels<-as.numeric(labels)
levels<-unique(labels)

write.csv(labels,'pancreatic_num_labels_filtered_ordered.csv')
original_indices<-list()

for(i in levels){
  original_indices[[i]]<-which(labels == i)
}

# Startified Downsampling -----------------------------------------------------------------

down_per<-c(0.5,0.25,0.1)

sampled_indices<- vector(mode = "list", length = length(down_per))
names(sampled_indices)<-down_per
j=0
set.seed(22)
for(per in down_per){
   samp=list()
   j=j+1
   for (i in 1:8) {
     s= round(per*length(original_indices[[i]]))
     
     samp[[i]]<-sample(original_indices[[i]], s, replace = FALSE, prob = NULL)
   }
   sampled_indices[[j]]=samp
}


pancreatic_50 = original_data[, unlist(sampled_indices$'0.5')]
pancreatic_50_labels = labels[unlist(sampled_indices$'0.5')]
write.csv(pancreatic_50, 'pancreatic_50_strat.csv')
write.csv(pancreatic_50_labels, 'pancreatic_50_strat_labels.csv')


pancreatic_25 = original_data[, unlist(sampled_indices$'0.25')]
pancreatic_25_labels = labels[unlist(sampled_indices$'0.25')]
write.csv(pancreatic_25, 'pancreatic_25_strat.csv')
write.csv(pancreatic_25_labels, 'pancreatic_25_strat_labels.csv')


pancreatic_10 = original_data[, unlist(sampled_indices$'0.1')]
pancreatic_10_labels = labels[unlist(sampled_indices$'0.1')]
write.csv(pancreatic_10, 'pancreatic_10_strat.csv')
write.csv(pancreatic_10_labels, 'pancreatic_10_strat_labels.csv')

# Seurat with PM processing clustering ----------------------------------------
set.seed(1)

pancreatic_50_seurat<-processing_prior_pm_clustering(dataset=pancreatic_50,
                                       LogNormalization=TRUE,
                                       var_cutoff=1)
pancreatic_25_seurat<-processing_prior_pm_clustering(dataset=pancreatic_25,
                                                     LogNormalization=TRUE,
                                                     var_cutoff=1)
pancreatic_10_seurat<-processing_prior_pm_clustering(dataset=pancreatic_10,
                                                     LogNormalization=TRUE,
                                                     var_cutoff=1)

data_list <- list(pancreatic_50_seurat, pancreatic_25_seurat, pancreatic_10_seurat)

labels_list <- list(pancreatic_50_labels, pancreatic_25_labels, pancreatic_10_labels)

q = 0

seurat_list <- list()

for(mydata in data_list){
 
   q = q + 1
  
  rownames(mydata) = 1:nrow(mydata)
  colnames(mydata) = paste0("Cell_", 1:ncol(mydata)) 
  
  Seurat_mydata <- CreateSeuratObject(counts = mydata, project = "pancreatic",
                                      min.cells = 0, 
                                      min.features = 0)
  
  Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(mydata))
  Seurat_mydata@assays$RNA@scale.data <- as.matrix(mydata)
  
  rownames(mydata) <- paste0("PC_", 1:min(nrow(mydata),40))
  Seurat_mydata[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(t(mydata)),
                                                 key = "PC_",
                                                 assay='RNA')
  
  Seurat_mydata <- FindNeighbors(Seurat_mydata, dims = 1:min(nrow(mydata),40))
  
  true_labels <- labels_list[[q]] 
  colnames(true_labels)=NULL
  Seurat_mydata<-AddMetaData(Seurat_mydata,true_labels ,
                             col.name = 'true_labels')
  
  seurat_list <- append(seurat_list, Seurat_mydata)
}

results <- matrix(0, nrow = 1, ncol = 5)
colnames(results) <- c('dataset', 'ARI', 'ECP', 'ECA', 'number_of_clusters')

# Pancreatic 50

seurat_list[[1]] <- FindClusters(seurat_list[[1]], resolution = 0.015)
Seurat_labels_50 <- seurat_list[[1]]$RNA_snn_res.0.015
eval<-clustering_evaluation(Seurat_labels_50,labels_list[[1]])
results <-rbind(results,c('seurat_pan50',unlist(eval)))


# Pancreatic 25
seurat_list[[2]] <- FindClusters(seurat_list[[2]], resolution = 0.1)
Seurat_labels_25 <- seurat_list[[2]]$RNA_snn_res.0.1
eval<-clustering_evaluation(Seurat_labels_25,labels_list[[2]])
results <-rbind(results,c('seurat_pan25',unlist(eval)))


# Pancreatic 10
seurat_list[[3]] <- FindClusters(seurat_list[[3]], resolution = 3)
Seurat_labels_10 <- seurat_list[[3]]$RNA_snn_res.3
eval<-clustering_evaluation(Seurat_labels_10,labels_list[[3]])
results <-rbind(results,c('seurat_pan10',unlist(eval)))

saveRDS(results, '../Results/downsampling_res.rds')

# Seurat default clustering ----------------------------------------------------

pan_50_def_seurat <- CreateSeuratObject(counts = pancreatic_50, 
                                               project = "pancreatic_50",
                                               min.cells = 0, 
                                               min.features = 0)

pan_25_def_seurat <- CreateSeuratObject(counts = pancreatic_25, 
                                               project = "pancreatic_25",
                                               min.cells = 0, 
                                               min.features = 0)

pan_10_def_seurat <- CreateSeuratObject(counts = pancreatic_10, 
                                               project = "pancreatic_10",
                                               min.cells = 0, 
                                               min.features = 0)


data_list <- list(pan_50_def_seurat, pan_25_def_seurat, pan_10_def_seurat)

q = 0
n_hvg=2000
dim=50

for(Seurat_mydata in data_list){
  
  q = q + 1
  
  Seurat_mydata[["percent.mt"]] <- PercentageFeatureSet(Seurat_mydata, 
                                                        pattern = "^MT-")
  Seurat_mydata <- NormalizeData(Seurat_mydata)
  Seurat_mydata <- FindVariableFeatures(Seurat_mydata, selection.method = "vst",
                                        nfeatures = n_hvg)
  Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(Seurat_mydata))
  Seurat_mydata<- RunPCA(Seurat_mydata, 
                         features = VariableFeatures(object = Seurat_mydata))
  Seurat_mydata <- JackStraw(Seurat_mydata, num.replicate = 100, dims = dim)
  Seurat_mydata <- ScoreJackStraw(Seurat_mydata, dims = 1:dim)
  
  
  pcid<-which(Seurat_mydata@reductions$pca@jackstraw@overall.p.values[,2]<0.001)
  
  Seurat_mydata <- FindNeighbors(Seurat_mydata, dims =pcid )
  
  true_labels <- labels_list[[q]]
  Seurat_mydata<-AddMetaData(Seurat_mydata, true_labels, col.name = 'true_labels')
  
  data_list[[q]] <- Seurat_mydata
}

  
# Pancreatic 50

data_list[[1]] <- FindClusters(data_list[[1]], resolution = 0.1)
Seurat_def_labels_50 <- data_list[[1]]$RNA_snn_res.0.1
eval <- clustering_evaluation(Seurat_def_labels_50, labels_list[[1]])
results <- rbind(results,c('seurat_def_pan50', unlist(eval)))

# Pancreatic 25

data_list[[2]] <- FindClusters(data_list[[2]], resolution = 0.43)
Seurat_def_labels_25 <- data_list[[2]]$RNA_snn_res.0.43
eval <- clustering_evaluation(Seurat_def_labels_25, labels_list[[2]])
results <- rbind(results,c('seurat_def_pan25', unlist(eval)))

# Pancreatic 10

data_list[[3]] <- FindClusters(data_list[[3]], resolution = 2)
Seurat_def_labels_10 <- data_list[[3]]$RNA_snn_res.2
eval <- clustering_evaluation(Seurat_def_labels_10, labels_list[[3]])
results <- rbind(results,c('seurat_def_pan10', unlist(eval)))

saveRDS(results, '../Results/downsampling_res.rds')

# Path metrics clustering -----------------------------------------------------

set.seed(1)

pancreatic_50_pm<-processing_prior_pm_clustering(dataset=pancreatic_50,
                                                     LogNormalization=TRUE,
                                                     var_cutoff=1)
pancreatic_25_pm<-processing_prior_pm_clustering(dataset=pancreatic_25,
                                                     LogNormalization=TRUE,
                                                     var_cutoff=1)
pancreatic_10_pm<-processing_prior_pm_clustering(dataset=pancreatic_10,
                                                     LogNormalization=TRUE,
                                                     var_cutoff=NULL)

# PM Clustering with p=1.5 

pan_50_pm_1.5 <- Path_Metrics_Clustering_R(dataset = pancreatic_50_pm,
                                     n_clust = 8,
                                     p = 1.5,
                                     Fast_mds = FALSE,
                                     silhouette = FALSE)

eval <- clustering_evaluation(pan_50_pm_1.5$pm_labels, labels_list[[1]])
results <- rbind(results,c('pan_50_pm_1.5', unlist(eval)))

pan_25_pm_1.5 <- Path_Metrics_Clustering_R(dataset = pancreatic_25_pm,
                                           n_clust = 8,
                                           p = 1.5,
                                           Fast_mds = FALSE,
                                           silhouette = FALSE)

eval <- clustering_evaluation(pan_25_pm_1.5$pm_labels, labels_list[[2]])
results <- rbind(results,c('pan_25_pm_1.5', unlist(eval)))

pan_10_pm_1.5 <- Path_Metrics_Clustering_R(dataset = pancreatic_10_pm,
                                           n_clust = 8,
                                           p = 1.5,
                                           Fast_mds = FALSE,
                                           silhouette = FALSE)

eval <- clustering_evaluation(pan_10_pm_1.5$pm_labels, labels_list[[3]])
results <- rbind(results,c('pan_10_pm_1.5', unlist(eval)))

# PM Clustering with p=2 

pan_50_pm_2 <- Path_Metrics_Clustering_R(dataset = pancreatic_50_pm,
                                           n_clust = 8,
                                           p = 2,
                                           Fast_mds = FALSE,
                                           silhouette = FALSE)

eval <- clustering_evaluation(pan_50_pm_2$pm_labels, labels_list[[1]])
results <- rbind(results,c('pan_50_pm_2', unlist(eval)))

pan_25_pm_2 <- Path_Metrics_Clustering_R(dataset = pancreatic_25_pm,
                                           n_clust = 8,
                                           p = 2,
                                           Fast_mds = FALSE,
                                           silhouette = FALSE)

eval <- clustering_evaluation(pan_25_pm_2$pm_labels, labels_list[[2]])
results <- rbind(results,c('pan_25_pm_2', unlist(eval)))

pan_10_pm_2 <- Path_Metrics_Clustering_R(dataset = pancreatic_10_pm,
                                           n_clust = 8,
                                           p = 2,
                                           run.adjusted.kmeans=FALSE,
                                           Fast_mds = FALSE,
                                           silhouette = FALSE)

eval <- clustering_evaluation(pan_10_pm_2$pm_labels, labels_list[[3]])
results <- rbind(results,c('pan_10_pm_2', unlist(eval)))
# num_cl=9 for pan_10_pm_2 two tiny clusters and sdjusted kmeans can't solve it

# PM Clustering with p=4 

pan_50_pm_4 <- Path_Metrics_Clustering_R(dataset = pancreatic_50_pm,
                                         n_clust = 8,
                                         p = 4,
                                         Fast_mds = FALSE,
                                         silhouette = FALSE)

eval <- clustering_evaluation(pan_50_pm_4$pm_labels, labels_list[[1]])
results <- rbind(results,c('pan_50_pm_4', unlist(eval)))

pan_25_pm_4 <- Path_Metrics_Clustering_R(dataset = pancreatic_25_pm,
                                         n_clust = 8,
                                         p = 4,
                                         Fast_mds = FALSE,
                                         silhouette = FALSE)

eval <- clustering_evaluation(pan_25_pm_4$pm_labels, labels_list[[2]])
results <- rbind(results,c('pan_25_pm_4', unlist(eval)))

pan_10_pm_4 <- Path_Metrics_Clustering_R(dataset = pancreatic_10_pm,
                                         n_clust = 8,
                                         p = 4,
                                         Fast_mds = FALSE,
                                         silhouette = FALSE)

eval <- clustering_evaluation(pan_10_pm_4$pm_labels, labels_list[[3]])
results <- rbind(results,c('pan_10_pm_4', unlist(eval)))

saveRDS(results, '../Results/downsampling_res.rds')
library(writexl)
write_xlsx(as.data.frame(results), '../Results/downsampling_res.xlsx')
