Seurat_def_clustering<-function(dataset=c('RNAMix1','RNAMix2','TMLung','Beta2',
                                'TMPanc','BaronPanc','PBMC4k','CellMixSng')){
  
#Required Libraies ----------------------------------------------------------
  s1<-proc.time()
require(Seurat)

#Input dataset ---------------------------------------------------------------


if( dataset =='RNAMix1'){
  mydata <- read.csv(file = '../Data_after_Imputation/RNAmix1_original.csv')
  true_labels <- read.csv(file = '../Data_after_Imputation/rnamix1_original_labels.csv')}


if( dataset =='RNAMix2'){
  mydata <- read.csv(file = '../Data_after_Imputation/RNAmix2_original.csv')
  true_labels <- read.csv(file = '../Data_after_Imputation/rnamix2_original_labels.csv')}


if( dataset =='TMLung'){
  mydata <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_saver.csv')
  true_labels <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_labels.csv')}

if( dataset =='Beta2'){
  mydata <- read.csv(file = '../Data_after_Imputation/beta_3_4_10_filtered_saver.csv')
  true_labels <- read.csv(file = '../Data_after_Imputation/beta_cell_groups_3_4_10_after_filtering.csv')
  true_labels<-true_labels[,-c(1,2)]}

if(dataset =='TMPanc'){
  mydata <- read.csv(file = '../Data_after_Imputation/pancreas_tabula_muris_saver.csv')
  true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_tabula_muris_labels.csv')}

if(dataset =='BaronPanc'){
  mydata <- read.csv(file = '../Data_after_Imputation/BaronPancSCT_filtered.csv')
  true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_num_labels_filtered.csv')}

if(dataset =='PBMC4k'){
  mydata <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k.csv')
  true_labels <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k_labels.csv')}

if(dataset =='CellMixSng'){
  mydata <- read.csv(file = '../Data_after_Imputation/Cellmix_sng_SCT.csv')
  true_labels <- read.csv('../Data_after_Imputation/cellmix_sng_labels.csv')}

# Prep data set ----------------------------------------------------------------
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]


#  Create Seurat Object --------------------------------------------------------
set.seed(20)
Seurat_mydata <- CreateSeuratObject(counts = mydata, project = dataset, min.cells = 0, min.features = 0)
n_hvg=2000
dim=50

# Normalization and Dim reduction of Seurat Object ------------------------------

if(dataset=='TMPanc' | dataset=='TMLung'){ 
  Seurat_mydata[["percent.mt"]] <- PercentageFeatureSet(Seurat_mydata, pattern = "^mt-")
}else 
{Seurat_mydata[["percent.mt"]] <-PercentageFeatureSet(Seurat_mydata, pattern = "^MT-")}

if(dataset!='BaronPanc' &&  dataset!='CellMixSng'){ 
  Seurat_mydata <- NormalizeData(Seurat_mydata)
  Seurat_mydata <- FindVariableFeatures(Seurat_mydata, selection.method = "vst", nfeatures = n_hvg)
  Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(Seurat_mydata))
  Seurat_mydata<- RunPCA(Seurat_mydata, features = VariableFeatures(object = Seurat_mydata))
  Seurat_mydata <- JackStraw(Seurat_mydata, num.replicate = 100, dims = dim)
  Seurat_mydata <- ScoreJackStraw(Seurat_mydata, dims = 1:dim)
}else{
  Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(mydata))
  Seurat_mydata@assays$RNA@scale.data <- as.matrix(mydata)
  Seurat_mydata<- RunPCA(Seurat_mydata, features = rownames(Seurat_mydata))
  Seurat_mydata <- JackStraw(Seurat_mydata, num.replicate = 100, dims = dim)
  Seurat_mydata <- ScoreJackStraw(Seurat_mydata, dims = 1:dim)
}

# find important PCs ------------------------------------------------------------
pcid<-which(Seurat_mydata@reductions$pca@jackstraw@overall.p.values[,2]<0.001)
#find neighbors
s2<-proc.time()
Seurat_mydata <- FindNeighbors(Seurat_mydata, dims =pcid )

Seurat_mydata<-AddMetaData(Seurat_mydata, true_labels, col.name = 'true_labels')

# Clustering ---------------------------------------------------------------------


if( dataset =='RNAMix1'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.5)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.5}

if( dataset =='RNAMix2'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.5)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.5}

if( dataset =='TMLung'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.5)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.5}


if( dataset =='Beta2'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.3)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.3}

if(dataset =='TMPanc'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.3)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.3}

if(dataset =='BaronPanc'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.12)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.12}


if(dataset =='PBMC4k'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.007)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.007}


if(dataset =='CellMixSng'){
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.01)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.01}

end<-proc.time()
runtime<-data.frame(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)

# Evaluation ------------------------------------------------------
eval<-clustering_evaluation(Seurat_labels,t(true_labels))
eval= t(data.frame(unlist(eval)))
rownames(eval)<-dataset
accuracy<-cbind(eval,runtime)
return(accuracy)
}