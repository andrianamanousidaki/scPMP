SC3_clustering<-function(dataset=c('RNAMix1','RNAMix2','TMLung','Beta2',
                                   'TMPanc','BaronPanc','PBMC4k','CellMixSng')){
  
  #Required Libraies ----------------------------------------------------------
   s1<-proc.time()
  require(SingleCellExperiment)
  require(SC3)
  require(scater)
  
  
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
  
  # Prep data set ------------------------------------------------------------
  rownames(mydata) = mydata[,1]
  mydata<-mydata[,-1]
  true_labels<-true_labels[,-1]
  
  
  #Clustering -----------------------------------------------------------------
  
  # Define the number of clusters 
  n_clusters=length(unique(true_labels))
  
  if( dataset!='BaronPanc' && dataset!='CellMixSng'){ 
    
  set.seed(190)
 
  # Normalization 
  mydata<-as.matrix(mydata)%*%diag(1/colSums(as.matrix(mydata)))*10000
  
  
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(mydata),
                                            logcounts = log2(as.matrix(mydata) + 1)))
  
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  sce <- sc3_prepare(sce,kmeans_nstart =50,kmeans_iter_max = 100, gene_filter=FALSE )
  sce <- sc3_calc_dists(sce)
  sce <- sc3_calc_transfs(sce)
  s2<-proc.time()
  t= (3+n_clusters)
  t
  sce <- sc3_kmeans(sce, ks = 2:t)
  sce <- sc3_calc_consens(sce)
  end<-proc.time()
  runtime<-data.frame(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
  
  
  # Evaluation 
  eval=clustering_evaluation(colData(sce)[,(n_clusters-1)] ,t(true_labels))
  eval= t(data.frame(unlist(eval)[1:3]))
  rownames(eval)<-dataset
  accuracy<-cbind(eval,runtime)
  return(accuracy)}else{
    
    set.seed(190)
    
    sce <- SingleCellExperiment(assays = list(counts = as.matrix(mydata),
                    logcounts = as.matrix(mydata)))
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
    sce <- sc3_prepare(sce,kmeans_nstart =50,kmeans_iter_max = 100, gene_filter=FALSE )
    sce <- sc3_calc_dists(sce)
    sce <- sc3_calc_transfs(sce)
    t= (3+n_clusters)
    t
    s2<-proc.time()
    sce <- sc3_kmeans(sce, ks = 2:t)
    sce <- sc3_calc_consens(sce)
    end<-proc.time()
    runtime<-data.frame(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
    
    
    # Evaluation 
    eval=clustering_evaluation(colData(sce)[,(n_clusters-1)] ,t(true_labels))
    eval= t(data.frame(unlist(eval)[1:3]))
    rownames(eval)<-dataset
    accuracy<-cbind(eval,runtime)
    return(accuracy)
  }
  
}