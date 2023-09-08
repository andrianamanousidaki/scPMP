SIMLR_clustering<-function(dataset=c('RNAMix1','RNAMix2','TMLung','Beta2',
                                   'TMPanc','BaronPanc','PBMC4k','CellMixSng')){
  
  #Required Libraies ----------------------------------------------------------
  s1<-proc.time()
  require(SIMLR)

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
  
  
  # Normalize Data    ----------------------------------------------------------------
  if( dataset!='BaronPanc' && dataset!='CellMixSng'){ 
  mydata<-as.matrix(mydata)%*%diag(1/colSums(as.matrix(mydata)))*10000
  mydata<-log(mydata+1)}
  
  # Run SIMLR -------------------------------------------------------------------------
  
  set.seed(175)
  s2<-proc.time()
  clustering<-SIMLR(X=mydata, c=length(unique(true_labels)),normalize =FALSE)
  end<-proc.time()
  runtime<-data.frame(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
  
  
  # Evaluate and save evaluation ------------------------------------------------------
  eval<-clustering_evaluation(clustering$y$cluster,t(true_labels))
  eval= t(data.frame(unlist(eval)[1:3]))
  rownames(eval)<-dataset
  
  accuracy<-cbind(eval,runtime)
  return(accuracy)
  
}