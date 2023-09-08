RaceID_clustering<-function(dataset=c('RNAMix1','RNAMix2','TMLung','Beta2',
                                   'TMPanc','BaronPanc','PBMC4k','CellMixSng')){
  
  #Required Libraies ----------------------------------------------------------
  
  require(RaceID)

  #Input dataset ---------------------------------------------------------------
  s1<-proc.time()
  
  if( dataset =='RNAMix1'){
    mydata <- as.matrix(readRDS('../Count_data/rnamix1_filtered.rds'))
    true_labels <- read.csv(file = '../Data_after_Imputation/rnamix1_original_labels.csv')  
    true_labels<-true_labels[,-1]}
 
  
  if( dataset =='RNAMix2'){
    mydata <- as.matrix(readRDS('../Count_data/rnamix2_filtered.rds'))
    true_labels <- read.csv(file = '../Data_after_Imputation/rnamix2_original_labels.csv')
    true_labels<-true_labels[,-1]}
  
  if( dataset =='TMLung'){
    mydata <- read.csv('../Count_data/lung_counts_updated.csv')
    true_labels <- read.csv('../Count_data/lung_counts_updated_labels.csv')
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]}
  
  if( dataset =='Beta2'){
    mydata <- readRDS('../Count_data/beta2_counts_for_raceid.rds')
    true_labels <- read.csv('../Data_after_Imputation/beta_cell_groups_3_4_10_after_filtering.csv')
    true_labels<-true_labels[,-c(1,2,3)]
    }
  
  if(dataset =='TMPanc'){
    mydata <- read.csv('../Count_data/pancreas_tm_counts_updated.csv')
    true_labels <- read.csv('../Count_data/pancreatic_tm_updated_labels.csv')
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-c(1,2)]}
  
  if(dataset =='BaronPanc'){
    mydata <- readRDS('../Count_data/barons_pancreatic_for_Raceid3.rds')
    true_labels <- read.csv("../Data_after_Imputation/pancreatic_num_labels_filtered.csv")
    true_labels<-true_labels[,-1]}
  
  if(dataset =='PBMC4k'){
    mydata <- read.csv('../Count_data/subset_pbmc4k_only_mt.csv')
    true_labels <- read.csv('../Count_data/pbmc4k_only_mt_5_MAIN_LABELS.csv')
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]}
  
  if(dataset =='CellMixSng'){
    mydata <- read.csv("../Count_data/cellmix_sng.csv")
    true_labels <- read.csv('../Data_after_Imputation/cellmix_sng_labels.csv')
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]}
  

  
# Run Raceid3 -------------------------------------------------------

set.seed(20)
sc <- SCseq(mydata)

# Filtering and Normalization --------------------------------------

if(dataset=='Beta2'){
  sc <- filterdata(sc,mintotal = 10,
                   minexpr = 0,
                   minnumber = 0,
                   LBatch = NULL,
                   knn = 10,
                   CGenes = NULL,
                   FGenes = NULL,
                   ccor = 0.4,
                   bmode = "RaceID",
                   verbose = TRUE)}else{
  sc <- filterdata(sc,mintotal = 1000,
                 minexpr = 0,
                 minnumber = 2,
                 LBatch = NULL,
                 knn = 10,
                 CGenes = NULL,
                 FGenes = NULL,
                 ccor = 0.4,
                 bmode = "RaceID",
                 verbose = TRUE)}

sc <-CCcorrect(sc, mode="pca", dimR=TRUE,logscale =TRUE)
#start clock here
s2<-proc.time()
sc <- compdist(sc)

k=length(unique(true_labels))
sc <- clustexp(sc,cln=k,sat=FALSE)
#sc<-findoutliers(sc)
end<-proc.time()
runtime<-data.frame(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)


eval=clustering_evaluation(sc@cluster$kpart ,t(true_labels))
eval= t(data.frame(unlist(eval)[1:3]))
rownames(eval)<-dataset
accuracy<-cbind(eval,runtime)
return(accuracy)
}