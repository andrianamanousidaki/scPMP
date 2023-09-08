setwd('C:/Users/Andriana/Dropbox/scRNAseq/scPMP_final_repo/')
source('R_functions/clustering.R')
source('R_functions/processing_prior_pm_clustering.R')
source("R_functions/clustering_evaluation.R")

# Ranmix1 --------------------------------------------------------------------
dataset ='RNAMix1'
mydata <- read.csv(file = 'Data_after_Imputation/RNAmix1_original.csv')
true_labels <- read.csv(file = 'Data_after_Imputation/rnamix1_original_labels.csv')
var_cutoff=0.3

# epsilon = 2.6
 epsilon_umap = 0.64

r = 7

rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]


mydata<-Linnorm::Linnorm(mydata)
mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
set.seed(20)
eps_res<-matrix(0,nrow=1,ncol=6)
for(epsilon_umap in seq(from =0.779, to=0.78, by=0.0001 ) ){
  
  
  umap_embedding<-umap::umap(t(mydata),n_components=r)
  
  db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
  db_labels<-db$cluster
  
  out_points<-which(db_labels==0)
  core_points<-which(db_labels!=0)
  
  if(length(out_points>0)){
    new_labels<-rep(0,times=length(db_labels))
    new_labels[core_points]<-db_labels[core_points]
    core_point_labels<-db_labels
    
    knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
    IDX=knn$nn.idx
    
    for(i in 1:length(out_points)){
      new_labels[out_points[i]]<-core_point_labels[IDX[i]]
    }
    db_labels<-new_labels}

  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
}

eps_res[which(eps_res[,4]==eps_res[,6] ),]
eps_res[which(eps_res[,1]==max(eps_res[,1] )),]

#use 0.78
# Ranmix2 --------------------------------------------------------------------
dataset ='RNAMix2'
mydata <- read.csv(file = 'Data_after_Imputation/RNAmix2_original.csv')
true_labels <- read.csv(file = 'Data_after_Imputation/rnamix2_original_labels.csv')
var_cutoff=0.2

epsilon_umap = 0.7
r= 6

rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]

mydata<-Linnorm::Linnorm(mydata)
mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
set.seed(20)
eps_res<-matrix(0,nrow=1,ncol=6)
for(epsilon_umap in seq(from =0.5, to=0.78, by=0.01 ) ){
  
  
  umap_embedding<-umap::umap(t(mydata),n_components=r)
  
  db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
  db_labels<-db$cluster
  
  out_points<-which(db_labels==0)
  core_points<-which(db_labels!=0)
  
  if(length(out_points>0)){
    new_labels<-rep(0,times=length(db_labels))
    new_labels[core_points]<-db_labels[core_points]
    core_point_labels<-db_labels
    
    knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
    IDX=knn$nn.idx
    
    for(i in 1:length(out_points)){
      new_labels[out_points[i]]<-core_point_labels[IDX[i]]
    }
    db_labels<-new_labels}
  
  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
}

eps_res[which(eps_res[,4]==eps_res[,6] ),]
eps_res[which(eps_res[,1]==max(eps_res[,1] )),]
eps_res[which(eps_res[,4]==eps_res[,6] ),]

# TMLUNG --------------------
dataset ='TMLung'
mydata <- read.csv(file = 'Data_after_Imputation/lung_tabula_muris_saver.csv')
true_labels <- read.csv(file = 'Data_after_Imputation/lung_tabula_muris_labels.csv')
var_cutoff=2
epsilon = 10
# epsilon_umap =0.52
r=5
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]


mydata<-Linnorm::Linnorm(mydata)
mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
set.seed(20)
eps_res<-matrix(0,nrow=1,ncol=6)
for(epsilon_umap in seq(from =0.66, to=0.67, by=0.001 ) ){
  
  
  umap_embedding<-umap::umap(t(mydata),n_components=r)
  
  db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
  db_labels<-db$cluster
  
  out_points<-which(db_labels==0)
  core_points<-which(db_labels!=0)
  
  if(length(out_points>0)){
    new_labels<-rep(0,times=length(db_labels))
    new_labels[core_points]<-db_labels[core_points]
    core_point_labels<-db_labels
    
    knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
    IDX=knn$nn.idx
    
    for(i in 1:length(out_points)){
      new_labels[out_points[i]]<-core_point_labels[IDX[i]]
    }
    db_labels<-new_labels}
  
  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
}
eps_res[which(eps_res[,4]==eps_res[,6] ),] 

# Beta2 -------------------------------------------------------------------
dataset ='Beta2'
mydata <- read.csv(file = '../Data_after_Imputation/beta_3_4_10_filtered_saver.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/beta_cell_groups_3_4_10_after_filtering.csv')
true_labels<-true_labels[,-c(1,2)]
var_cutoff=0.3
epsilon =2.42
# epsilon_umap = 0.626
r=3
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]


mydata<-Linnorm::Linnorm(mydata)
mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
set.seed(20)
eps_res<-matrix(0,nrow=1,ncol=6)
for(epsilon_umap in seq(from =0.5, to=1, by=0.01 ) ){
  
  
  umap_embedding<-umap::umap(t(mydata),n_components=r)
  
  db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
  db_labels<-db$cluster
  
  out_points<-which(db_labels==0)
  core_points<-which(db_labels!=0)
  
  if(length(out_points>0)){
    new_labels<-rep(0,times=length(db_labels))
    new_labels[core_points]<-db_labels[core_points]
    core_point_labels<-db_labels
    
    knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
    IDX=knn$nn.idx
    
    for(i in 1:length(out_points)){
      new_labels[out_points[i]]<-core_point_labels[IDX[i]]
    }
    db_labels<-new_labels}
  
  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
}
eps_res[which(eps_res[,4]==eps_res[,6] ),] 

#TM_Panc ---------------------------------------------------------
dataset ='TMPanc'
mydata <- read.csv(file = 'Data_after_Imputation/pancreas_tabula_muris_saver.csv')
true_labels <- read.csv(file = 'Data_after_Imputation/pancreatic_tabula_muris_labels.csv')
var_cutoff=1
# epsilon = 6.237
# epsilon_umap = 0.46
r=5
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]



mydata<-Linnorm::Linnorm(mydata)
mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
set.seed(20)
eps_res<-matrix(0,nrow=1,ncol=6)
for(epsilon_umap in seq(from =0.6, to=1, by=0.01 ) ){
  
  
  umap_embedding<-umap::umap(t(mydata),n_components=r)
  
  db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
  db_labels<-db$cluster
  
  out_points<-which(db_labels==0)
  core_points<-which(db_labels!=0)
  
  if(length(out_points>0)){
    new_labels<-rep(0,times=length(db_labels))
    new_labels[core_points]<-db_labels[core_points]
    core_point_labels<-db_labels
    
    knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
    IDX=knn$nn.idx
    
    for(i in 1:length(out_points)){
      new_labels[out_points[i]]<-core_point_labels[IDX[i]]
    }
    db_labels<-new_labels}
  
  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
}
eps_res[which(eps_res[,4]==eps_res[,6] ),] 

# Baron Panc -------------------------------------------------------------------
dataset ='BaronPanc'
mydata <- read.csv(file = 'Data_after_Imputation/pancreatic_saver_no_normalization_filtered.csv')
true_labels <- read.csv(file = 'Data_after_Imputation/pancreatic_num_labels_filtered.csv')
var_cutoff=NULL
# epsilon = 4.5
# epsilon_umap = 2.5
r=9
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1] 

  
mydata<-Linnorm::Linnorm(mydata)
mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
set.seed(20)
eps_res<-matrix(0,nrow=1,ncol=6)
for(epsilon_umap in seq(from =0.1, to=3, by=0.1 ) ){
  
  
  umap_embedding<-umap::umap(t(mydata),n_components=r)
  
  db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
  db_labels<-db$cluster
  
  out_points<-which(db_labels==0)
  core_points<-which(db_labels!=0)
  
  if(length(out_points>0)){
    new_labels<-rep(0,times=length(db_labels))
    new_labels[core_points]<-db_labels[core_points]
    core_point_labels<-db_labels
    
    knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
    IDX=knn$nn.idx
    
    for(i in 1:length(out_points)){
      new_labels[out_points[i]]<-core_point_labels[IDX[i]]
    }
    db_labels<-new_labels}
  
  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
}
eps_res[which(eps_res[,4]==eps_res[,6] ),] 

# PBMC4K ---------------------------------------------------------------------

dataset ='PBMC4k'
mydata <- read.csv(file = 'Data_after_Imputation/subset_pbmc4k.csv')
true_labels <- read.csv(file = 'Data_after_Imputation/subset_pbmc4k_labels.csv')
var_cutoff=2
# epsilon = 5
# epsilon_umap = 5.7
r=3
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]


mydata<-Linnorm::Linnorm(mydata)
mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
set.seed(20)
eps_res<-matrix(0,nrow=1,ncol=6)
for(epsilon_umap in seq(from =1, to=7, by=1 ) ){
  
  
  umap_embedding<-umap::umap(t(mydata),n_components=r)
  
  db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
  db_labels<-db$cluster
  
  out_points<-which(db_labels==0)
  core_points<-which(db_labels!=0)
  
  if(length(out_points>0)){
    new_labels<-rep(0,times=length(db_labels))
    new_labels[core_points]<-db_labels[core_points]
    core_point_labels<-db_labels
    
    knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
    IDX=knn$nn.idx
    
    for(i in 1:length(out_points)){
      new_labels[out_points[i]]<-core_point_labels[IDX[i]]
    }
    db_labels<-new_labels}
  
  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
}
eps_res[which(eps_res[,4]==eps_res[,6] ),] 

# Cellmix --------------------------------------------------------------------
dataset ='CellMixSng'
  mydata <- read.csv(file = 'Count_data/cellmix_sng.csv')
  true_labels <- read.csv('Data_after_Imputation/cellmix_sng_labels.csv')
  var_cutoff=NULL
  # epsilon = 8
  # epsilon_umap = 5
  r=4
  rownames(mydata) = mydata[,1]
  mydata<-mydata[,-1]
  true_labels<-true_labels[,-1]
  
  
  mydata<-Linnorm::Linnorm(mydata)
  mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)
  set.seed(20)
  eps_res<-matrix(0,nrow=1,ncol=6)
  for(epsilon_umap in seq(from =1, to=7, by=1 ) ){
    
    
    umap_embedding<-umap::umap(t(mydata),n_components=r)
    
    db<-fpc::dbscan(data=umap_embedding$layout,eps=epsilon_umap,MinPts = 10)
    db_labels<-db$cluster
    
    out_points<-which(db_labels==0)
    core_points<-which(db_labels!=0)
    
    if(length(out_points>0)){
      new_labels<-rep(0,times=length(db_labels))
      new_labels[core_points]<-db_labels[core_points]
      core_point_labels<-db_labels
      
      knn=nn2(data= t(mydata[,core_points]),query=t(mydata[,out_points]), k = 1)
      IDX=knn$nn.idx
      
      for(i in 1:length(out_points)){
        new_labels[out_points[i]]<-core_point_labels[IDX[i]]
      }
      db_labels<-new_labels}
    
    
    eval<-clustering_evaluation(db_labels,t(true_labels))
    eval= t(data.frame(unlist(eval)))
    rownames(eval)<-dataset
    eval<-cbind(eval,epsilon_umap,length(unique(true_labels)))
    eps_res<-rbind(eps_res,eval)
  }
  eps_res[which(eps_res[,4]==eps_res[,6] ),] 
