setwd('C:/Users/Andriana/Dropbox/scRNAseq/scPMP_final_repo/R_functions')
source('../R_functions/clustering.R')
source('../R_functions/processing_prior_pm_clustering.R')
source("../R_functions/clustering_evaluation.R")

# Ranmix1 --------------------------------------------------------------------
dataset ='RNAMix1'
  mydata <- read.csv(file = 'Data_after_Imputation/RNAmix1_original.csv')
  true_labels <- read.csv(file = 'Data_after_Imputation/rnamix1_original_labels.csv')
  var_cutoff=0.3
  # epsilon = 2.6
  # epsilon_umap = 0.64
  
  rownames(mydata) = mydata[,1]
  mydata<-mydata[,-1]
  true_labels<-true_labels[,-1]
  set.seed(1)
  s1<-proc.time()
  mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)
  
  s2<-proc.time()
  
  eps_res<-matrix(0,nrow=1,ncol=6)
  for(epsilon in seq(from =1, to=3.5, by= 0.2 ) ){
    
    db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts = 10)
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
    end<-proc.time()
    runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
    
    eval<-clustering_evaluation(db_labels,t(true_labels))
    eval= t(data.frame(unlist(eval)))
    rownames(eval)<-dataset
    eval<-cbind(eval,epsilon,length(unique(true_labels)))
    eps_res<-rbind(eps_res,eval)
  }
  
  eps_res[which(eps_res[,4]==eps_res[,6] ),]

# Ranmix2 --------------------------------------------------------------------
 dataset ='RNAMix2'
  mydata <- read.csv(file = 'Data_after_Imputation/RNAmix2_original.csv')
  true_labels <- read.csv(file = 'Data_after_Imputation/rnamix2_original_labels.csv')
  var_cutoff=0.2
  
  #epsilon_umap = 0.7
  #if(two_dimensions){r=2}else{r= 6}
  rownames(mydata) = mydata[,1]
  mydata<-mydata[,-1]
  true_labels<-true_labels[,-1]
  
  set.seed(1)
  s1<-proc.time()
  mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)
  
  s2<-proc.time()
  
  eps_res<-matrix(0,nrow=1,ncol=6)
  for(epsilon in seq(from =1, to=3.5, by= 0.2 ) ){
  
    db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts = 10)
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
  end<-proc.time()
  runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
  
  eval<-clustering_evaluation(db_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  eval<-cbind(eval,epsilon,length(unique(true_labels)))
  eps_res<-rbind(eps_res,eval)
  }
  
  eps_res[which(eps_res[,4]==eps_res[,6] ),]
# eps=1.6 is fine
  
# TMLUNG --------------------
  dataset ='TMLung'
  mydata <- read.csv(file = 'Data_after_Imputation/lung_tabula_muris_saver.csv')
  true_labels <- read.csv(file = 'Data_after_Imputation/lung_tabula_muris_labels.csv')
  var_cutoff=2
  epsilon = 10
  # epsilon_umap =0.52
  # if(two_dimensions){r=2}else{r=5}
  rownames(mydata) = mydata[,1]
  mydata<-mydata[,-1]
  true_labels<-true_labels[,-1]
  
  
  set.seed(1)
  s1<-proc.time()
  mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)
  
  s2<-proc.time()
  
  eps_res<-matrix(0,nrow=1,ncol=6)
  for(epsilon in seq(from =11, to=15, by= 0.001 ) ){
    
    db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts =10)
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
    end<-proc.time()
    runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
    
    eval<-clustering_evaluation(db_labels,t(true_labels))
    eval= t(data.frame(unlist(eval)))
    rownames(eval)<-dataset
    eval<-cbind(eval,epsilon,length(unique(true_labels)))
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
  
  rownames(mydata) = mydata[,1]
  mydata<-mydata[,-1]
  true_labels<-true_labels[,-1]
  
  set.seed(1)
  s1<-proc.time()
  mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)
  
  s2<-proc.time()
  
  eps_res<-matrix(0,nrow=1,ncol=6)
  for(epsilon in seq(from =2.4, to=3, by= 0.01 ) ){
    
    db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts = 10,method='raw')
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
    end<-proc.time()
    runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
    
    eval<-clustering_evaluation(db_labels,t(true_labels))
    eval= t(data.frame(unlist(eval)))
    rownames(eval)<-dataset
    eval<-cbind(eval,epsilon,length(unique(true_labels)))
    eps_res<-rbind(eps_res,eval)
  }
  
  eps_res[which(eps_res[,4]==eps_res[,6] ),] 
  
  #TM_Panc ---------------------------------------------------------
  dataset ='TMPanc'
    mydata <- read.csv(file = 'Data_after_Imputation/pancreas_tabula_muris_saver.csv')
    true_labels <- read.csv(file = 'Data_after_Imputation/pancreatic_tabula_muris_labels.csv')
    var_cutoff=1
    # epsilon = 6
    # epsilon_umap = 0.46
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
   
    
    set.seed(1)
    s1<-proc.time()
    mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)
    
    s2<-proc.time()
    
    eps_res<-matrix(0,nrow=1,ncol=6)
    for(epsilon in seq(from =6.25, to=6.26, by= 0.001 ) ){
      
      db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts = 10)
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
      end<-proc.time()
      runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
      
      eval<-clustering_evaluation(db_labels,t(true_labels))
      eval= t(data.frame(unlist(eval)))
      rownames(eval)<-dataset
      eval<-cbind(eval,epsilon,length(unique(true_labels)))
      eps_res<-rbind(eps_res,eval)
    }
    
    eps_res[which(eps_res[,4]==eps_res[,6] ),] 
    
# Baron Panc -------------------------------------------------------------------
   dataset ='BaronPanc'
    set.seed(1)
    s1<-proc.time()
    
    mydata <- read.csv(file = 'Data_after_Imputation/BaronPancSCT_filtered.csv')
    true_labels <- read.csv(file = 'Data_after_Imputation/pancreatic_num_labels_filtered.csv')
    var_cutoff=NULL
    # epsilon = 5
    # epsilon_umap = 2.5
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]

    mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff)
      s2<-proc.time()
      
      eps_res<-matrix(0,nrow=1,ncol=6)
      for(epsilon in seq(from =4.5, to=5.5, by= 0.1 ) ){
        
        db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts = 10)
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
        end<-proc.time()
        runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
        
        eval<-clustering_evaluation(db_labels,t(true_labels))
        eval= t(data.frame(unlist(eval)))
        rownames(eval)<-dataset
        eval<-cbind(eval,epsilon,length(unique(true_labels)))
        eps_res<-rbind(eps_res,eval)
      }
      
      eps_res[which(eps_res[,4]==eps_res[,6] ),]  
      
# PBMC4K ---------------------------------------------------------------------
      
      dataset ='PBMC4k'
      set.seed(1)
        mydata <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k.csv')
        true_labels <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k_labels.csv')
        var_cutoff=2
        # epsilon = 5
        # epsilon_umap = 5.7
      
        rownames(mydata) = mydata[,1]
        mydata<-mydata[,-1]
        true_labels<-true_labels[,-1]
       
       
         mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)
        
        eps_res<-matrix(0,nrow=1,ncol=6)
        for(epsilon in seq(from =5, to=10, by= 1) ){
          
          db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts = 10)
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
          eval<-cbind(eval,epsilon,length(unique(true_labels)))
          eps_res<-rbind(eps_res,eval)
        }
        
        eps_res[which(eps_res[,4]==eps_res[,6] ),]  
        
# Cellmix --------------------------------------------------------------------
        dataset ='CellMixSng'
        mydata <- read.csv(file = 'Data_after_Imputation/Cellmix_sng_SCT.csv')
        true_labels <- read.csv('Data_after_Imputation/cellmix_sng_labels.csv')
        var_cutoff=NULL
        # epsilon = 9
        # epsilon_umap = 5
     
        rownames(mydata) = mydata[,1]
        mydata<-mydata[,-1]
        true_labels<-true_labels[,-1]
        print('I loaded CellmixSCT')
        mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff)
          
          
          eps_res<-matrix(0,nrow=1,ncol=6)
          for(epsilon in seq(from =2, to=15, by= 1 ) ){
            
            db<-fpc::dbscan(data=t(mydata),eps=epsilon,MinPts = 10)
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
            eval<-cbind(eval,epsilon,length(unique(true_labels)))
            eps_res<-rbind(eps_res,eval)
          }
          
          eps_res[which(eps_res[,4]==eps_res[,6] ),]  
          