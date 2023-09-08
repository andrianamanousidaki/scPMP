# This function reproduces clustering results 

clustering<-function(dataset=c('RNAMix1','RNAMix2','TMLung','Beta2','TMPanc','BaronPanc',
                                 'PBMC4k','CellMixSng','Balls','Enlongated_With_Bridge',
                               'Swiss','SO(3)'),
                     method=c('RaceID_clustering','SC3_clustering','Seurat_clustering_default_res_0.8',
                              'Seurat_Clustering_PM_processing','Seurat_def_clustering',
                              'SIMLR_clustering','Path_metrics_clustering','kmeans',
                              'dbscan',	'umap+dbscan',	'tsne+kmeans'),
                     p=NULL,
                     geometric_perturbation=FALSE,
                     two_dimensions=FALSE,
                     predict_pm_clusters=FALSE){
  
  
  #SOURCE CLUSTERING FUNCTION
 
  source('../R_functions/clustering_evaluation.R')
  source('../R_functions/processing_prior_pm_clustering.R')
  source('../R_functions/geometric_perturbation.R')
  require('RANN')
  require(ClusterR)
  
  adjusted_kmeans<-function(U,k){
    #this function needs library(pracma)
    km<-KMeans_rcpp(U,clusters=k,num_init = 20)
    pm_labels<-km$clusters
    H=table(pm_labels)
    
    if(min(H)< sqrt(nrow(U)/2)){
      num_tiny_clusters=length(which(H<sqrt(nrow(U)/2)))
      km2<- KMeans_rcpp(U,clusters=k+num_tiny_clusters,num_init = 20)
      pm_labels_original<-pm_labels
      pm_labels<-km2$clusters
      C<-km2$centroids
      H<-table(pm_labels)
      h2<-sort.int(H,index.return=TRUE)
      H2<-h2$x
      idxH<-h2$ix
      tiny_cluster_idx = idxH[1:num_tiny_clusters]
      CD =squareform(as.vector(dist(C, method = "euclidean")))# P CENTROID DISTANCE MATRIX 
      options = CD[,tiny_cluster_idx]
      
      if(length(tiny_cluster_idx)==1){
        v=options
        merge_idx = which.min(v[v>0])
        pm_labels[pm_labels==tiny_cluster_idx] = merge_idx
      }else{
        for (i in 1:length(tiny_cluster_idx)){
          v=options[,i]
          merge_idx = which( options[,i] == min(v[v>0]))
          pm_labels[pm_labels==tiny_cluster_idx[i]] = merge_idx}
      }
      
      
      remaining_labels = unique(pm_labels)
      remaining_labels = remaining_labels[order(remaining_labels)]
      for (i in 1:length(remaining_labels)){
        pm_labels[pm_labels==remaining_labels[i]]=i
      }
      
    }
    
    klabels=pm_labels
    return(klabels)
  }
  
  
  # Clustering --------------------------------------------------------------
  if(method=='RaceID_clustering'){
    
    source('../R_functions/RaceID_clustering.R')
    
    accuracy<-RaceID_clustering(dataset=dataset)
    accuracy<-cbind(accuracy,method)
    accuracy<-cbind(accuracy,dataset)
  }
  if(method=='SC3_clustering'){
   
    source('../R_functions/SC3_clustering.R')
    
    accuracy<-SC3_clustering(dataset=dataset)
    accuracy<-cbind(accuracy,method)
    accuracy<-cbind(accuracy,dataset)
  }
  
  if(method=='Seurat_clustering_default_res_0.8'){
   
    source('../R_functions/Seurat_clustering_default_res_0.8.R')
    
    accuracy<-Seurat_clustering_default_res_0.8(dataset=dataset)
    accuracy<-list(accuracy=accuracy,dataset=dataset,method=method)
  }
  
  if(method=='Seurat_Clustering_PM_processing'){
    
    source('../R_functions/Seurat_Clustering_PM_processing.R')
    
    accuracy<-Seurat_clustering_PM_processing(dataset=dataset)
    accuracy<-cbind(accuracy,method)
    accuracy<-cbind(accuracy,dataset)}
  
  if(method=='Seurat_def_clustering'){
  
    source('../R_functions/Seurat_def_clustering.R')
    
    accuracy<-Seurat_def_clustering(dataset=dataset)
    accuracy<-cbind(accuracy,method)
    accuracy<-cbind(accuracy,dataset)
  }
  
  if(method=='SIMLR_clustering'){
   
    source('../R_functions/SIMLR_clustering.R')
    
    accuracy<-SIMLR_clustering(dataset=dataset)
    accuracy<-cbind(accuracy,method)
    accuracy<-cbind(accuracy,dataset)
  }
  
  if(method=='Path_metrics_clustering' ||method=='kmeans'||method=='dbscan'||method=='umap+dbscan' ||method=='tsne+kmeans'){
    
    if( dataset =='RNAMix1'){
      mydata <- read.csv(file = '../Data_after_Imputation/RNAmix1_original.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/rnamix1_original_labels.csv')
      var_cutoff=0.3
      epsilon = 2.6
      epsilon_umap = 0.78
      if(two_dimensions){r=2}else{r = 7}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      fast=FALSE
      }
    
    
    if( dataset =='RNAMix2'){
      mydata <- read.csv(file = '../Data_after_Imputation/RNAmix2_original.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/rnamix2_original_labels.csv')
      var_cutoff=0.2
      epsilon = 1.6
      epsilon_umap = 0.64
      if(two_dimensions){r=2}else{r= 6}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      fast=FALSE
      }
    
    
    if( dataset =='TMLung'){
      mydata <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_saver.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_labels.csv')
      var_cutoff=2
      epsilon = 8.4
      epsilon_umap =0.66
      if(two_dimensions){r=2}else{r=5}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      fast=FALSE
      }
    
    if( dataset =='Beta2'){
      mydata <- read.csv(file = '../Data_after_Imputation/beta_3_4_10_filtered_saver.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/beta_cell_groups_3_4_10_after_filtering.csv')
      true_labels<-true_labels[,-c(1,2)]
      var_cutoff=0.5
      epsilon = 2.6
      epsilon_umap = 0.89
      if(two_dimensions){r=2}else{r=3}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      fast=FALSE
      }
    
    if(dataset =='TMPanc'){
      mydata <- read.csv(file = '../Data_after_Imputation/pancreas_tabula_muris_saver.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_tabula_muris_labels.csv')
      var_cutoff=1
      epsilon = 6.237
      epsilon_umap = 0.69
      if(two_dimensions){r=2}else{r=5}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      fast=FALSE
      }
    
    if(dataset =='BaronPanc' && method !='umap+dbscan' && method!='tsne+kmeans'){
      mydata <- read.csv(file = '../Data_after_Imputation/BaronPancSCT_filtered.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_num_labels_filtered.csv')
      var_cutoff=NULL
      epsilon = 4.5
      epsilon_umap = 0.7
      if(two_dimensions){r=2}else{r=9}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      fast=FALSE
    }else{
      if(dataset =='BaronPanc'){
      mydata <- read.csv(file = '../Data_after_Imputation/pancreatic_saver_no_normalization_filtered.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_num_labels_filtered.csv')
      var_cutoff=NULL
      epsilon = 4.5
      epsilon_umap = 0.7
      if(two_dimensions){r=2}else{r=9}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1] }
      
      }
    
    if(dataset =='PBMC4k'){
      mydata <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k.csv')
      true_labels <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k_labels.csv')
      var_cutoff=2
      epsilon = 5
      epsilon_umap = 4
      if(two_dimensions){r=2}else{r=3}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      fast=TRUE
      }
    
    if(dataset =='CellMixSng'&& method !='umap+dbscan' && method!='tsne+kmeans'){
      mydata <- read.csv(file = '../Data_after_Imputation/Cellmix_sng_SCT.csv')
      true_labels <- read.csv('../Data_after_Imputation/cellmix_sng_labels.csv')
      var_cutoff=NULL
      epsilon = 8
      epsilon_umap = 6
      if(two_dimensions){r=2}else{r=4}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      print('I loaded CellmixSCT')
      fast=TRUE
    }else{
      if(dataset =='CellMixSng'){
      mydata <- read.csv(file = '../Count_data/cellmix_sng.csv')
      true_labels <- read.csv('../Data_after_Imputation/cellmix_sng_labels.csv')
      var_cutoff=NULL
      epsilon = 8
      epsilon_umap = 6
      if(two_dimensions){r=2}else{r=4}
      rownames(mydata) = mydata[,1]
      mydata<-mydata[,-1]
      true_labels<-true_labels[,-1]
      }
      }
    
    # Manifold datasets
    
    if(dataset =='SO(3)'){
      mydata <- read.csv(file = '../Count_data/GLmanifold_d9_N3000_k3_sig0075.csv',header = FALSE)
      true_labels <- read.csv('../Count_data/GLmanifold_d9_N3000_k3_sig0075_Labels.csv',header = FALSE)
      mydata <-t(mydata)
      var_cutoff=NULL
      epsilon = 0.15
      epsilon_umap = 3
      if(two_dimensions){r=2}else{r=7}
      fast=FALSE}
    
    if(dataset =='Balls'){
      mydata <- read.csv(file = '../Count_data/Balls_n400_r53.csv',header = FALSE)
      true_labels <- read.csv('../Count_data/Balls_n400_r53_Labels.csv',header = FALSE)
      mydata <-t(mydata)
      var_cutoff=NULL
      epsilon = 0.08687
      epsilon_umap = 0.64
      if(two_dimensions){r=2}else{r=5}
      fast=FALSE}
    
    if(dataset =='Enlongated_With_Bridge'){
      mydata <- read.csv(file = '../Count_data/ElongatedGaussiansWithBridge3.csv',header = FALSE)
      true_labels <- read.csv('../Count_data/ElongatedGaussiansWithBridge3_Labels.csv',header = FALSE)
      mydata <-t(mydata)
      var_cutoff=NULL
      epsilon = 0.2
      epsilon_umap = 0.8
      if(two_dimensions){r=2}else{r=5}
      fast=FALSE}
    
    if(dataset =='Swiss'){
      mydata <- read.csv(file = '../Count_data/SwissRoll1.csv',header = FALSE)
      true_labels <- read.csv('../Count_data/SwissRoll1_Labels.csv',header = FALSE)
      mydata <-t(mydata)
      var_cutoff=NULL
      epsilon = 3
      epsilon_umap = 5
      if(two_dimensions){r=2}else{r=5}
      fast=FALSE}
    
    
    
    
    if(method=='Path_metrics_clustering'){
      source('../R_functions/Path_Metrics_Clustering_R.R')
      
      # Process data set : 
      set.seed(1)
      s1<-proc.time()
      if(dataset=='BaronPanc' || dataset=='CellMixSng'){ 
        
        mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff)}else{
          
          if(dataset!='Balls'&& dataset!='Enlongated_With_Bridge'&&dataset!='Swiss' && dataset!='SO(3)'){
            mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)}else{
              mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,EliminateGenes =FALSE,var_cutoff=var_cutoff)}
          
        }

      s2<-proc.time()
      pm_output<-Path_Metrics_Clustering_R(dataset=mydata,n_clust=length(unique(unlist(true_labels))),p=p,Fast_mds=fast,silhouette=predict_pm_clusters)
      end<-proc.time()
      
      eval<-clustering_evaluation(pm_output$pm_labels,t(true_labels))
      eval= t(data.frame(unlist(eval)[1:3]))
      rownames(eval)<-dataset
      accuracy<-eval
      
      runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
      acc<-list(accuracy=cbind(accuracy,method))
      accuracy<-append(pm_output,acc)
      accuracy<-append(accuracy,runtime)
      
      if(geometric_perturbation){
        if(two_dimensions){
          Emb=pm_output$U1
          gp<-geometric_perturbation(X=mydata,U=Emb[1:2,],Labels=true_labels)
        }else{
        gp<-geometric_perturbation(X=mydata,U=pm_output$U1,Labels=true_labels)}
        accuracy<-append(accuracy,gp)
      }
      }
    
    if(method=='kmeans'){
      set.seed(1)
      s1<-proc.time()
      # Process data set : 
      if(dataset=='BaronPanc' || dataset=='CellMixSng'){ 
        
        mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff)}else{
          
          if(dataset!='Balls'&& dataset!='Enlongated_With_Bridge'&&dataset!='Swiss' && dataset!='SO(3)'){
            mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)}else{
              mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,EliminateGenes =FALSE,var_cutoff=var_cutoff)}
          
        }
      s2<-proc.time()
      km<-KMeans_rcpp(t(mydata),clusters=length(unique(unlist(true_labels))),num_init = 20)
      end<-proc.time()
      runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
      km_labels<-km$clusters
      eval<-clustering_evaluation(km_labels,t(true_labels))
      eval= t(data.frame(unlist(eval)[1:3]))
      rownames(eval)<-dataset
      accuracy<-eval
      accuracy<-list(accuracy=accuracy,method=method,labels=km_labels,runtime=runtime)
   
    }
    
    if(method=='dbscan'){
      set.seed(1)
      s1<-proc.time()
      # Process data set : 
      if(dataset=='BaronPanc' || dataset=='CellMixSng'){ 
        
        mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff)}else{
          
          if(dataset!='Balls'&& dataset!='Enlongated_With_Bridge'&&dataset!='Swiss' && dataset!='SO(3)'){
            mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)}else{
              mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,EliminateGenes =FALSE,var_cutoff=var_cutoff)}
          
        }
      s2<-proc.time()
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
      eval= t(data.frame(unlist(eval)[1:3]))
      rownames(eval)<-dataset
      accuracy<-eval
      accuracy<-list(accuracy=accuracy,method=method,labels=db_labels,runtime=runtime)
    }
    
    
    if(method=='umap+dbscan'){
      set.seed(20)
      s1<-proc.time()
      if(dataset!='Balls'&& dataset!='Enlongated_With_Bridge'&&dataset!='Swiss' && dataset!='SO(3)'){
      mydata<-Linnorm::Linnorm(mydata)
      mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)}
      
      umap_embedding<-umap::umap(t(mydata),n_components=r)
      
      s2<-proc.time()
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
      end<-proc.time()
      runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
      
      eval<-clustering_evaluation(db_labels,t(true_labels))
      eval= t(data.frame(unlist(eval)[1:3]))
      rownames(eval)<-dataset
      accuracy<-eval
      umap_out<-list(umap_em=umap_embedding$layout,labels=db_labels)
      acc<-list(accuracy=cbind(accuracy,method))
      accuracy<-append(umap_out,acc)
      accuracy<-append(accuracy,runtime)
      if(geometric_perturbation){
          gp<-geometric_perturbation(X=mydata,U=t(umap_embedding$layout),Labels=true_labels)
          accuracy<-append(accuracy,gp)
      }
    }
    
    if(method=='tsne+kmeans'){
      set.seed(20)
      s1<-proc.time()
      if(dataset!='Balls'&& dataset!='Enlongated_With_Bridge'&&dataset!='Swiss' && dataset!='SO(3)'){
        mydata<-Linnorm::Linnorm(mydata)
        mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff,LocalAvgAllPts = FALSE)}
      if(!two_dimensions){
      tsne_embedding<-Rtsne::Rtsne(t(mydata),dims=3)}else{
        tsne_embedding<-Rtsne::Rtsne(t(mydata),dims=r)
      }
      s2<-proc.time()
      akmeans_labels<-adjusted_kmeans(tsne_embedding$Y,k=length(unique(unlist(true_labels))))
      end<-proc.time()
      runtime<-list(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
      
      eval<-clustering_evaluation(akmeans_labels,t(true_labels))
      eval= t(data.frame(unlist(eval)[1:3]))
      rownames(eval)<-dataset
      accuracy<-eval
      tsne_out<-list(tsne_em=tsne_embedding$Y,labels=akmeans_labels)
      acc<-list(accuracy=cbind(accuracy,method))
      accuracy<-append(tsne_out,acc)
      accuracy<-append(accuracy,runtime)
      if(geometric_perturbation){
        gp<-geometric_perturbation(X=mydata,U=t(tsne_embedding$Y),Labels=true_labels)
        accuracy<-append(accuracy,gp)
      }
    }
  
  }
  
  
  
  
  
 
  return(accuracy)
  
}